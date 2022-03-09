'''
@author: niko.popitsch@imba.oeaw.ac.at, tobias.neumann@imp.ac.at
'''
from collections import *
import csv, datetime, time, logging, sys, os, json, pickle
import pysam
from genomic_iterators import Location, get_chrom_dicts_from_bam, get_chrom_dicts
from genomic_iterators.bam_iterators import ReadIterator
from genomic_iterators.grouped_iterators import BlockedIterator, BlockLocationIterator, AnnotationOverlapIterator
from splice_sim.utils import *
from splice_sim.model import *
from splice_sim.iterator import *
import glob

def basename(f):
    return os.path.splitext(os.path.basename(f))[0]

def calc_coverage(true_chr, read_chr, true_pos, read_pos):
    """ Calculate the coverage (number of bases) of a read and its simulation origin """
    if true_chr != read_chr:
        return 0
    cov = len(true_pos & read_pos) / len(true_pos)
    return cov
def blocks_overlap(a,b):
    return a[0] <= b[1] and b[0] <= a[1]
def get_overlapping_isoforms(read, tid, m):
    """ gets all ids of isoforms of the passed transcript that overlap the passed read """
    t=m.transcripts[tid]
    ret=[]
    for iso_id, iso in t.isoforms.items():
        overlaps=False
        for ablock in iso.aln_blocks:
            for rblock in read.get_blocks():
                if blocks_overlap(ablock, rblock):
                    overlaps=True
        if overlaps:
            ret+=[iso_id]
    return ret

def classify_read(m, read, overlapping_tids, cr, mapper, dict_chr2idx, performance, out_reads, sam_out, overlap_cutoff = 0.8):
    """ Classifies a read """
    read_name = read.query_name
    is_converted_bam=cr>0 # are there any conversions?
    true_tid, true_strand, true_isoform, read_tag, true_chr, true_start, true_cigar, n_seqerr, n_converted, is_converted_read = read.query_name.split('_')
    read_chr = 'NA' if read.is_unmapped else read.reference_name
    read_coord = 'NA' if read.is_unmapped else '%s:%i-%i' % (read.reference_name, read.reference_start, read.reference_end)

    read_pos=set([b for a,b in read.get_aligned_pairs() if a is not None and b is not None])

    # recreate true read
    true_read = pysam.AlignedSegment()
    true_read.cigarstring=true_cigar
    true_read.reference_start=int(true_start)
    true_pos = set([b for a,b in true_read.get_aligned_pairs() if a is not None and b is not None])

    overlap = calc_coverage(true_chr, read_chr, true_pos, read_pos) #overlap with transcript region
    true_coord = "%s:%i-%i" % ( true_chr, min(true_pos),  max(true_pos))

    # read strongly overlaps with real location; FIXME: check strand!
    if overlap >= overlap_cutoff:
        performance[true_tid, true_isoform, 'TP'] += 1
    else:
        performance[true_tid, true_isoform, 'FN'] += 1
        if out_reads is not None:
            print('\t'.join([str(x) for x in (read_coord, true_coord, 'FN', 'NA',
                                              1 if is_converted_bam else 0, mapper,
                                              cr,
                                              overlap, true_tid,
                                              true_strand, true_isoform, read_tag, true_chr,n_seqerr, n_converted,
                                              is_converted_read,
                                              ','.join(overlapping_tids) if len(overlapping_tids) > 0 else 'NA',
                                              read_name)]), file=out_reads)
        # ignore case where overlap with true transcript is too small.
        overlapping_tids=[t for t in overlapping_tids if t != true_tid]
        # count for each isoform
        for tid in overlapping_tids:
            # check overlap with all isoforms:
            iso_ids=get_overlapping_isoforms(read, tid, m)
            for iso_id in iso_ids:
                performance[tid, iso_id, 'FP'] += 1/len(overlapping_tids)
                if out_reads is not None:
                    print('\t'.join([str(x) for x in (read_coord, true_coord, 'FP', tid,
                                                      1 if is_converted_bam else 0, mapper,
                                                      cr,
                                                      overlap, iso_id,
                                                      true_strand, true_isoform, read_tag, true_chr,n_seqerr, n_converted,
                                                      is_converted_read,
                                                      ','.join(overlapping_tids) if len(overlapping_tids) > 0 else 'NA',
                                                      read_name)]), file=out_reads)
        # write FN/FP reads
        if sam_out is not None:
            true_read.query_name = read_name
            is_correct_strand = ( read.is_reverse and true_strand=='-' ) or ((not read.is_reverse) and true_strand=='+')
            true_read.query_sequence = read.query_sequence if is_correct_strand else reverse_complement(read.query_sequence)
            true_read.query_qualities = read.query_qualities
            true_read.reference_id = dict_chr2idx[true_chr]
            true_read.tags = (("NM", n_seqerr + n_converted ),)
            true_read.flag = 0 if true_strand=='+' else 16
            true_read.set_tag(tag='YC', value='255,255,255', value_type="Z")
            sam_out.write(true_read)
            read.set_tag(tag='YC', value='255,0,0', value_type="Z")
            sam_out.write(read)
    return performance

def evaluate_bam(bam_file, bam_out, m, mapper, cr, out_performance, out_reads):
    """ evaluate the passed BAM file """
    logging.info("Evaluating %s" % bam_file)
    performance = Counter()
    dict_chr2idx, dict_idx2chr, dict_chr2len=get_chrom_dicts_from_bam(bam_file)
    n_reads=0
    is_converted_bam=cr>0
    samin = pysam.AlignmentFile(bam_file, "rb")
    samout = pysam.AlignmentFile(bam_out+'.tmp.bam', "wb", template=samin ) if bam_out is not None else None
    chrom2trans=OrderedDict()
    for c, c_len in dict_chr2len.items():
        chrom2trans[c]=[]
    for tid,t in m.transcripts.items():
        chrom2trans[t.chromosome].append((Location(dict_chr2idx[t.chromosome], t.start, t.end, strand=t.strand),tid))
    for c, c_len in dict_chr2len.items():
        if c not in samin.references:
            continue # no data on this chrom
        if len(chrom2trans[c])==0: # no transcript on this chrom: all reads are FN
            rit = ReadIterator(bam_file, dict_chr2idx, reference=c, start=1, end=c_len, max_span=None, flag_filter=0) # max_span=m.max_ilen
            for loc, read in rit:
                n_reads+=1
                performance = classify_read(m, read, [], cr, mapper, dict_chr2idx, performance, out_reads, samout)
        else:
            aits = [BlockLocationIterator(iter(chrom2trans[c]))]
            rit = ReadIterator(bam_file, dict_chr2idx, reference=c, start=1, end=c_len, max_span=None, flag_filter=0) # max_span=m.max_ilen
            it = AnnotationOverlapIterator(rit, aits, check_read_alignment=False)
            for loc, (read, annos) in it:
                n_reads+=1
                overlapping_tids = [t[0] for (_, (_, t)) in annos[0]]
                performance = classify_read(m, read, overlapping_tids, cr, mapper, dict_chr2idx, performance, out_reads, samout)
    # add unmapped reads
    for read in samin.fetch(contig=None, until_eof=True):
        if not read.is_unmapped:
            continue
        n_reads+=1
        performance = classify_read(m, read, [], cr, mapper, dict_chr2idx, performance, out_reads, samout)
    # print("%s reads:  %i %i %i" % (bam_file, n_reads, samin.mapped, samin.unmapped))
    # write tid performance table
    tids= set([x for x, y, z in performance.keys()])
    isos= set([y for x, y, z in performance.keys()])
    for tid in tids:
        for iso in isos:
             print("\t".join([str(x) for x in [
                 1 if is_converted_bam else 0,
                 mapper,
                 cr,
                 iso,
                 tid,
                 performance[tid, iso, 'TP'],
                 performance[tid, iso, 'FP'],
                 performance[tid, iso, 'FN'] ]]), file=out_performance)

    samin.close()
    if samout is not None:
        samout.close()

def evaluate_bam_performance(config, m, bam_file, out_dir):
    # parse mapper and cr from bam file name
    assert os.path.exists(bam_file) and os.path.exists(bam_file+'.bai'), "Could not find bam file (or idx) for %s" % (bam_file)
    assert '.cr' in bam_file, "Cannot parse bam file name %s" % (bam_file)
    tmp = bam_file[:-10] if bam_file.endswith('.final.bam') else bam_file[:-4] # file ends with .bam or .final.bam
    cr,mapper=tmp[tmp.find('.cr')+3:].rsplit('.', 1)
    out_file_prefix=out_dir+'/'+config['dataset_name']+'.cr'+cr+'.'+mapper
    cr=float(cr)
    is_truth=mapper=='truth'
    # write read data?
    write_reads=config['write_reads'] if 'write_reads' in config else False
    if write_reads:
        out_reads=open(out_file_prefix+'.reads.tsv', 'w')
        print("mapped_coords\ttrue_coords\tclassification\ttid\tis_converted_bam\tmapper\tcondition_id\toverlap\ttrue_tid\ttrue_strand\ttrue_isoform\ttag\ttrue_chr\tn_true_seqerr\tn_tc_pos\tis_converted_read\toverlapping_tids\tread_name", file=out_reads)
    else:
        out_reads=None
    logging.info("write_reads: %s" % str(write_reads))
    # evaluate
    with open(out_file_prefix+'.tid_performance.tsv', 'w') as out_performance:
        print("is_converted_bam\tmapper\tcondition_id\tiso\ttid\tTP\tFP\tFN", file=out_performance)
        if is_truth:
            # do not write unmapped bam or reads.tsv
            evaluate_bam(bam_file, None, m, mapper, cr, out_performance, None)
        else:
            bam_out_file=out_file_prefix+".mismapped.bam"
            evaluate_bam(bam_file, bam_out_file, m, mapper, cr, out_performance, out_reads)
    if write_reads:
        out_reads.close()


def evaluate_splice_sites_performance(config, m, bam_file, out_dir):
    """ evaluate the passed BAM file """
    logging.info("Evaluating splice sites in %s" % bam_file)
    # write intron bam?
    write_intron_bam=config['write_intron_bam'] if 'write_intron_bam' in config else False
    # parse mapper and cr from bam file name
    assert os.path.exists(bam_file) and os.path.exists(bam_file+'.bai'), "Could not find bam file (or idx) for %s" % (bam_file)
    assert '.cr' in bam_file, "Cannot parse bam file name %s" % (bam_file)
    tmp = bam_file[:-10] if bam_file.endswith('.final.bam') else bam_file[:-4] # file ends with .bam or .final.bam
    cr,mapper=tmp[tmp.find('.cr')+3:].rsplit('.', 1)
    out_file_prefix=out_dir+'/'+config['dataset_name']+'.cr'+cr+'.'+mapper
    is_converted_bam=float(cr)>0 # are there any conversions?
    is_truth=mapper=='truth'
    bam_out_file=out_file_prefix+".intron.bam"
    # build transcript intervaltree
    tiv = m.build_transcript_iv()
    # evaluate
    performance = Counter()
    with open(out_file_prefix+'.splice_site_performance.tsv', 'w') as out_performance:
        print("is_converted_bam\tmapper\tcondition\ttid\tintron_id\tdonor_rc_splicing_TP\tdonor_rc_splicing_FP\tdonor_rc_overlapping_TP\tdonor_rc_overlapping_FP\tdonor_rc_splicing_wrong\tacceptor_rc_splicing_TP\tacceptor_rc_splicing_FP\tacceptor_rc_overlapping_TP\tacceptor_rc_overlapping_FP\tacceptor_rc_splicing_wrong", file=out_performance)
        n_reads = 0
        samin = pysam.AlignmentFile(bam_file, "rb")
        samout = pysam.AlignmentFile(bam_out_file+'.tmp.bam', "wb", template=samin) if write_intron_bam else None
        performance = Counter()
        for tid, t in m.transcripts.items():
            for intron in t.introns:
                intronID = tid + "_" + str(intron['exon_number'])
                iv = Interval(intron['start'] - 1, intron['end'] + 1)
                # get tids of overlapping features
                annots=[x.data.tid for x in tiv[t.chromosome].overlap(iv.begin, iv.end)]
                # Donor
                donorIterator = SimulatedReadIterator(samin.fetch(contig=t.chromosome, start=intron['start'] - 1, stop=intron['start']), flagFilter=0)
                for read in donorIterator:
                    start = read.start
                    end = read.end
                    bamRead = read.bamRead
                    color = ""
                    if start <= iv.begin:
                        if read.splicing and not read.hasSpliceSite(iv) and end > iv.begin:
                            # mapped read contains a splice-site but it doesn't cover the intron.
                            print(read.splicing, iv, read.spliceSites)
                            performance[tid, intronID, 'donor_rc_splicing_wrong'] += 1
                            bamRead.setTag('YC', colors["exonexonfalse"])
                            if (samout):
                                samout.write(bamRead)
                        elif read.splicing and read.hasSpliceSite(iv) and end > iv.begin:
                            if read.trueTid == tid and read.trueSpliced:
                                performance[tid, intronID, 'donor_rc_splicing_TP'] += 1
                                bamRead.setTag('YC', colors["exonexontrue"])
                                if (samout):
                                    samout.write(bamRead)
                            elif not read.trueTid in annots or not read.trueSpliced:
                                performance[tid, intronID, 'donor_rc_splicing_FP'] += 1
                                bamRead.setTag('YC', colors["exonexonfalse"])
                                if (samout):
                                    samout.write(bamRead)
                        elif end > iv.begin:
                            if read.trueTid == tid and not read.trueSpliced:
                                performance[tid, intronID, 'donor_rc_overlapping_TP'] += 1
                                bamRead.setTag('YC', colors["exonintron"])
                                if (samout):
                                    samout.write(bamRead)
                            elif not read.trueTid in annots or read.trueSpliced:
                                performance[tid, intronID, 'donor_rc_overlapping_FP'] += 1
                                bamRead.setTag('YC', colors["exonexonfalse"])
                                if (samout):
                                    samout.write(bamRead)
                        else:
                            bamRead.setTag('YC', colors["intron"])
                            if (samout):
                                samout.write(bamRead)
                    n_reads += 1
                # Acceptor
                acceptorIterator = SimulatedReadIterator(samin.fetch(contig=t.chromosome, start=intron['end'], stop=intron['end'] + 1), flagFilter=0)
                for read in acceptorIterator:
                    start = read.start
                    end = read.end
                    bamRead = read.bamRead
                    if end >= iv.end:
                        if read.splicing and not read.hasSpliceSite(iv) and start < iv.end: # simulated (true) read is splicing, mapped read is not splicing
                            performance[tid, intronID, 'acceptor_rc_splicing_wrong'] += 1
                            bamRead.setTag('YC', colors["exonexonfalse"])
                            if (samout):
                                samout.write(bamRead)
                        elif read.splicing and read.hasSpliceSite(iv) and start < iv.end:
                            if read.trueTid == tid and read.trueSpliced:
                                performance[tid, intronID, 'acceptor_rc_splicing_TP'] += 1
                                bamRead.setTag('YC', colors["exonexontrue"])
                                if (samout):
                                    samout.write(bamRead)
                            elif not read.trueTid in annots or not read.trueSpliced:
                                performance[tid, intronID, 'acceptor_rc_splicing_FP'] += 1
                                bamRead.setTag('YC', colors["exonexonfalse"])
                                if (samout):
                                    samout.write(bamRead)
                        elif start < iv.end:
                            if read.trueTid == tid and not read.trueSpliced:
                                performance[tid, intronID, 'acceptor_rc_overlapping_TP'] += 1
                                bamRead.setTag('YC', colors["exonintron"])
                                if (samout):
                                    samout.write(bamRead)
                            elif not read.trueTid in annots or read.trueSpliced:
                                performance[tid, intronID, 'acceptor_rc_overlapping_FP'] += 1
                                bamRead.setTag('YC', colors["exonexonfalse"])
                                if (samout):
                                    samout.write(bamRead)
                        else:
                            # Start ends at first exon position of acceptor
                            pass
                            # bamRead.setTag('YC', colors["intron"])
                            # if (samout):
                            #     samout.write(bamRead)
                    n_reads += 1
        print("%s reads:  %i %i %i" % (bam_file, n_reads, samin.mapped, samin.unmapped))
        samin.close()
        prevIntron = ""
        prevTid = ""
        intronBuffer = Counter()
        for tid, intron, classification in performance.keys():
            if intron != prevIntron and prevIntron != "":
                print("\t".join([str(x) for x in [
                    1 if is_converted_bam else 0,
                    mapper,
                    cr,
                    prevTid,
                    prevIntron,
                    intronBuffer['donor_rc_splicing_TP'],
                    intronBuffer['donor_rc_splicing_FP'],
                    intronBuffer['donor_rc_overlapping_TP'],
                    intronBuffer['donor_rc_overlapping_FP'],
                    intronBuffer['donor_rc_splicing_wrong'],
                    intronBuffer['acceptor_rc_splicing_TP'],
                    intronBuffer['acceptor_rc_splicing_FP'],
                    intronBuffer['acceptor_rc_overlapping_TP'],
                    intronBuffer['acceptor_rc_overlapping_FP'],
                    intronBuffer['acceptor_rc_splicing_wrong']
                ]]), file=out_performance)
                intronBuffer = Counter()
            intronBuffer[classification] = performance[tid, intron, classification]
            prevTid = tid
            prevIntron = intron
        print("\t".join([str(x) for x in [
            1 if is_converted_bam else 0,
            mapper,
            cr,
            prevTid,
            prevIntron,
            intronBuffer['donor_rc_splicing_TP'],
            intronBuffer['donor_rc_splicing_FP'],
            intronBuffer['donor_rc_overlapping_TP'],
            intronBuffer['donor_rc_overlapping_FP'],
            intronBuffer['donor_rc_splicing_wrong'],
            intronBuffer['acceptor_rc_splicing_TP'],
            intronBuffer['acceptor_rc_splicing_FP'],
            intronBuffer['acceptor_rc_overlapping_TP'],
            intronBuffer['acceptor_rc_overlapping_FP'],
            intronBuffer['acceptor_rc_splicing_wrong'],
        ]]), file=out_performance)
        if samout is not None:
            samout.close()


def calc_splicesite_performance_OLD(config, m, bam_file, out_dir):
    """ Count and classify reads overlapping splice_sites """
    logging.info("Evaluating splice sites in %s" % bam_file)
    # write intron bam?
    write_intron_bam=config['write_intron_bam'] if 'write_intron_bam' in config else False
    # parse mapper and cr from bam file name
    assert os.path.exists(bam_file) and os.path.exists(bam_file+'.bai'), "Could not find bam file (or idx) for %s" % (bam_file)
    assert '.cr' in bam_file, "Cannot parse bam file name %s" % (bam_file)
    tmp = bam_file[:-10] if bam_file.endswith('.final.bam') else bam_file[:-4] # file ends with .bam or .final.bam
    cr,mapper=tmp[tmp.find('.cr')+3:].rsplit('.', 1)
    out_file_prefix=out_dir+'/'+config['dataset_name']+'.cr'+cr+'.'+mapper
    is_converted_bam=float(cr)>0 # are there any conversions?
    is_truth=mapper=='truth'
    bam_out_file=out_file_prefix+".intron.bam"
    # build transcript intervaltree
    tiv = m.build_transcript_iv()
    # evaluate
    performance = Counter()
    with open(out_file_prefix+'.splice_site_performance.tsv', 'w') as out_performance:
        print("is_converted_bam\tmapper\tcondition\ttid\tintron_id\tdonor_rc_splicing_TP\tdonor_rc_splicing_FP\tdonor_rc_overlapping_TP\tdonor_rc_overlapping_FP\tdonor_rc_splicing_wrong\tacceptor_rc_splicing_TP\tacceptor_rc_splicing_FP\tacceptor_rc_overlapping_TP\tacceptor_rc_overlapping_FP\tacceptor_rc_splicing_wrong", file=out_performance)
        n_reads = 0
        samin = pysam.AlignmentFile(bam_file, "rb")
        samout = pysam.AlignmentFile(bam_out_file+'.tmp.bam', "wb", template=samin) if write_intron_bam else None
        performance = Counter()
        for tid, t in m.transcripts.items():
            for intron in t.introns:
                intronID = tid + "_" + str(intron['exon_number'])
                iv = Interval(intron['start'] - 1, intron['end'] + 1)
                # get tids of overlapping features
                annots=[x.data.tid for x in tiv[t.chromosome].overlap(iv.begin, iv.end)]
                # Iterate reads overlapping Donor read
                donorIterator = SimulatedReadIterator(samin.fetch(contig=t.chromosome, start=intron['start'] - 1, stop=intron['start']), flagFilter=0)
                for read in donorIterator:
                    start = read.start
                    end = read.end
                    bamRead = read.bamRead
                    color = ""
                    if start <= iv.begin:
                        if read.splicing and not read.hasSpliceSite(iv) and end > iv.begin:
                            # mapped read contains a splice-site but it doesn't cover the intron.
                            print(read.splicing, iv, read.spliceSites)
                            performance[tid, intronID, 'donor_rc_splicing_wrong'] += 1
                            bamRead.setTag('YC', colors["exonexonfalse"])
                            if (samout):
                                samout.write(bamRead)
                        elif read.splicing and read.hasSpliceSite(iv) and end > iv.begin:
                            if read.trueTid == tid and read.trueSpliced:
                                performance[tid, intronID, 'donor_rc_splicing_TP'] += 1
                                bamRead.setTag('YC', colors["exonexontrue"])
                                if (samout):
                                    samout.write(bamRead)
                            elif not read.trueTid in annots or not read.trueSpliced:
                                performance[tid, intronID, 'donor_rc_splicing_FP'] += 1
                                bamRead.setTag('YC', colors["exonexonfalse"])
                                if (samout):
                                    samout.write(bamRead)
                        elif end > iv.begin:
                            if read.trueTid == tid and not read.trueSpliced:
                                performance[tid, intronID, 'donor_rc_overlapping_TP'] += 1
                                bamRead.setTag('YC', colors["exonintron"])
                                if (samout):
                                    samout.write(bamRead)
                            elif not read.trueTid in annots or read.trueSpliced:
                                performance[tid, intronID, 'donor_rc_overlapping_FP'] += 1
                                bamRead.setTag('YC', colors["exonexonfalse"])
                                if (samout):
                                    samout.write(bamRead)
                        else:
                            bamRead.setTag('YC', colors["intron"])
                            if (samout):
                                samout.write(bamRead)
                    n_reads += 1
                # Acceptor
                acceptorIterator = SimulatedReadIterator(samin.fetch(contig=t.chromosome, start=intron['end'], stop=intron['end'] + 1), flagFilter=0)
                for read in acceptorIterator:
                    start = read.start
                    end = read.end
                    bamRead = read.bamRead
                    if end >= iv.end:
                        if read.splicing and not read.hasSpliceSite(iv) and start < iv.end: # simulated (true) read is splicing, mapped read is not splicing
                            performance[tid, intronID, 'acceptor_rc_splicing_wrong'] += 1
                            bamRead.setTag('YC', colors["exonexonfalse"])
                            if (samout):
                                samout.write(bamRead)
                        elif read.splicing and read.hasSpliceSite(iv) and start < iv.end:
                            if read.trueTid == tid and read.trueSpliced:
                                performance[tid, intronID, 'acceptor_rc_splicing_TP'] += 1
                                bamRead.setTag('YC', colors["exonexontrue"])
                                if (samout):
                                    samout.write(bamRead)
                            elif not read.trueTid in annots or not read.trueSpliced:
                                performance[tid, intronID, 'acceptor_rc_splicing_FP'] += 1
                                bamRead.setTag('YC', colors["exonexonfalse"])
                                if (samout):
                                    samout.write(bamRead)
                        elif start < iv.end:
                            if read.trueTid == tid and not read.trueSpliced:
                                performance[tid, intronID, 'acceptor_rc_overlapping_TP'] += 1
                                bamRead.setTag('YC', colors["exonintron"])
                                if (samout):
                                    samout.write(bamRead)
                            elif not read.trueTid in annots or read.trueSpliced:
                                performance[tid, intronID, 'acceptor_rc_overlapping_FP'] += 1
                                bamRead.setTag('YC', colors["exonexonfalse"])
                                if (samout):
                                    samout.write(bamRead)
                        else:
                            # Start ends at first exon position of acceptor
                            pass
                            # bamRead.setTag('YC', colors["intron"])
                            # if (samout):
                            #     samout.write(bamRead)
                    n_reads += 1
        print("%s reads:  %i %i %i" % (bam_file, n_reads, samin.mapped, samin.unmapped))
        samin.close()
        prevIntron = ""
        prevTid = ""
        intronBuffer = Counter()
        for tid, intron, classification in performance.keys():
            if intron != prevIntron and prevIntron != "":
                print("\t".join([str(x) for x in [
                    1 if is_converted_bam else 0,
                    mapper,
                    cr,
                    prevTid,
                    prevIntron,
                    intronBuffer['donor_rc_splicing_TP'],
                    intronBuffer['donor_rc_splicing_FP'],
                    intronBuffer['donor_rc_overlapping_TP'],
                    intronBuffer['donor_rc_overlapping_FP'],
                    intronBuffer['donor_rc_splicing_wrong'],
                    intronBuffer['acceptor_rc_splicing_TP'],
                    intronBuffer['acceptor_rc_splicing_FP'],
                    intronBuffer['acceptor_rc_overlapping_TP'],
                    intronBuffer['acceptor_rc_overlapping_FP'],
                    intronBuffer['acceptor_rc_splicing_wrong']
                ]]), file=out_performance)
                intronBuffer = Counter()
            intronBuffer[classification] = performance[tid, intron, classification]
            prevTid = tid
            prevIntron = intron
        print("\t".join([str(x) for x in [
            1 if is_converted_bam else 0,
            mapper,
            cr,
            prevTid,
            prevIntron,
            intronBuffer['donor_rc_splicing_TP'],
            intronBuffer['donor_rc_splicing_FP'],
            intronBuffer['donor_rc_overlapping_TP'],
            intronBuffer['donor_rc_overlapping_FP'],
            intronBuffer['donor_rc_splicing_wrong'],
            intronBuffer['acceptor_rc_splicing_TP'],
            intronBuffer['acceptor_rc_splicing_FP'],
            intronBuffer['acceptor_rc_overlapping_TP'],
            intronBuffer['acceptor_rc_overlapping_FP'],
            intronBuffer['acceptor_rc_splicing_wrong'],
        ]]), file=out_performance)
        if samout is not None:
            samout.close()

def supports_splice_site(unal, intron_iv, max_diff=5):
        """ Checks overlap of passed intron interval and read unaligned blocks.
            Overlap check is fuzzy to account for some small alignment errors.
            Total diff between endpoints must be <= max_diff
        """
        for ss in unal:
            diff=abs(ss[0]-intron_iv.begin)+abs(ss[1]-intron_iv.end)
            if diff <= max_diff:
                return True
        return False

def aligns_to_iv(alblocks, iv, slop=0):
        """ Checks overlap of passed interval +/- slop and read blocks.""" 
        iv=(iv.begin-slop, iv.end+slop)
        for ss in alblocks:
            if (ss[0] <= iv[1] and iv[0] <= ss[1]):
                return True
        return False
def aligns_to_SJ(alblocks, intron_iv, slop=10):
        """ Checks overlap of passed intron interval +/- slop and read blocks. 
        """ 
        sj1=(intron_iv.begin-slop, intron_iv.begin+slop)
        sj2=(intron_iv.end-slop, intron_iv.end+slop)
        for ss in alblocks:
            if (ss[0] <= sj1[1] and sj1[0] <= ss[1]) or (ss[0] <= sj2[1] and sj2[0] <= ss[1]):
                return True
        return False
    
def set_compare(truth, found):
    """ Compares two sets """
    TP = truth & found
    FP = found.difference(truth)
    FN = truth.difference(found)
    return TP, FP, FN

def evaluate_splice_site_performance(config, m, bam_file, truth_bam_dir, out_dir):
    """ Count and classify reads overlapping splice_sites """
    logging.info("Evaluating splice sites in %s" % bam_file)
    dict_chr2idx, dict_idx2chr, dict_chr2len=get_chrom_dicts_from_bam(bam_file)
    # write intron bam?
    write_intron_bam=config['write_intron_bam'] if 'write_intron_bam' in config else False
    # parse mapper and cr from bam file name
    assert os.path.exists(bam_file) and os.path.exists(bam_file+'.bai'), "Could not find bam file (or idx) for %s" % (bam_file)
    assert '.cr' in bam_file, "Cannot parse bam file name %s" % (bam_file)
    tmp = bam_file[:-10] if bam_file.endswith('.final.bam') else bam_file[:-4] # file ends with .bam or .final.bam
    cr,mapper=tmp[tmp.find('.cr')+3:].rsplit('.', 1)
    out_file_prefix=out_dir+'/'+config['dataset_name']+'.cr'+cr+'.'+mapper
    is_converted_bam=float(cr)>0 # are there any conversions?
    bam_out_file=out_file_prefix+".intron.bam"
    truth_bam_file =  truth_bam_dir+'/'+config['dataset_name']+'.cr'+cr+'.truth.bam'
    truthBam = pysam.AlignmentFile(truth_bam_file, "rb")
    # build transcript intervaltree
    tiv = m.build_transcript_iv()
    # open BAMs
    truth_bam = pysam.AlignmentFile(truth_bam_file, "rb")
    bam = pysam.AlignmentFile(bam_file, "rb")
    samout = pysam.AlignmentFile(bam_out_file+'.tmp.bam', "wb", template=bam) if write_intron_bam else None
    # iterate transcripts/introns
    with open(out_file_prefix+'.splice_site_performance.tsv', 'w') as out_performance:
        print("is_converted_bam\tmapper\tcondition\ttid\tintron_id\tspl_TP\tspl_FP\tspl_FN\tdon_TP\tdon_FP\tdon_FN\tacc_TP\tacc_FP\tacc_FN\tchromosome\tstart\tend\tstrand", file=out_performance)
        for tid, t in m.transcripts.items():
            if 'debug_tid' in config and config['debug_tid'] != tid:
                continue
            for intron in t.introns:
                intronID = tid + "_" + str(intron['exon_number'])
                iv = Interval(intron['start'] - 1, intron['end'] + 1)
                # get all true reads that splice/overlap an intron.
                true_overlapping_don=set()
                true_overlapping_acc=set()
                true_splicing=set()
                rit = ReadIterator(truth_bam, dict_chr2idx, reference=t.chromosome, start=intron['start'] - 1, end=intron['end'], max_span=None, flag_filter=0)
                for _,read in rit:
                    #  skip reads that do not align to [SJ+/- slop] 
                    if not aligns_to_SJ(read.get_blocks(), iv):
                        continue                    
                    true_tid, true_strand, true_isoform, read_tag, true_chr, true_start, true_cigar, n_seqerr, n_converted, is_converted_read = read.query_name.split('_')
                    # NOTE: read may also stem from different (overlapping) tid, but we ignore for this analysis
                    # recreate true read
                    true_read = pysam.AlignedSegment()
                    true_read.cigarstring=true_cigar
                    true_read.reference_start=int(true_start)
                    # check SJ support
                    true_spliced = supports_splice_site(get_unaligned_blocks(true_read), iv)
                    true_overlaps_don = true_read.reference_start <= iv.begin and true_read.reference_end >= iv.begin and true_read.reference_end<=iv.end
                    true_overlaps_acc = true_read.reference_start <= iv.end and true_read.reference_end >= iv.end and true_read.reference_start>iv.begin
                    if true_spliced:
                        true_splicing.add(read.query_name)
                    elif true_overlaps_don:
                        true_overlapping_don.add(read.query_name)
                    elif true_overlaps_acc:
                        true_overlapping_acc.add(read.query_name)
                # calc mapped overlap
                mapped_overlapping_don=set()
                mapped_overlapping_acc=set()
                mapped_splicing=set()
                rit = ReadIterator(bam, dict_chr2idx, reference=t.chromosome, start=intron['start'] - 1, end=intron['end'], max_span=None, flag_filter=0)
                for _, read in rit:
                    # skip reads that do not align to [SJ+/- slop] 
                    if not aligns_to_SJ(read.get_blocks(), iv):
                        continue                    
                    true_tid, true_strand, true_isoform, read_tag, true_chr, true_start, true_cigar, n_seqerr, n_converted, is_converted_read = read.query_name.split('_')
                    spliced = supports_splice_site(get_unaligned_blocks(read), iv)
                    overlaps_don = read.reference_start <= iv.begin and read.reference_end >= iv.begin and read.reference_end<=iv.end
                    overlaps_acc = read.reference_start <= iv.end and read.reference_end >= iv.end and read.reference_start>iv.begin
                    if spliced:
                        mapped_splicing.add(read.query_name)
                    elif overlaps_don:
                        mapped_overlapping_don.add(read.query_name)
                    elif overlaps_acc:
                        mapped_overlapping_acc.add(read.query_name)
                # if transcript is on negative strand: switch donor/acceptor sets
                if t.strand=='-':
                    true_overlapping_don,true_overlapping_acc=true_overlapping_acc,true_overlapping_don
                    mapped_overlapping_don,mapped_overlapping_acc=mapped_overlapping_acc,mapped_overlapping_don
                # compare sets and calc TP,FP,FN
                spl_TP,spl_FP,spl_FN=set_compare(true_splicing, mapped_splicing)
                don_TP,don_FP,don_FN=set_compare(true_overlapping_don, mapped_overlapping_don)
                acc_TP,acc_FP,acc_FN=set_compare(true_overlapping_acc, mapped_overlapping_acc)
                # write results
                print("\t".join([str(x) for x in [
                    1 if is_converted_bam else 0,
                    mapper,
                    cr,
                    tid,
                    intronID,
                    len(spl_TP), len(spl_FP), len(spl_FN),
                    len(don_TP), len(don_FP), len(don_FN),
                    len(acc_TP), len(acc_FP), len(acc_FN),
                    t.chromosome, intron['start'], intron['end'], t.strand]]), file=out_performance)
                if samout is not None:
                    # FIXME: slow. can we improve?
                    FPs=spl_FP | don_FP | acc_FP
                    FNs=spl_FN | don_FN | acc_FN
                    rit = ReadIterator(bam, dict_chr2idx, reference=t.chromosome, start=intron['start'] - 1, end=intron['end'], max_span=None, flag_filter=0)
                    for _, read in rit:
                        if read.query_name in FPs:
                           read.set_tag(tag='YC', value='255,0,0', value_type="Z")
                           samout.write(read)
                    rit = ReadIterator(truth_bam, dict_chr2idx, reference=t.chromosome, start=intron['start'] - 1, end=intron['end'], max_span=None, flag_filter=0)
                    for _, read in rit:
                        if read.query_name in FNs:
                           read.set_tag(tag='YC', value='255,255,255', value_type="Z")
                           samout.write(read)
    if samout is not None:
            samout.close()

def get_spliced_fraction(bam, tid, iv, chromosome, start, end, annots, donor = True) :
    """calc spliced fraction for SJ"""
    itor = SimulatedReadIterator(bam.fetch(contig=chromosome, start=start, stop=end), flagFilter=0)
    splicedTP = 0
    unsplicedTP = 0
    splicedFP = 0
    unsplicedFP = 0
    for read in itor:
        start = read.start
        end = read.end
        if donor:
            if start <= iv.begin:
                if read.splicing and end > iv.begin:
                    if read.trueTid == tid and read.trueSpliced:
                        splicedTP += 1
                    elif not read.trueTid in annots or not read.trueSpliced:
                        splicedFP += 1
                elif end > iv.begin:
                    if read.trueTid == tid and not read.trueSpliced:
                        unsplicedTP += 1
                    elif not read.trueTid in annots or read.trueSpliced:
                        unsplicedFP += 1
                else:
                    # Intronic
                    pass
        else :
            if end >= iv.end:
                if read.splicing and start < iv.end:
                    if read.trueTid == tid and read.trueSpliced:
                        splicedTP += 1
                    elif not read.trueTid in annots or not read.trueSpliced:
                        splicedFP += 1
                elif start < iv.end:
                    if read.trueTid == tid and not read.trueSpliced:
                        unsplicedTP += 1
                    elif not read.trueTid in annots or read.trueSpliced:
                        unsplicedFP += 1
    return splicedTP, unsplicedTP, splicedFP, unsplicedFP

def calculate_splice_site_mappability(config, m, bam_file, truth_bam_dir, out_dir):
    """ evaluate the passed BAM file """
    logging.info("Calculating mappability for %s" % bam_file)
    # parse mapper and cr from bam file name
    assert os.path.exists(bam_file) and os.path.exists(bam_file+'.bai'), "Could not find bam file (or idx) for %s" % (bam_file)
    assert '.cr' in bam_file, "Cannot parse bam file name %s" % (bam_file)
    tmp = bam_file[:-10] if bam_file.endswith('.final.bam') else bam_file[:-4] # file ends with .bam or .final.bam
    cr,mapper=tmp[tmp.find('.cr')+3:].rsplit('.', 1)
    out_file_prefix=out_dir+'/'+config['dataset_name']+'.cr'+cr+'.'+mapper
    truth_bam_file =  truth_bam_dir+'/'+config['dataset_name']+'.cr'+cr+'.truth.bam'
    truthBam = pysam.AlignmentFile(truth_bam_file, "rb")
    mappedBam = pysam.AlignmentFile(bam_file, "rb")
    is_converted_bam=float(cr)>0 # are there any conversions?
    readlen=config['readlen']
    genomeMappability = pysam.TabixFile(config['genome_mappability'], mode="r")
    # build transcript intervaltree
    tiv = m.build_transcript_iv()
    tmpBed = out_file_prefix + ".conv.SJ_mappability.bedGraph.tmp"
    with open(tmpBed, 'w') as out_bed:
        with open(out_file_prefix+'.mappability.tsv', 'w') as out_mappability:
            print("""is_converted_bam\tmapper\tcondition\ttid\tintron_id\tchromosome\tstart\tend\tstrand\t""" +
                                              """don_sj_mappability_TP\tdon_sj_mappability_FP\tdon_TP_reads\tdon_FP_reads\t""" +
                                              """acc_sj_mappability_TP\tacc_sj_mappability_FP\tacc_TP_reads\tacc_FP_reads""",
                                              file=out_mappability)
            for tid, t in m.transcripts.items():
                for intron in t.introns:
                    intronID = tid + "_" + str(intron['exon_number'])
                    iv = Interval(intron['start'] - 1, intron['end'] + 1)
                    ivTree = IntervalTree()
                    ivTree.add(iv)
                    # get tids of overlapping features
                    annots=[x.data.tid for x in tiv[t.chromosome].overlap(iv.begin, iv.end)]
                    # Donor
                    donorTruthSplicedTP, donorTruthUnsplicedTP, donorTruthSplicedFP, donorTruthUnsplicedFP = get_spliced_fraction(truthBam, tid, iv, t.chromosome, intron['start'] - 1, intron['start'], annots, True)
                    donorMappedSplicedTP, donorMappedUnsplicedTP, donorMappedSplicedFP, donorMappedUnsplicedFP = get_spliced_fraction(mappedBam, tid, iv, t.chromosome, intron['start'] - 1, intron['start'], annots, True)
                    if donorTruthSplicedTP + donorTruthUnsplicedTP > 0:
                        donorTruthFractionTP = donorTruthSplicedTP / (donorTruthSplicedTP + donorTruthUnsplicedTP)
                    else :
                        donorTruthFractionTP = 1
                    if donorTruthSplicedFP + donorTruthUnsplicedFP > 0:
                        donorTruthFractionFP = donorTruthSplicedFP / (donorTruthSplicedFP + donorTruthUnsplicedFP)
                    else :
                        donorTruthFractionFP = 1
                    if donorMappedSplicedTP + donorMappedUnsplicedTP > 0:
                        donorMappedFractionTP = donorMappedSplicedTP / (donorMappedSplicedTP + donorMappedUnsplicedTP)
                    else :
                        donorMappedFractionTP = 1
                    if donorMappedSplicedFP + donorMappedUnsplicedFP > 0:
                        donorMappedFractionFP = donorMappedSplicedFP / (donorMappedSplicedFP + donorMappedUnsplicedFP)
                    else:
                        donorMappedFractionFP = 1
                    donorMappabilityTP = 1 - abs(donorTruthFractionTP - donorMappedFractionTP)
                    donorMappabilityFP = 1 - abs(donorTruthFractionFP - donorMappedFractionFP)
                    # Acceptor
                    acceptorTruthSplicedTP, acceptorTruthUnsplicedTP, acceptorTruthSplicedFP, acceptorTruthUnsplicedFP = get_spliced_fraction(truthBam, tid, iv, t.chromosome, intron['end'], intron['end'] + 1, annots, False)
                    acceptorMappedSplicedTP, acceptorMappedUnsplicedTP, acceptorMappedSplicedFP, acceptorMappedUnsplicedFP = get_spliced_fraction(mappedBam, tid, iv, t.chromosome, intron['end'], intron['end'] + 1, annots, False)
                    if acceptorTruthSplicedTP + acceptorTruthUnsplicedTP > 0:
                        acceptorTruthFractionTP = acceptorTruthSplicedTP / (acceptorTruthSplicedTP + acceptorTruthUnsplicedTP)
                    else :
                        acceptorTruthFractionTP = 1
                    if acceptorTruthSplicedFP + acceptorTruthUnsplicedFP > 0:
                        acceptorTruthFractionFP = acceptorTruthSplicedFP / (acceptorTruthSplicedFP + acceptorTruthUnsplicedFP)
                    else :
                        acceptorTruthFractionFP = 1
                    if acceptorMappedSplicedTP + acceptorMappedUnsplicedTP > 0:
                        acceptorMappedFractionTP = acceptorMappedSplicedTP / (acceptorMappedSplicedTP + acceptorMappedUnsplicedTP)
                    else :
                        acceptorMappedFractionTP = 1
                    if acceptorMappedSplicedFP + acceptorMappedUnsplicedFP > 0:
                        acceptorMappedFractionFP = acceptorMappedSplicedFP / (acceptorMappedSplicedFP + acceptorMappedUnsplicedFP)
                    else :
                        acceptorMappedFractionFP = 1
                    acceptorMappabilityTP = 1 - abs(acceptorTruthFractionTP - acceptorMappedFractionTP)
                    acceptorMappabilityFP = 1 - abs(acceptorTruthFractionFP - acceptorMappedFractionFP)
                    # calculate genomic mappability
                    don_win_genomic = [intron['start'] - readlen, intron['start'] + readlen] if intron['strand'] == "+" else [intron['end'] - readlen, intron['end'] + readlen]
                    acc_win_genomic = [intron['start'] - readlen, intron['start'] + readlen] if intron['strand'] == "-" else [intron['end'] - readlen, intron['end'] + readlen]
                    gmap_prefix=config['gmap_prefix'] if 'gmap_prefix' in config else 'chr'
                    don_win_map = [float(x[3]) for x in genomeMappability.fetch(gmap_prefix+t.chromosome, don_win_genomic[0], don_win_genomic[1], parser=pysam.asTuple())]
                    don_win_min_map=min(don_win_map) if len(don_win_map)>0 else None
                    don_win_max_map=max(don_win_map) if len(don_win_map)>0 else None
                    acc_win_map = [float(x[3]) for x in genomeMappability.fetch(gmap_prefix+t.chromosome, acc_win_genomic[0], acc_win_genomic[1], parser=pysam.asTuple())]
                    acc_win_min_map=min(acc_win_map) if len(acc_win_map)>0 else None  
                    acc_win_max_map=max(acc_win_map) if len(acc_win_map)>0 else None                       
                    # get RNA subseq from transcript seq. NOTE: may be shorter than 2xreadlen
                    rna_seq=reverse_complement(t.transcript_seq) if intron['strand'] == "-" else t.transcript_seq 
                    l=len(rna_seq)
                    # note can be shorter than readlen if we run out of sequence!
                    don_seq_rna_ex = rna_seq[max(0,t.end-intron['end']-readlen):min(t.end-intron['end'],l)] if intron['strand'] == "-" else rna_seq[max(0,intron['start']-readlen-t.start):min(l,intron['start']-t.start)]
                    don_seq_rna_in = rna_seq[max(0,t.end-intron['end']):min(t.end-intron['end']+readlen,l)] if intron['strand'] == "-" else rna_seq[max(0,intron['start']-t.start):min(l,intron['start']+readlen-t.start)]
                    acc_seq_rna_in = rna_seq[max(0,t.end-intron['start']-readlen+1):min(t.end-intron['start']+1,l)] if intron['strand'] == "-" else rna_seq[max(0,intron['end']-readlen-t.start+1):min(l,intron['end']-t.start+1)]
                    acc_seq_rna_ex = rna_seq[max(0,t.end-intron['start']+1):min(t.end-intron['start']+readlen+1,l)] if intron['strand'] == "-" else rna_seq[max(0,intron['end']-t.start+1):min(l,intron['end']+readlen-t.start+1)]
                    # write bedgraph file with SJ mappability
                    print("\t".join([t.chromosome, str(don_win_genomic[0]), str(don_win_genomic[1]), intronID + "_donor", str(donorMappabilityTP), intron['strand']]), file = out_bed)
                    print("\t".join([t.chromosome, str(acc_win_genomic[0]), str(acc_win_genomic[1]), intronID + "_acceptor", str(acceptorMappabilityTP), intron['strand']]), file = out_bed)
                    print("\t".join([str(x) for x in [
                        1 if is_converted_bam else 0,
                        mapper,
                        cr,
                        tid,
                        intronID,
                        t.chromosome,
                        intron['start'],
                        intron['end'],
                        intron['strand'],
                        # donor
                        donorMappabilityTP,
                        donorMappabilityFP,
                        donorMappedSplicedTP + donorMappedUnsplicedTP,
                        donorMappedSplicedFP + donorMappedUnsplicedFP,
                        # acceptor
                        acceptorMappabilityTP,
                        acceptorMappabilityFP,
                        acceptorMappedSplicedTP + acceptorMappedUnsplicedTP,
                        acceptorMappedSplicedFP + acceptorMappedUnsplicedFP
                    ]]), file=out_mappability)
    # unmerged = pybedtools.BedTool(tmpBed).sort()
    # merged = unmerged.genome_coverage(bg=True,g=config['genome_chromosome_sizes']).map(unmerged, o="max")
    # f = open(out_bedGraph, "w")
    # for line in open(merged.fn):
    #     chr, start, end, ovlp, score = line.strip().split()
    #     print("\t".join([chr, start, end, score]), file = f)
    # f.close()
    # os.remove(tmpBed)

def calc_feature_overlap(config, m, mismapped_bam_file, out_dir):
    """ Calculate FP/FN overlap with transcript feature """
    # expected filename: big3_slamseq_nf.cr0.HISAT3N.mismapped.bam
    # parse mapper and cr from bam file name
    assert os.path.exists(mismapped_bam_file) and os.path.exists(mismapped_bam_file+'.bai'), "Could not find bam file (or idx) for %s" % (mismapped_bam_file)
    assert '.cr' in mismapped_bam_file, "Cannot parse bam file name %s" % (mismapped_bam_file)
    assert mismapped_bam_file.endswith('.mismapped.bam')
    tmp = mismapped_bam_file[:-len('.mismapped.bam')]
    cr,mapper=tmp[tmp.find('.cr')+3:].rsplit('.', 1)
    out_file_prefix=out_dir+'/'+config['dataset_name']+'.cr'+cr+'.'+mapper
    # calculate list of features
    dict_chr2idx, dict_idx2chr, dict_chr2len=get_chrom_dicts_from_bam(mismapped_bam_file)
    chrom2feat=OrderedDict()
    
    for tid,t in m.transcripts.items():
        c=t.chromosome
        if c not in chrom2feat:
            chrom2feat[c]=[]
        exon_len=0
        for exon in t.exons:
            fid = tid + "_ex" + str(exon['exon_number'])
            loc=Location(dict_chr2idx[t.chromosome], exon['start'], exon['end'], strand=t.strand)
            chrom2feat[c].append((loc,(tid, 'exon', fid)))
            exon_len+=len(loc)
        intron_len=0
        for intron in t.introns:
            loc=Location(dict_chr2idx[t.chromosome], intron['start'], intron['end'], strand=t.strand)
            fid = tid + "_" + str(intron['exon_number'])
            chrom2feat[c].append((loc,(tid, 'intron', fid)))
            intron_len+=len(loc)
    # count FP/FN
    feat_counter_fp=Counter()
    feat_counter_fn=Counter()
    samin = pysam.AlignmentFile(mismapped_bam_file, "rb")
    for c, c_len in dict_chr2len.items():
        if c not in samin.references:
            continue # no data on this chrom 
        if c not in chrom2feat or len(chrom2feat[c])==0:
            continue # no annos on this chrom   
        aits = [BlockLocationIterator(iter(sorted(chrom2feat[c], key=lambda x: x[0]) ))]
        rit = ReadIterator(samin, dict_chr2idx, reference=c, start=1, end=c_len, max_span=None, flag_filter=0) # max_span=m.max_ilen
        it = AnnotationOverlapIterator(rit, aits, check_read_alignment=True)
        for loc, (read, annos) in it:
            true_tid, true_strand, true_isoform, read_tag, true_chr, true_start, true_cigar, n_seqerr, n_converted, is_converted_read = read.query_name.split('_')
            #annos [[('ENSMUST00000029124.7', 'exon', 'ENSMUST00000029124.7_ex1'), ('ENSMUST00000029124.7', 'intron', 'ENSMUST00000029124.7_1')]]
            overlapping_annos = [t for (_, (_, t)) in annos[0]]
            if len(overlapping_annos)==0:
                continue
            if read.get_tag("YC")=='255,255,255':
                # count FN for true tid annos
                for tid, ftype, fid in overlapping_annos[0]:
                    if tid == true_tid:
                        feat_counter_fn[(tid, c, ftype, fid)]+=1
            elif read.get_tag("YC")=='255,0,0':
                # count FP for all overlapping annos
                for tid, ftype, fid in overlapping_annos[0]:
                    feat_counter_fp[(tid, c, ftype, fid)]+=1
    with open(out_file_prefix+'.feature_counts.tsv', 'w') as out:
        print("tid\tftype\tfid\tchromosome\tstart\tend\tFP\tFN", file=out)
        for c, c_len in dict_chr2len.items():
            if c not in chrom2feat or len(chrom2feat[c])==0:
                continue # no annos on this chrom   
            for loc, (tid, ftype, fid) in sorted(chrom2feat[c], key=lambda x: x[0]):
                print('\t'.join([str(x) for x in [tid,ftype,fid,c,loc.start,loc.end,
                                                  feat_counter_fp[(tid, c, ftype, fid)],
                                                  feat_counter_fn[(tid, c, ftype, fid)]
                                                  ]]), file=out)
        # exon length
        # overlaps_other_transcripts
        # FP/FN per feature (exon, intron)

def extract_splice_site_features(config, m, out_dir):
    """ extract splice-junction features into metadata file """
    out_file=out_dir+'/'+config['dataset_name']+'.SJ.metadata.tsv'
    genomeMappability = pysam.TabixFile(config['genome_mappability'], mode="r")
    readlen=config['readlen']

    with open(out_file, 'w') as out:
        print("""tid\tintron_id\tchromosome\tstart\tend\tstrand\t""" +
              """don_ex_A\tdon_ex_C\tdon_ex_T\tdon_ex_G\tdon_in_A\tdon_in_C\tdon_in_T\tdon_in_G\tdon_win_min_map\tdon_win_max_map\t""" +
              """acc_ex_A\tacc_ex_C\tacc_ex_T\tacc_ex_G\tacc_in_A\tacc_in_C\tacc_in_T\tacc_in_G\tacc_win_min_map\tacc_win_max_map""",
                                          file=out)
        for tid, t in m.transcripts.items():
            for intron in t.introns:
                intronID = tid + "_" + str(intron['exon_number'])

                # calculate genomic mappability
                don_win_genomic = [intron['start'] - readlen, intron['start'] + readlen] if intron['strand'] == "+" else [intron['end'] - readlen, intron['end'] + readlen]
                acc_win_genomic = [intron['start'] - readlen, intron['start'] + readlen] if intron['strand'] == "-" else [intron['end'] - readlen, intron['end'] + readlen]
                don_win_map = [float(x[3]) for x in genomeMappability.fetch('chr'+t.chromosome, don_win_genomic[0], don_win_genomic[1], parser=pysam.asTuple())]
                don_win_min_map=min(don_win_map) if len(don_win_map)>0 else None
                don_win_max_map=max(don_win_map) if len(don_win_map)>0 else None
                acc_win_map = [float(x[3]) for x in genomeMappability.fetch('chr'+t.chromosome, acc_win_genomic[0], acc_win_genomic[1], parser=pysam.asTuple())]
                acc_win_min_map=min(acc_win_map) if len(acc_win_map)>0 else None
                acc_win_max_map=max(acc_win_map) if len(acc_win_map)>0 else None
                # get RNA subseq from transcript seq. NOTE: may be shorter than 2xreadlen
                rna_seq=reverse_complement(t.transcript_seq) if intron['strand'] == "-" else t.transcript_seq
                l=len(rna_seq)
                # note can be shorter than readlen if we run out of sequence!
                don_seq_rna_ex = rna_seq[max(0,t.end-intron['end']-readlen):min(t.end-intron['end'],l)] if intron['strand'] == "-" else rna_seq[max(0,intron['start']-readlen-t.start):min(l,intron['start']-t.start)]
                don_seq_rna_in = rna_seq[max(0,t.end-intron['end']):min(t.end-intron['end']+readlen,l)] if intron['strand'] == "-" else rna_seq[max(0,intron['start']-t.start):min(l,intron['start']+readlen-t.start)]
                acc_seq_rna_in = rna_seq[max(0,t.end-intron['start']-readlen+1):min(t.end-intron['start']+1,l)] if intron['strand'] == "-" else rna_seq[max(0,intron['end']-readlen-t.start+1):min(l,intron['end']-t.start+1)]
                acc_seq_rna_ex = rna_seq[max(0,t.end-intron['start']+1):min(t.end-intron['start']+readlen+1,l)] if intron['strand'] == "-" else rna_seq[max(0,intron['end']-t.start+1):min(l,intron['end']+readlen-t.start+1)]

                # write metadata file
                print("\t".join([str(x) for x in [
                    tid,
                    intronID,
                    t.chromosome,
                    intron['start'],
                    intron['end'],
                    intron['strand'],
                    # donor
                    don_seq_rna_ex.count("A"),
                    don_seq_rna_ex.count("C"),
                    don_seq_rna_ex.count("T"),
                    don_seq_rna_ex.count("G"),
                    don_seq_rna_in.count("A"),
                    don_seq_rna_in.count("C"),
                    don_seq_rna_in.count("T"),
                    don_seq_rna_in.count("G"),
                    don_win_min_map,
                    don_win_max_map,
                    # acceptor
                    acc_seq_rna_ex.count("A"),
                    acc_seq_rna_ex.count("C"),
                    acc_seq_rna_ex.count("T"),
                    acc_seq_rna_ex.count("G"),
                    acc_seq_rna_in.count("A"),
                    acc_seq_rna_in.count("C"),
                    acc_seq_rna_in.count("T"),
                    acc_seq_rna_in.count("G"),
                    acc_win_min_map,
                    acc_win_max_map
                ]]), file=out)


def extract_transcript_features(config, m, out_dir):
    """ Calculate FP/FN overlap with transcript feature """
    out_file=out_dir+'/'+config['dataset_name']+'.transcript.metadata.tsv'
    # calculate list of features
    dict_chr2idx, dict_idx2chr, dict_chr2len=get_chrom_dicts(config["genome_fa"])
    tid2info={}
    # build transcript intervaltree
    tiv = m.build_transcript_iv()
    for c, c_len in dict_chr2len.items():
        for tid,t in m.transcripts.items():
            exon_len=0
            for exon in t.exons:
                loc=Location(dict_chr2idx[t.chromosome], exon['start'], exon['end'], strand=t.strand)
                exon_len+=len(loc)
            intron_len=0
            for intron in t.introns:
                loc=Location(dict_chr2idx[t.chromosome], intron['start'], intron['end'], strand=t.strand)
                intron_len+=len(loc)
            overlapping=[x.data.tid for x in tiv[t.chromosome].overlap(t.start, t.end) if x.data.tid != tid]
            tid2info[tid]=(exon_len, intron_len, overlapping)
    # write result table
    with open(out_file, 'w') as out:
        print("tid\texon_len\tintron_len\tn_overlapping\toverapping_tids", file=out)
        for tid,(exon_len, intron_len, overlapping) in tid2info.items():
            print('\t'.join([str(x) for x in [tid,exon_len,intron_len,len(overlapping),','.join(overlapping)]]), file=out)

def write_parquet_table(tsv_file, partition_cols, out_dir, add_columns=None, add_values=None, ren_columns=None, final_columns=None):
    """ Convert TSV to parquet. Constant columns can be added by providing a list of column names (add_columns) and values (add_values)"""
    tab = pd.read_csv(tsv_file,
                          delimiter='\t',
                          comment='#',
                          encoding='utf-8',
                          float_precision='high',
                          quoting=csv.QUOTE_NONE,
                          low_memory=False,
                          error_bad_lines=True,
                          dtype={'condition':  'object'}
                          )  
    if add_columns is not None:
        for i,col in enumerate(add_columns):
            tab.insert(0, col, [add_values[i]] * len(tab.index), True)
    if ren_columns is not None:
        tab.rename(columns=ren_columns, inplace=True)
    if final_columns is not None:
        tab.columns=final_columns
    tab.to_parquet(out_dir, partition_cols=partition_cols)

def write_parquet_db(config, in_dir, out_dir):  
    """ Create a result parquet database """
    for data_file in glob.glob("%s/eva/overall_performance/*.tid_performance.tsv" % in_dir):
        write_parquet_table(data_file, ['mapper'], out_dir+'/overall_performance')
    for data_file in glob.glob("%s/eva/splice_site_mappability/*.mappability.tsv" % in_dir):
        write_parquet_table(data_file, ['mapper'], out_dir+'/splice_site_mappability')
    for data_file in glob.glob("%s/eva/splice_site_performance/*.splice_site_performance.tsv" % in_dir):
        write_parquet_table(data_file, ['mapper'], out_dir+'/splice_site_performance')
    for data_file in glob.glob("%s/eva/feature_overlap/*.feature_counts.tsv" % in_dir):
        tmp=data_file[len(config['dataset_name'])+1:-len('.feature_counts.tsv')]
        cr,mapper=tmp[tmp.find('.cr')+3:].rsplit('.', 1)
        write_parquet_table(data_file, ['mapper', 'chromosome', 'ftype'], out_dir+'/feature_counts_mismapped', add_columns=['condition_id', 'mapper'], add_values=[cr, mapper])
    # write featureCounts output
    for data_file in glob.glob("%s/eva/feature_counts/*.featureCounts.txt" % (in_dir)):
        tmp=data_file[len("%s/eva/feature_counts/"%in_dir + config['dataset_name'])+1:-len('.featureCounts.txt')]
        tmp,ftype=tmp.rsplit('.', 1)
        tmp,mapper=tmp.rsplit('.', 1)
        if mapper=='final':
            tmp,mapper=tmp.rsplit('.', 1)
        cr=tmp[len('cr'):]
        write_parquet_table(data_file, ['mapper', 'ftype'], 
                            out_dir+'/feature_counts', 
                            add_columns=['condition_id', 'mapper', 'ftype'], 
                            add_values=[cr, mapper, ftype],
                            final_columns=['ftype', 'mapper', 'condition_id', 
                                           'Geneid', 'Chr', 'Start', 'End',
                                           'Strand', 'Length','read_count']) # TODO this is ugly, improve!
        
    write_parquet_table("%s/eva/meta/%s.SJ.metadata.tsv" % (in_dir, config['dataset_name']), [], out_dir+'/sj_metadata')
    write_parquet_table("%s/eva/meta/%s.transcript.metadata.tsv" % (in_dir, config['dataset_name']), [], out_dir+'/tx_metadata')
    write_parquet_table("%s/sim/reference_model/gene_anno.tsv.gz" % in_dir, [], out_dir+'/gene_anno')
    
    print("All done.")
    

# if __name__ == '__main__':
#     config_file='/Volumes/groups/ameres/Niko/projects/Ameres/splicing/splice_sim/testruns/big3_slamseq_nf/splice_sim.config.json'
#     model_file='/Volumes/groups/ameres/Niko/projects/Ameres/splicing/splice_sim/testruns/big3_slamseq_nf/sim/reference_model/big3_slamseq_nf.model'
#     out_dir='/Volumes/groups/ameres/Niko/projects/Ameres/splicing/splice_sim/testruns/big3_slamseq_nf/sim/reference_model/delme/'
#     bam_file='/Volumes/groups/ameres/Niko/projects/Ameres/splicing/splice_sim/testruns/big3_slamseq_nf/sim/final_bams/big3_slamseq_nf.cr0.1.STAR.final.bam'
#     config=json.load(open(config_file), object_pairs_hook=OrderedDict)
#     config=localize_config(config)
#     m = Model.load_from_file(model_file)
#     print(m)
    #evaluate_splice_sites_performance(config, m, bam_file,out_dir)
    #truth_bam_dir='/Volumes/groups/ameres/Niko/projects/Ameres/splicing/splice_sim/testruns/nf3/sim/bams_truth'
    #calculate_splice_site_mappability(config, m, bam_file, truth_bam_dir, out_dir)
