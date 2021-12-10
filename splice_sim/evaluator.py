'''
@author: niko.popitsch@imba.oeaw.ac.at, tobias.neumann@imp.ac.at
'''
from collections import *
import csv, datetime, time, logging, sys, os, json, pickle
import pysam
from genomic_iterators import Location, get_chrom_dicts_from_bam
from genomic_iterators.bam_iterators import ReadIterator
from genomic_iterators.grouped_iterators import BlockedIterator, BlockLocationIterator, AnnotationOverlapIterator
from splice_sim.utils import *
from splice_sim.model import *
from splice_sim.iterator import *


def basename(f):
    return os.path.splitext(os.path.basename(f))[0]

def calc_coverage(true_chr, read_chr, true_pos, read_pos):
    """ Calculate the coverage (number of bases) of a read and its simulation origin """
    if true_chr != read_chr:
        return 0
    cov = len(true_pos & read_pos) / len(true_pos)
    return cov

def classify_read(read, overlapping_tids, cr, mapper, dict_chr2idx, performance, out_reads, sam_out, overlap_cutoff = 0.8):
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
    
    overlap = calc_coverage(true_chr, read_chr, true_pos, read_pos)
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
        for tid in overlapping_tids:
            performance[tid, true_isoform, 'FP'] += 1/len(overlapping_tids)
            if out_reads is not None:
                print('\t'.join([str(x) for x in (read_coord, true_coord, 'FP', tid, 
                                                  1 if is_converted_bam else 0, mapper, 
                                                  cr, 
                                                  overlap, true_tid, 
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
    chrom2trans={}
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
                performance = classify_read(read, [], cr, mapper, dict_chr2idx, performance, out_reads, samout)
        else:
            aits = [BlockLocationIterator(iter(chrom2trans[c]))]
            rit = ReadIterator(bam_file, dict_chr2idx, reference=c, start=1, end=c_len, max_span=None, flag_filter=0) # max_span=m.max_ilen
            it = AnnotationOverlapIterator(rit, aits)
            for loc, (read, annos) in it:
                n_reads+=1
                overlapping_tids = [t[0] for (_, (_, t)) in annos[0]]
                performance = classify_read(read, overlapping_tids, cr, mapper, dict_chr2idx, performance, out_reads, samout)             
    # add unmapped reads
    for read in samin.fetch(contig=None, until_eof=True):
        if not read.is_unmapped:
            continue
        n_reads+=1
        performance = classify_read(read, [], cr, mapper, dict_chr2idx, performance, out_reads, samout)
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
                        if read.splicing and not read.hasSpliceSite(iv) and start < iv.end:
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
    genomeMappability = pysam.TabixFile(config['genome_mappability'], mode="r")
    is_converted_bam=float(cr)>0 # are there any conversions?
    readlen=config['readlen']
    # build transcript intervaltree
    tiv = m.build_transcript_iv()
    tmpBed = out_file_prefix + ".conv.SJ_mappability.bedGraph.tmp"
    with open(tmpBed, 'w') as out_bed:
        with open(out_file_prefix+'.mappability.tsv', 'w') as out_mappability:
            print("""is_converted_bam\tmapper\tcondition\ttid\tintron_id\tchromosome\tstart\tend\tstrand\t""" +
                                              """don_sj_mappability_TP\tdon_sj_mappability_FP\tdon_TP_reads\tdon_FP_reads\tdon_ex_A\tdon_ex_C\tdon_ex_T\tdon_ex_G\tdon_in_A\tdon_in_C\tdon_in_T\tdon_in_G\tdon_win_min_map\tdon_win_max_map\t""" +
                                              """acc_sj_mappability_TP\tacc_sj_mappability_FP\tacc_TP_reads\tacc_FP_reads\tacc_ex_A\tacc_ex_C\tacc_ex_T\tacc_ex_G\tacc_in_A\tacc_in_C\tacc_in_T\tacc_in_G\tacc_win_min_map\tacc_win_max_map""", 
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
                        acceptorMappabilityTP,
                        acceptorMappabilityFP,
                        acceptorMappedSplicedTP + acceptorMappedUnsplicedTP,
                        acceptorMappedSplicedFP + acceptorMappedUnsplicedFP,
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
                    ]]), file=out_mappability)
    # unmerged = pybedtools.BedTool(tmpBed).sort()
    # merged = unmerged.genome_coverage(bg=True,g=config['genome_chromosome_sizes']).map(unmerged, o="max")
    # f = open(out_bedGraph, "w")
    # for line in open(merged.fn):
    #     chr, start, end, ovlp, score = line.strip().split()
    #     print("\t".join([chr, start, end, score]), file = f)
    # f.close()
    # os.remove(tmpBed)           
# if __name__ == '__main__':
#     config_file='/Volumes/groups/ameres/Niko/projects/Ameres/splicing/splice_sim/testruns/nf3/splice_sim.config.json'
#     model_file='/Volumes/groups/ameres/Niko/projects/Ameres/splicing/splice_sim/testruns/nf3/tmp/nf2.model'
#     out_dir='/Volumes/groups/ameres/Niko/projects/Ameres/splicing/splice_sim/testruns/nf3/tmp/'
#     bam_file='/Volumes/groups/ameres/Niko/projects/Ameres/splicing/splice_sim/testruns/nf3/sim/final_bams/nf2.cr0.STAR.final.bam'
#     config=json.load(open(config_file), object_pairs_hook=OrderedDict)
#     config=localize_config(config)
#     m = Model.load_from_file(model_file)
#     #evaluate_splice_sites_performance(config, m, bam_file,out_dir)
#     truth_bam_dir='/Volumes/groups/ameres/Niko/projects/Ameres/splicing/splice_sim/testruns/nf3/sim/bams_truth'
#     calculate_splice_site_mappability(config, m, bam_file, truth_bam_dir, out_dir)