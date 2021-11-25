'''
@author: niko.popitsch@imba.oeaw.ac.at, tobias.neumann@imp.ac.at
'''

from collections import *
import csv, datetime, time, logging, sys, os, json
import pysam
from genomic_iterators import Location, get_chrom_dicts_from_bam
from genomic_iterators.bam_iterators import ReadIterator
from genomic_iterators.grouped_iterators import BlockedIterator, BlockLocationIterator, AnnotationOverlapIterator
from splice_sim import *

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
        chrom2trans[t.Chromosome].append((Location(dict_chr2idx[t.Chromosome], t.Start, t.End, strand=t.Strand),tid))
    
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
    
def evaluate_overall_performance(config, m, final_bam_dir, truth_bam_dir, out_dir):  
    # write read data?
    write_reads=config['write_reads'] if 'write_reads' in config else False
    if write_reads:
        out_reads=open(out_dir+'/reads.tsv', 'w') 
        print("mapped_coords\ttrue_coords\tclassification\ttid\tis_converted_bam\tmapper\tcondition_id\toverlap\ttrue_tid\ttrue_strand\ttrue_isoform\ttag\ttrue_chr\tn_true_seqerr\tn_tc_pos\tis_converted_read\toverlapping_tids\tread_name", file=out_reads)
    else:
        out_reads=None
    logging.info("write_reads: %s" % str(write_reads))
    with open(out_dir+'/tid_performance.tsv', 'w') as out_performance:
        print("is_converted_bam\tmapper\tcondition_id\tiso\ttid\tTP\tFP\tFN", file=out_performance)
        for cr in [0]+m.condition.conversion_rates:
            # write truth
            bam_file=truth_bam_dir+'/'+config["dataset_name"]+".cr"+str(cr)+".truth.bam"
            evaluate_bam(bam_file, None, m, 'truth', cr, out_performance, None)
            # write mapper stats
            for mapper in config['mappers']:
                bam_file=final_bam_dir+'/'+config["dataset_name"]+".cr"+str(cr)+"."+mapper+".final.bam"
                assert os.path.exists(bam_file) and os.path.exists(bam_file+'.bai'), "Could not find bam file (or idx) for mapper %s/cr %f: %s" % (mapper, cr, bam_file)                
                bam_out_file=out_dir+'/'+config["dataset_name"]+".cr"+str(cr)+"."+mapper+".mismapped.bam"
                evaluate_bam(bam_file, bam_out_file, m, mapper, cr, out_performance, out_reads)
    if write_reads:
        out_reads.close()
