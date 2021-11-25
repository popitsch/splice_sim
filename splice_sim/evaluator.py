'''
@author: niko.popitsch@imba.oeaw.ac.at, tobias.neumann@imp.ac.at
'''

from collections import *
import csv, datetime, time, logging, sys, os, json
import pysam
import pyranges as pr
import pandas as pd
import numpy as np
from scipy import stats, special
from sklearn.metrics import mean_squared_error
from utils import *
import random
import math
import gzip
import pybedtools
from intervaltree import Interval, IntervalTree
import logging
from model import Model, Condition, Isoform, Transcript
from iterator import *
from Bio.Seq import reverse_complement
from genomic_iterators import Location

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
    is_converted_bam=cr>0 # are there any concversions?
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

def calc_overall_performance(bam_file, bam_out, is_converted_bam, m, mapper, cr, out_reads, out_performance):
    """ evaluate the passed BAM file """
    logging.info("Evaluating %s" % bam_file)
    performance = Counter()
    dict_chr2idx, dict_idx2chr, dict_chr2len=get_chrom_dicts_from_bam(bam_file)
    n_reads=0
    samin = pysam.AlignmentFile(bam_file, "rb") 
    samout = pysam.AlignmentFile(bam_out+'.tmp.bam', "wb", template=samin ) if bam_out is not None else None
    
    chrom2trans={}
    for c, c_len in dict_chr2len.items():
        chrom2trans[c]=[]
    for t in m.transcripts.values():
        chrom2trans[t.Chromosome].append(Location(dict_chr2idx[t.Chromosome], t.Start, t.End, strand=t.Strand))
    
    for c, c_len in dict_chr2len.items():
        if c not in samin.references:
            continue # no data on this chrom 
        if len(chrom2trans[c])==0: # no transcript on this chrom: all reads are FN
            rit = ReadIterator(bam_file, dict_chr2idx, reference=c, start=1, end=c_len, max_span=None, flag_filter=0) # max_span=m.max_ilen
            for loc, read in rit:
                n_reads+=1
                performance = classify_read(read, [], is_converted_bam, mapper, cr, dict_chr2idx, performance, out_reads, samout)       
        else:
            aits = [BlockLocationIterator(iter(chrom2trans[c]))]
            rit = ReadIterator(bam_file, dict_chr2idx, reference=c, start=1, end=c_len, max_span=None, flag_filter=0) # max_span=m.max_ilen
            it = AnnotationOverlapIterator(rit, aits)
            for loc, (read, annos) in it:
                n_reads+=1
                overlapping_tids = [t[0] for (_, (_, t)) in annos[0]]
                performance = classify_read(read, overlapping_tids, is_converted_bam, mapper, cr, dict_chr2idx, performance, out_reads, samout)
             
    # add unmapped reads
    for read in samin.fetch(contig=None, until_eof=True):
        if not read.is_unmapped:
            continue
        n_reads+=1
        performance = classify_read(read, [], is_converted_bam, mapper, cr, dict_chr2idx, performance, out_reads, samout)
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
        sambambasortbam(bam_out+'.tmp.bam', 
                        bam_out, 
                        index=True, 
                        override=True, 
                        delinFile=True)

    
    
def evaluate_overall_performance(m, final_bam_dir, out_dir):   
    for cr in [0]+m.condition.conversion_rates:
        # convert truth files
        for mapper in m.config['mappers']:
            bam_file=final_bam_dir+'/'+config["dataset_name"]+".cr"+str(cr)+"."+mapper+".bam"
            print(bam_file)
            evaluate_bam(bam_ori, bam_ori_out, cr, m, mapper, out, out2)


model_file='/Volumes/groups/ameres/Niko/projects/Ameres/splicing/splice_sim/testruns/nf1/sim/reference_model/nf1.model'
final_bam_dir='/Volumes/groups/ameres/Niko/projects/Ameres/splicing/splice_sim/testruns/nf1/sim/final_bams/'
out_dir='/Volumes/groups/ameres/Niko/projects/Ameres/splicing/splice_sim/testruns/nf1/eva/
m = Model.load_from_file(model_file)
evaluate_overall_performance(n, final_bam_dir, out_dir)