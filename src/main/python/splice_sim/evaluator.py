'''
@author: niko.popitsch@imba.oeaw.ac.at
'''

from collections import *
import csv, datetime, time, logging, sys, os, json
import pysam
import pyranges as pr
import pandas as pd
import numpy as np
from utils import *
import random
import math
import gzip
import logging
from model import Model, Condition, Isoform, Transcript
from iterator import *

# config = json.load(open('/Volumes/groups/ameres/Niko/projects/Ameres/splicing/splice_sim/testruns/big3/config.json'), object_pairs_hook=OrderedDict)
# config["transcripts"] = json.load(open('/Volumes/groups/ameres/Niko/projects/Ameres/splicing/splice_sim/testruns/big3/' + config['transcript_data']), object_pairs_hook=OrderedDict)
#  
# for k in ['gene_gff', 'genome_fa']:
#     config[k] = '/Volumes' + config[k] 
# m = Model(config)
# bam_file = '/Volumes/groups/ameres/Niko/projects/Ameres/splicing/splice_sim/testruns/big3/big3/sim/bam_tc/STAR/big3.c_01.STAR.TC.bam'
# # anno='/Volumes/groups/ameres/Niko/projects/Ameres/splicing/splice_sim/testruns/small4/gencode.SMALL.gff3.gz'
#  
# fasta = pysam.FastaFile(config["genome_fa"])
# chromosomes = fasta.references
# dict_chr2idx = {k: v for v, k in enumerate(chromosomes)}
# dict_idx2chr = {v: k for v, k in enumerate(chromosomes)}
#  
# # dict_idx2chr = {v: k for v, k in enumerate(chromosomes)}
# dict_chr2len = {c: fasta.get_reference_length(c) for c in fasta.references}
# x=Counter()
# with open('/Volumes/groups/ameres/Niko/projects/Ameres/splicing/splice_sim/testruns/big3/delme.tsv', 'w') as out:
#     for c, c_len in dict_chr2len.items():
#         print(c)
#         rit = ReadIterator(bam_file, dict_chr2idx, reference=c, start=1, end=c_len, max_span=None, flag_filter=0)
#         for _,r in rit:
#             x['reads']+=1
# print(x)
#         
#         df = m.df[(m.df.Feature == 'transcript') & (m.df.Chromosome == c)]
#         if not df.empty:
#             print("chrom", c)
#             aits = [BlockLocationIterator(PyrangeIterator(df, dict_chr2idx, 'transcript_id'))]
#             rit = ReadIterator(bam_file, dict_chr2idx, reference=c, start=1, end=l, max_span=100000)
#             it = AnnotationOverlapIterator(rit, aits)
#             for loc, (read, annos) in it:
#                 tids = [t[0] for (_, (_, t)) in annos[0]]
#                 print(loc.to_str(dict_idx2chr), read.query_name, tids, file=out)
# BlockLocationIterator()
 
# [(14:142903115-142906754, ([14:142903115-142906754], ['ENSMUST00000100497.10'])), 
#   (14:142903501-142906702, ([14:142903501-142906702], ['ENSMUST00000167721.7']))]

def classify_read(read, overlapping_tids, is_converted, mapper, condition, dict_chr2idx, performance, out_reads, sam_out):
    """ Classifies a read """
    overlap_cutoff = 10 #0.8 * m.readlen 
    read_name = read.query_name
    true_tid, true_strand, true_isoform, tag, true_chr, true_read_cigar, true_seqerr, tc_pos = read_name.split("_")
    n_true_seqerr = len(true_seqerr.split(',')) if true_seqerr != 'NA' else 0
    n_tc_pos = len(tc_pos.split(',')) if tc_pos != 'NA' else 0                
    read_chr = 'NA' if read.is_unmapped else read.reference_name
    read_coord = 'NA' if read.is_unmapped else '%s:%i-%i' % (read.reference_name, read.reference_start, read.reference_end)
    true_tuples = [tuple([int(y) for y in x.split('-')]) for x in true_read_cigar.split(',')]
    read_tuples = read.get_blocks()
    overlap = calc_coverage(true_chr, read_chr, true_tuples, read_tuples)
    true_coord = "%s:%i-%i" % ( true_chr, true_tuples[0][0],  true_tuples[-1][1])
    
    # read strongly overlaps with real location; FIXME: check strand!
    if overlap >= overlap_cutoff:
        performance[true_tid, true_isoform, 'TP'] += 1
    else:
        performance[true_tid, true_isoform, 'FN'] += 1
        print('\t'.join([str(x) for x in (read_coord, true_coord, 'FN', 'NA', 
                                          1 if is_converted else 0, mapper, condition, overlap, true_tid, 
                                          true_strand, true_isoform, tag, true_chr,n_true_seqerr, n_tc_pos, 
                                          ','.join(overlapping_tids) if len(overlapping_tids) > 0 else 'NA', 
                                          read_name)]), file=out_reads)
        for tid in overlapping_tids:
            performance[tid, true_isoform, 'FP'] += 1/len(overlapping_tids)
            print('\t'.join([str(x) for x in (read_coord, true_coord, 'FP', tid, 
                                              1 if is_converted else 0, mapper, condition, overlap, true_tid, 
                                              true_strand, true_isoform, tag, true_chr,n_true_seqerr, n_tc_pos, 
                                              ','.join(overlapping_tids) if len(overlapping_tids) > 0 else 'NA', 
                                              read_name)]), file=out_reads)
        # write FN/FP reads
        if sam_out is not None:
            true_read = pysam.AlignedSegment()
            true_read.query_name = read_name
            true_read.query_sequence = read.query_sequence
            true_read.query_qualities = read.query_qualities
            true_read.reference_start = true_tuples[0][0]-1
            true_read.cigartuples=create_cigatuples(true_tuples)
            true_read.reference_id = dict_chr2idx[true_chr]
            true_read.tags = (("NM", n_true_seqerr + n_tc_pos ),)
            true_read.flag = 0 if true_strand=='+' else 16
            true_read.set_tag(tag='YC', value='255,255,255', value_type="Z")
            sam_out.write(true_read)
            read.set_tag(tag='YC', value='255,0,0', value_type="Z")
            sam_out.write(read)
    return performance

def evaluate_bam(bam_file, bam_out, is_converted, m, mapper, condition, out_reads, out_performance):
    """ evaluate the passed BAM file """
    logging.info("Evaluating %s" % bam_file)
    performance = Counter()
    chromosomes = m.genome.references
    dict_chr2idx = {k: v for v, k in enumerate(chromosomes)}
    dict_idx2chr = {v: k for v, k in enumerate(chromosomes)}
    dict_chr2len = {c: m.genome.get_reference_length(c) for c in m.genome.references}
    n_reads=0
    samin = pysam.AlignmentFile(bam_file, "rb") 
    samout = pysam.AlignmentFile(bam_out, "wb", template=samin ) if bam_out is not None else None
    
    for c, c_len in dict_chr2len.items():
        df = m.df[(m.df.Feature == 'transcript') & (m.df.Chromosome == c)]  # get annotations
        #print("Processing chromosome %s of %s/%s/%s"  % (c,  is_converted, mapper, condition) )
        if df.empty: # no transcript on this chrom: all reads are FN
            rit = ReadIterator(bam_file, dict_chr2idx, reference=c, start=1, end=c_len, max_span=None, flag_filter=0) # max_span=m.max_ilen
            for loc, read in rit:
                n_reads+=1
                performance = classify_read(read, [], is_converted, mapper, condition, dict_chr2idx, performance, out_reads, samout)       
        else:
            aits = [BlockLocationIterator(PyrangeIterator(df, dict_chr2idx, 'transcript_id'))]
            rit = ReadIterator(bam_file, dict_chr2idx, reference=c, start=1, end=c_len, max_span=None, flag_filter=0) # max_span=m.max_ilen
            it = AnnotationOverlapIterator(rit, aits)
            for loc, (read, annos) in it:
                n_reads+=1
                overlapping_tids = [t[0] for (_, (_, t)) in annos[0]]
                performance = classify_read(read, overlapping_tids, is_converted, mapper, condition, dict_chr2idx, performance, out_reads, samout)
    tids= set([x for x, y, z in performance.keys()])
    isos= set([y for x, y, z in performance.keys()])
    for tid in tids:
        for iso in isos:
             print("\t".join([str(x) for x in [
                 1 if is_converted else 0, 
                 mapper, 
                 condition, 
                 iso, 
                 tid, 
                 performance[tid, iso, 'TP'], 
                 performance[tid, iso, 'FP'], 
                 performance[tid, iso, 'FN'] ]]), file=out_performance) 
             
    # add unmapped reads
    for read in samin.fetch(contig=None, until_eof=True):
        n_reads+=1
        performance = classify_read(read, [], is_converted, mapper, condition, dict_chr2idx, performance, out_reads, samout)
    print("%s reads:  %i %i %i" % (bam_file, n_reads, samin.mapped, samin.unmapped))
    samin.close()
    if samout is not None:
        samout.close()
    
# #bam='/Volumes/groups/ameres/Niko/projects/Ameres/splicing/splice_sim/testruns/small4/small4/sim/bam_tc/STAR/small4.5min.STAR.TC.bam'
# bam='/Volumes/groups/ameres/Niko/projects/Ameres/splicing/splice_sim/testruns/small4/small4/sim/bam_tc/HISAT2_TLA/small4.5min.HISAT2_TLA.TC.bam'
# evaluate_bam(bam, 'test', 'map', 'cond', None)

     
def evaluate_dataset(config, config_dir, simdir, outdir, overwrite=False):
    """ Evaluates a dataset. """
    startTime = time.time()

    startTime = time.time()
    logging.info("Evaluating dataset %s" % config['dataset_name'])
    
    # output dirs
    if not os.path.exists(outdir):
        print("Creating dir " + outdir)
        os.makedirs(outdir)
    

    # read transcript data from external file if not in config
    if 'transcripts' in config:
        logging.info("Reading transcript configuration from config file.")
    else:
        assert "transcript_data" in config, "Transcript data needs to be configured either in config file ('transcripts' section) or in an external file referenced via 'transcript_data'"
        tfile = config['transcript_data']
        if not os.path.isabs(tfile):
            tfile = config_dir + "/" + tfile
        assert files_exist(tfile), "transcript_data file not found: %s" % tdata
        logging.info("Reading transcript configuration from external config file %s" % tfile)
        tdata = json.load(open(tfile), object_pairs_hook=OrderedDict)
        config["transcripts"] = tdata
        
    # instantiate model
    m = Model(config)
    logging.info("Configured conditions: %s" % (", ".join(str(x) for x in m.conditions)))
     
    fout = outdir + 'reads.tsv'
    fout2 = outdir + 'tid_performance.tsv'
    with open(fout, 'w') as out:
        with open(fout2, 'w') as out2:
            print("mapped_coords\ttrue_coords\tclassification\ttid\tis_converted\tmapper\tcondition\toverlap\ttrue_tid\ttrue_strand\ttrue_isoform\ttag\ttrue_chr\tn_true_seqerr\tn_tc_pos\toverlapping_tids\tread_name", file=out)
            print("is_converted\tmapper\tcondition\tiso\ttid\tTP\tFP\tFN", file=out2)
            for cond in m.conditions:
                for mapper in config['mappers'].keys():
                    bam_ori = simdir + "bam_ori/" + mapper + "/" + config['dataset_name'] + "." + cond.id + "." + mapper + ".bam"
                    bam_tc = simdir + "bam_tc/" + mapper + "/" + config['dataset_name'] + "." + cond.id + "." + mapper + ".TC.bam"
                    bam_out_dir = outdir + "bam_ori/" + mapper + "/"
                    if not os.path.exists(bam_out_dir):
                        print("Creating dir " + bam_out_dir)
                        os.makedirs(bam_out_dir) 
                    bam_ori_out = bam_out_dir + config['dataset_name'] + "." + cond.id + "." + mapper + ".mismapped.bam"
                    bam_tc_out = bam_out_dir + config['dataset_name'] + "." + cond.id + "." + mapper + ".TC.mismapped.bam"
                    evaluate_bam(bam_ori, bam_ori_out, False, m, mapper, cond.id, out, out2)    
                    evaluate_bam(bam_tc, bam_tc_out, True, m, mapper, cond.id, out, out2)    
                # eval truth
                bam_ori_truth  = simdir + "bam_ori/TRUTH/" + config['dataset_name'] + "." + cond.id + ".TRUTH.bam"
                bam_tc_truth   = simdir + "bam_tc/TRUTH/"  + config['dataset_name'] + "." + cond.id + ".TRUTH.TC.bam"
                evaluate_bam(bam_ori_truth, None, False, m, 'NA', cond.id, out, out2)    
                evaluate_bam(bam_tc_truth, None, True, m, 'NA', cond.id, out, out2)    
    bgzip(fout, delinFile=True, override=True)
    logging.info("All done in %s" % str(datetime.timedelta(seconds=time.time() - startTime)))
    print("All done in", datetime.timedelta(seconds=time.time() - startTime))
    
