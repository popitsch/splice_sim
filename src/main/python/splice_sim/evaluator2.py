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

def basename(f):
    return os.path.splitext(os.path.basename(f))[0]

def calc_coverage(true_chr, read_chr, true_pos, read_pos):
    """ Calculate the coverage (number of bases) of a read and its simulation origin """
    if true_chr != read_chr:
        return 0
    cov = len(true_pos & read_pos) / len(true_pos)
    return cov

def classify_read(read, overlapping_tids, is_converted_bam, mapper, condition, dict_chr2idx, performance, out_reads, sam_out, overlap_cutoff = 0.8):
    """ Classifies a read """
    read_name = read.query_name
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
                                          condition.id,  
                                          overlap, true_tid, 
                                          true_strand, true_isoform, read_tag, true_chr,n_seqerr, n_converted,
                                          is_converted_read, 
                                          ','.join(overlapping_tids) if len(overlapping_tids) > 0 else 'NA', 
                                          read_name)]), file=out_reads)
        for tid in overlapping_tids:
            performance[tid, true_isoform, 'FP'] += 1/len(overlapping_tids)
            print('\t'.join([str(x) for x in (read_coord, true_coord, 'FP', tid, 
                                              1 if is_converted_bam else 0, mapper, 
                                              condition.id, 
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

def evaluate_bam(bam_file, bam_out, is_converted_bam, m, mapper, condition, out_reads, out_performance):
    """ evaluate the passed BAM file """
    logging.info("Evaluating %s" % bam_file)
    performance = Counter()
    dict_chr2idx, dict_idx2chr, dict_chr2len=get_chrom_dicts_from_bam(bam_file)
    n_reads=0
    samin = pysam.AlignmentFile(bam_file, "rb") 
    samout = pysam.AlignmentFile(bam_out+'.tmp.bam', "wb", template=samin ) if bam_out is not None else None
    for c, c_len in dict_chr2len.items():
        if c not in samin.references:
            continue # no data on this chrom 
        df = m.df[(m.df.Feature == 'transcript') & (m.df.Chromosome == c)]# get annotations
        # check whether df is sorted! 
        l=[x for x in df.Start]
        assert all(l[i] <= l[i+1] for i in range(len(l)-1)), "DF is not sorted by start coord!"
        #print("Processing chromosome %s of %s/%s/%s"  % (c,  is_converted_bam, mapper, condition) )
        if df.empty: # no transcript on this chrom: all reads are FN
            rit = ReadIterator(bam_file, dict_chr2idx, reference=c, start=1, end=c_len, max_span=None, flag_filter=0) # max_span=m.max_ilen
            for loc, read in rit:
                n_reads+=1
                performance = classify_read(read, [], is_converted_bam, mapper, condition, dict_chr2idx, performance, out_reads, samout)       
        else:
            aits = [BlockLocationIterator(PyrangeIterator(df, dict_chr2idx, 'transcript_id'))]
            rit = ReadIterator(bam_file, dict_chr2idx, reference=c, start=1, end=c_len, max_span=None, flag_filter=0) # max_span=m.max_ilen
            it = AnnotationOverlapIterator(rit, aits)
            for loc, (read, annos) in it:
                n_reads+=1
                overlapping_tids = [t[0] for (_, (_, t)) in annos[0]]
#                 if read.query_name=="ENSMUST00000119947.1_+_mat_4-5_16_37883024-37883123_NA_19,31,70,87":
#                     print(read, mapper, condition)
#                     print(annos)
#                     print(overlapping_tids)
#                     print(df)

                performance = classify_read(read, overlapping_tids, is_converted_bam, mapper, condition, dict_chr2idx, performance, out_reads, samout)
             
    # add unmapped reads
    for read in samin.fetch(contig=None, until_eof=True):
        if not read.is_unmapped:
            continue
        n_reads+=1
        performance = classify_read(read, [], is_converted_bam, mapper, condition, dict_chr2idx, performance, out_reads, samout)
    print("%s reads:  %i %i %i" % (bam_file, n_reads, samin.mapped, samin.unmapped))

    # write tid performance table
    tids= set([x for x, y, z in performance.keys()])
    isos= set([y for x, y, z in performance.keys()])
    for tid in tids:
        for iso in isos:
             print("\t".join([str(x) for x in [
                 1 if is_converted_bam else 0, 
                 mapper, 
                 condition.id,
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

    

def evaluate_splice_sites(bam_file, bam_out, is_converted_bam, m, mapper, condition, out_performance):
    """ evaluate the passed BAM file """
    logging.info("Evaluating splice sites in %s" % bam_file)
    performance = Counter()

    chromosomes = m.genome.references
    dict_chr2idx = {k: v for v, k in enumerate(chromosomes)}

    n_reads = 0
    samin = pysam.AlignmentFile(bam_file, "rb")
    samout = pysam.AlignmentFile(bam_out + '.tmp.bam', "wb", template=samin) if bam_out is not None else None

    transcripts = m.transcripts

    performance = Counter()

    for tid in transcripts:

        t = transcripts[tid]

        #print(tid + "\t" + str(len(t.introns)), file = ftxid)

        for i in range(0, len(t.introns)):

            # readBuffer = list()

            txid = t.introns.iloc[i]['transcript_id']
            exonnumber = t.introns.iloc[i]['exon_number']
            chromosome = t.introns.iloc[i]['Chromosome']
            intronstart = t.introns.iloc[i]['Start']
            intronend = t.introns.iloc[i]['End']

            intronID = txid + "_" + str(exonnumber)

            iv = Interval(intronstart - 1, intronend + 1)
            ivTree = IntervalTree()
            ivTree.add(iv)

            performance[tid, intronID, 'acceptor_rc_splicing_TP'] += 0
            performance[tid, intronID, 'acceptor_rc_splicing_FP'] += 0
            performance[tid, intronID, 'acceptor_rc_overlapping_TP'] += 0
            performance[tid, intronID, 'acceptor_rc_overlapping_FP'] += 0
            performance[tid, intronID, 'acceptor_rc_splicing_wrong'] += 0

            performance[tid, intronID, 'donor_rc_splicing_TP'] = 0
            performance[tid, intronID, 'donor_rc_splicing_FP'] = 0
            performance[tid, intronID, 'donor_rc_overlapping_TP'] = 0
            performance[tid, intronID, 'donor_rc_overlapping_FP'] = 0
            performance[tid, intronID, 'donor_rc_splicing_wrong'] = 0

            df = m.df[(m.df.Feature == 'transcript') & (m.df.Chromosome == chromosome) & (m.df.End >= intronstart - 1) & (m.df.Start <= intronend + 1)]  # get annotations
            annots = df['transcript_id'].tolist()

            # Donor
            donorIterator = SimulatedReadIterator(samin.fetch(contig=chromosome, start=intronstart - 1, stop=intronstart), flagFilter=0)

            for read in donorIterator:

                # if read.name in readBuffer:
                #     continue

                start = read.start
                end = read.end

                bamRead = read.bamRead

                color = ""

                if start <= iv.begin:

                    if read.splicing and not read.hasSpliceSite(ivTree) and end > iv.begin:

                        performance[tid, intronID, 'donor_rc_splicing_wrong'] += 1
                        bamRead.setTag('YC', colors["exonexonfalse"])
                        if (samout):
                            samout.write(bamRead)

                    elif read.splicing and read.hasSpliceSite(ivTree) and end > iv.begin:

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
                # readBuffer.append(read.name)

            # Acceptor
            acceptorIterator = SimulatedReadIterator(samin.fetch(contig=chromosome, start=intronend, stop=intronend + 1), flagFilter=0)

            for read in acceptorIterator:

                # if read.name in readBuffer:
                #     continue

                start = read.start
                end = read.end

                bamRead = read.bamRead

                if end >= iv.end:

                    if read.splicing and not read.hasSpliceSite(ivTree) and start < iv.end:

                        performance[tid, intronID, 'acceptor_rc_splicing_wrong'] += 1
                        bamRead.setTag('YC', colors["exonexonfalse"])
                        if (samout):
                            samout.write(bamRead)

                    elif read.splicing and read.hasSpliceSite(ivTree) and start < iv.end:

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
                # readBuffer.append(read.name)

    #ftxid.close()
    print("%s reads:  %i %i %i" % (bam_file, n_reads, samin.mapped, samin.unmapped))
    samin.close()

    #tids, intron, classes = performance.keys()
    #tids = set([x for x, y, z in performance.keys()])
    #introns = set([y for x, y, z in performance.keys()])
    prevIntron = ""
    prevTid = ""
    intronBuffer = dict()

    for tid, intron, classification in performance.keys():

        if intron != prevIntron and prevIntron != "":

            print("\t".join([str(x) for x in [
                1 if is_converted_bam else 0,
                mapper,
                condition,
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

            intronBuffer = dict()

        intronBuffer[classification] = performance[tid, intron, classification]
        prevTid = tid
        prevIntron = intron

    print("\t".join([str(x) for x in [
        1 if is_converted_bam else 0,
        mapper,
        condition,
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
        sambambasortbam(bam_out + '.tmp.bam',
                        bam_out,
                        index=True,
                        override=True,
                        delinFile=True)


def get_spliced_fraction(bam, tid, iv, chromosome, start, end, annots, donor = True) :

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


def calculate_splice_site_mappability(config, bam_file, truth_file, is_converted_bam, m, mapper, condition, out_mappability, out_bedGraph, fasta, readLength):
    """ evaluate the passed BAM file """
    logging.info("Calculating mappability for %s" % bam_file)

    genomeMappabilityFile = config['genome_mappability']

    chromosomes = m.genome.references
    dict_chr2idx = {k: v for v, k in enumerate(chromosomes)}

    truthBam = pysam.AlignmentFile(truth_file, "rb")
    mappedBam = pysam.AlignmentFile(bam_file, "rb")

    transcripts = m.transcripts

    tmpBed = out_bedGraph + ".tmp"
    f = open(tmpBed, "w")

    for tid in transcripts:

        t = transcripts[tid]

        for i in range(0, len(t.introns)):

            txid = t.introns.iloc[i]['transcript_id']
            exonnumber = t.introns.iloc[i]['exon_number']
            chromosome = t.introns.iloc[i]['Chromosome']
            intronstart = t.introns.iloc[i]['Start']
            intronend = t.introns.iloc[i]['End']
            intronstrand = t.introns.iloc[i]['Strand']

            intronID = txid + "_" + str(exonnumber)

            iv = Interval(intronstart - 1, intronend + 1)
            ivTree = IntervalTree()
            ivTree.add(iv)

            df = m.df[(m.df.Feature == 'transcript') & (m.df.Chromosome == chromosome) & (m.df.End >= intronstart - 1) & (m.df.Start <= intronend + 1)]  # get annotations
            annots = df['transcript_id'].tolist()

            # Donor
            donorTruthSplicedTP, donorTruthUnsplicedTP, donorTruthSplicedFP, donorTruthUnsplicedFP = get_spliced_fraction(truthBam, txid, iv, chromosome, intronstart - 1, intronstart, annots, True)
            donorMappedSplicedTP, donorMappedUnsplicedTP, donorMappedSplicedFP, donorMappedUnsplicedFP = get_spliced_fraction(mappedBam, txid, iv, chromosome, intronstart - 1, intronstart, annots, True)

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
            acceptorTruthSplicedTP, acceptorTruthUnsplicedTP, acceptorTruthSplicedFP, acceptorTruthUnsplicedFP = get_spliced_fraction(truthBam, txid, iv, chromosome, intronend, intronend + 1, annots, False)
            acceptorMappedSplicedTP, acceptorMappedUnsplicedTP, acceptorMappedSplicedFP, acceptorMappedUnsplicedFP = get_spliced_fraction(mappedBam, txid, iv, chromosome, intronend, intronend + 1, annots, False)

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
            donorExonGenomeMappability = np.zeros(readLength, dtype=float)
            donorExonBaseContent = ""
            
            it = BlockedPerPositionIterator([
                LocationTabixPerPositionIterator(
                    genomeMappabilityFile,
                    dict_chr2idx, reference=chromosome, start=intronstart - readLength,
                    end=intronstart - 1, chr_prefix_operator='add', fill=True, fill_value=(0,),
                    is_zero_based=True),
                LocationFastaIterator(fasta, dict_chr2idx, chromosome, intronstart - readLength, intronstart - 1)
                ])
            
            arrayIndex = 0
            for loc, (genomeMappability, referenceBase) in it:
                donorExonGenomeMappability[arrayIndex] = genomeMappability[0]
                donorExonBaseContent += referenceBase
                arrayIndex += 1
            
            donorIntronGenomeMappability = np.zeros(readLength, dtype=float)
            donorIntronBaseContent = ""
            
            it = BlockedPerPositionIterator([
                LocationTabixPerPositionIterator(
                    genomeMappabilityFile,
                    dict_chr2idx, reference=chromosome, start=intronstart,
                    end=intronstart + readLength - 1, chr_prefix_operator='add', fill=True, fill_value=(0,),
                    is_zero_based=True),
                LocationFastaIterator(fasta, dict_chr2idx, chromosome, intronstart, intronstart + readLength - 1)
            ])
            
            arrayIndex = 0
            for loc, (genomeMappability, referenceBase) in it:
                donorIntronGenomeMappability[arrayIndex] = genomeMappability[0]
                donorIntronBaseContent += referenceBase
                arrayIndex += 1
            
            acceptorIntronGenomeMappability = np.zeros(readLength, dtype=float)
            acceptorIntronBaseContent = ""
            
            it = BlockedPerPositionIterator([
                LocationTabixPerPositionIterator(
                    genomeMappabilityFile,
                    dict_chr2idx, reference=chromosome, start=intronend - readLength + 1,
                    end=intronend, chr_prefix_operator='add', fill=True, fill_value=(0,),
                    is_zero_based=True),
                LocationFastaIterator(fasta, dict_chr2idx, chromosome, intronend - readLength + 1, intronend)
            ])
            
            arrayIndex = 0
            for loc, (genomeMappability, referenceBase) in it:
                acceptorIntronGenomeMappability[arrayIndex] = genomeMappability[0]
                acceptorIntronBaseContent += referenceBase
                arrayIndex += 1
            
            acceptorExonGenomeMappability = np.zeros(readLength, dtype=float)
            acceptorExonBaseContent = ""
            
            it = BlockedPerPositionIterator([
                LocationTabixPerPositionIterator(
                    genomeMappabilityFile,
                    dict_chr2idx, reference=chromosome, start=intronend + 1,
                    end=intronend + readLength, chr_prefix_operator='add', fill=True, fill_value=(0,),
                    is_zero_based=True),
                LocationFastaIterator(fasta, dict_chr2idx, chromosome, intronend + 1, intronend + readLength)
            ])
            
            arrayIndex = 0
            for loc, (genomeMappability, referenceBase) in it:
                acceptorExonGenomeMappability[arrayIndex] = genomeMappability[0]
                acceptorExonBaseContent += referenceBase
                arrayIndex += 1

            if intronstrand == "-":
                print("\t".join([chromosome,str(intronstart - readLength), str(intronstart + readLength), intronID + "_acceptor", str(donorMappabilityTP), intronstrand]), file = f)
                print("\t".join([chromosome, str(intronend - readLength), str(intronend + readLength), intronID + "_donor", str(acceptorMappabilityTP), intronstrand]), file = f)

                print("\t".join([str(x) for x in [
                    1 if is_converted_bam else 0,
                    mapper,
                    condition,
                    tid,
                    intronID,
                    chromosome,
                    str(intronstart),
                    str(intronend),
                    intronstrand,
                    str(acceptorMappabilityTP),
                    str(acceptorMappabilityFP),
                    str(acceptorMappedSplicedTP + acceptorMappedUnsplicedTP),
                    str(acceptorMappedSplicedFP + acceptorMappedUnsplicedFP),
                    str(np.nanmean(acceptorExonGenomeMappability)),
                    str(acceptorExonBaseContent.upper().count("A")),
                    str(acceptorExonBaseContent.upper().count("C")),
                    str(acceptorExonBaseContent.upper().count("T")),
                    str(acceptorExonBaseContent.upper().count("G")),
                    str(np.nanmean(acceptorIntronGenomeMappability)),
                    str(acceptorIntronBaseContent.upper().count("A")),
                    str(acceptorIntronBaseContent.upper().count("C")),
                    str(acceptorIntronBaseContent.upper().count("T")),
                    str(acceptorIntronBaseContent.upper().count("G")),
                    str(donorMappabilityTP),
                    str(donorMappabilityFP),
                    str(donorMappedSplicedTP + donorMappedUnsplicedTP),
                    str(donorMappedSplicedFP + donorMappedUnsplicedFP),
                    str(np.nanmean(donorExonGenomeMappability)),
                    str(donorExonBaseContent.upper().count("A")),
                    str(donorExonBaseContent.upper().count("C")),
                    str(donorExonBaseContent.upper().count("T")),
                    str(donorExonBaseContent.upper().count("G")),
                    str(np.nanmean(donorIntronGenomeMappability)),
                    str(donorIntronBaseContent.upper().count("A")),
                    str(donorIntronBaseContent.upper().count("C")),
                    str(donorIntronBaseContent.upper().count("T")),
                    str(donorIntronBaseContent.upper().count("G"))
                ]]), file=out_mappability)
            else :
                print("\t".join([chromosome, str(intronstart - readLength), str(intronstart + readLength), intronID + "_donor", str(donorMappabilityTP), intronstrand]), file = f)
                print("\t".join([chromosome, str(intronend - readLength), str(intronend + readLength), intronID + "_acceptor", str(acceptorMappabilityTP), intronstrand]), file = f)

                print("\t".join([str(x) for x in [
                    1 if is_converted_bam else 0,
                    mapper,
                    condition,
                    tid,
                    intronID,
                    chromosome,
                    str(intronstart),
                    str(intronend),
                    intronstrand,
                    str(donorMappabilityTP),
                    str(donorMappabilityFP),
                    str(donorMappedSplicedTP + donorMappedUnsplicedTP),
                    str(donorMappedSplicedFP + donorMappedUnsplicedFP),
                    str(np.nanmean(donorExonGenomeMappability)),
                    str(donorExonBaseContent.upper().count("A")),
                    str(donorExonBaseContent.upper().count("C")),
                    str(donorExonBaseContent.upper().count("T")),
                    str(donorExonBaseContent.upper().count("G")),
                    str(np.nanmean(donorIntronGenomeMappability)),
                    str(donorIntronBaseContent.upper().count("A")),
                    str(donorIntronBaseContent.upper().count("C")),
                    str(donorIntronBaseContent.upper().count("T")),
                    str(donorIntronBaseContent.upper().count("G")),
                    str(acceptorMappabilityTP),
                    str(acceptorMappabilityFP),
                    str(acceptorMappedSplicedTP + acceptorMappedUnsplicedTP),
                    str(acceptorMappedSplicedFP + acceptorMappedUnsplicedFP),
                    str(np.nanmean(acceptorExonGenomeMappability)),
                    str(acceptorExonBaseContent.upper().count("A")),
                    str(acceptorExonBaseContent.upper().count("C")),
                    str(acceptorExonBaseContent.upper().count("T")),
                    str(acceptorExonBaseContent.upper().count("G")),
                    str(np.nanmean(acceptorIntronGenomeMappability)),
                    str(acceptorIntronBaseContent.upper().count("A")),
                    str(acceptorIntronBaseContent.upper().count("C")),
                    str(acceptorIntronBaseContent.upper().count("T")),
                    str(acceptorIntronBaseContent.upper().count("G"))
                ]]), file=out_mappability)

    f.close()

    unmerged = pybedtools.BedTool(tmpBed).sort()
    merged = unmerged.genome_coverage(bg=True,g=config['genome_chromosome_sizes']).map(unmerged, o="max")
    #merged = unmerged.genome_coverage(bg=True, g='/Volumes/groups/ameres/Niko/ref/genomes/mm10/Mus_musculus.GRCm38.dna.primary_assembly.fa.chrom.sizes').map(unmerged, o = "max")
    #merged = unmerged.genome_coverage(bg=True, g='/groups/ameres/Niko/ref/genomes/mm10/Mus_musculus.GRCm38.dna.primary_assembly.fa.chrom.sizes').map(unmerged, o="max")

    f = open(out_bedGraph, "w")

    for line in open(merged.fn):
        chr, start, end, ovlp, score = line.strip().split()
        print("\t".join([chr, start, end, score]), file = f)

    f.close()

    os.remove(tmpBed)


def get_mapq_histogram(bam_file, is_converted_bam, m, mapper, condition, out_hist):
    dict_chr2idx, dict_idx2chr, dict_chr2len=get_chrom_dicts_from_bam(bam_file)
    hist=Counter()
    samin = pysam.AlignmentFile(bam_file, "rb") 
    for c, c_len in dict_chr2len.items():
        if c not in samin.references:
            continue # no data on this chrom 
        rit = ReadIterator(bam_file, dict_chr2idx, reference=c, start=1, end=c_len, max_span=None, flag_filter=0)
        for _, r in rit:
            hist[r.mapping_quality]+=1  
    for mq in hist:
        print('\t'.join([str(x) for x in (
            bam_file,
            1 if is_converted_bam else 0,
            mapper,
            condition.id,
            mq,
            hist[mq])]), file=out_hist)
        
def uniformityStats(observed, truth):

    if truth is None:
        truthWaviness = np.nan
        truthMean = np.nan
        truthStdev = np.nan
        truthGini = np.nan
        truthKs = np.nan
        truthChisq = np.nan
    else:
        truthWaviness = np.subtract(*np.percentile(np.cumsum(truth), [75, 25]))
        truthMean = np.mean(truth)
        truthStdev = np.std(truth)
        truthGini = gini(truth)
        truthKs = stats.kstest(truth, 'uniform')[1]
        truthChisq = stats.chisquare(truth)[1]

    if observed is None:
        observedWaviness = np.nan
        observedMean = np.nan
        observedStdev = np.nan
        observedGini = np.nan
        observedKs = np.nan
        observedChisq = np.nan
    else:
        observedWaviness = np.subtract(*np.percentile(np.cumsum(observed), [75, 25]))
        observedMean = np.mean(observed)
        observedStdev = np.std(observed)
        observedGini = gini(observed)
        observedKs = stats.kstest(observed, 'uniform')[1]
        observedChisq = stats.chisquare(observed)[1]

    if (not truth is None and not observed is None) and (np.trapz(truth)!=0.0) and (np.trapz(abs(observed - truth))!=0):
        area = np.trapz(abs(observed - truth)) / np.trapz(truth)
        rmse = mean_squared_error(observed, truth, squared=False)
        KlDiv = stats.entropy(observed, truth)
        Chisq = stats.chisquare(observed, truth)[1]
        KS = stats.kstest(observed, truth)[1]
    else :
        area = np.nan
        rmse = np.nan
        KlDiv = np.nan
        Chisq = np.nan
        KS = np.nan

    return {
        'truth' : [truthMean, truthStdev, truthWaviness, truthGini, truthKs, truthChisq],
        'observed': [observedMean, observedStdev, observedWaviness, observedGini, observedKs, observedChisq],
        'global' : [area, rmse, KlDiv, Chisq, KS],
    }

def evaluate_coverage_uniformity(bam_file, truth_file, is_converted_bam, m, mapper, condition, out_performance, out_performance10bp, out_performance100bp):
    """ evaluate the passed BAM file """
    logging.info("Evaluating coverage uniformity in %s" % bam_file)
    np.seterr(divide='ignore', invalid='ignore')

    chromosomes = m.genome.references
    dict_chr2idx = {k: v for v, k in enumerate(chromosomes)}

    transcripts = m.transcripts

    performance = Counter()

    for tid in transcripts:

        t = transcripts[tid]

        intronMappedCoverage = None
        intronTruthCoverage = None

        for i in range(0, len(t.introns)):

            txid = t.introns.iloc[i]['transcript_id']
            chromosome = t.introns.iloc[i]['Chromosome']
            start = t.introns.iloc[i]['Start']
            end = t.introns.iloc[i]['End']

            mappedCoverage = np.zeros(end - start + 1, dtype=float)
            truthCoverage = np.zeros(end - start + 1, dtype=float)

            it = BlockedPerPositionIterator([
                LocationPileupIterator(bam_file, dict_chr2idx, chromosome, start, end),
                LocationPileupIterator(truth_file, dict_chr2idx, chromosome, start, end)])

            for loc, (mapped, truth) in it:

                truthPositionPileup = 0
                mappedPositionPileup = 0

                if truth:
                    for r in truth.pileups:
                        if not r.is_del and not r.is_refskip:
                            truthPositionPileup += 1

                if mapped:
                    for r in mapped.pileups:
                        if not r.is_del and not r.is_refskip:
                            mappedPositionPileup += 1

                truthCoverage[loc.start - start] = truthPositionPileup
                mappedCoverage[loc.start - start] = mappedPositionPileup

            if intronMappedCoverage is None:
                intronMappedCoverage = mappedCoverage
            else:
                intronMappedCoverage = np.concatenate((intronMappedCoverage, mappedCoverage))

            if intronTruthCoverage is None:
                intronTruthCoverage = truthCoverage
            else:
                intronTruthCoverage = np.concatenate((intronTruthCoverage, truthCoverage))

        intronStats = uniformityStats(intronMappedCoverage, intronTruthCoverage)

        bs = 10

        if not intronMappedCoverage is None and not intronTruthCoverage is None:

            if intronMappedCoverage.size % bs != 0:
                intronMappedCoverage10bp = np.concatenate([intronMappedCoverage, np.repeat(np.nan, bs - (intronMappedCoverage.size % bs))]).reshape(-1, bs)
                intronTruthCoverage10bp = np.concatenate([intronTruthCoverage, np.repeat(np.nan, bs - (intronTruthCoverage.size % bs))]).reshape(-1, bs)
            else:
                intronMappedCoverage10bp = intronMappedCoverage.reshape(-1, bs)
                intronTruthCoverage10bp = intronTruthCoverage.reshape(-1, bs)

            intronMappedCoverage10bp = np.nanmean(intronMappedCoverage10bp, axis = 1)
            intronTruthCoverage10bp = np.nanmean(intronTruthCoverage10bp, axis = 1)

        else :

            intronMappedCoverage10bp = None
            intronTruthCoverage10bp = None

        intronStats10bp = uniformityStats(intronMappedCoverage10bp, intronTruthCoverage10bp)

        bs = 100

        if not intronMappedCoverage is None and not intronTruthCoverage is None:

            if intronMappedCoverage.size % bs != 0:
                intronMappedCoverage100bp = np.concatenate([intronMappedCoverage, np.repeat(np.nan, bs - (intronMappedCoverage.size % bs))]).reshape(-1, bs)
                intronTruthCoverage100bp = np.concatenate([intronTruthCoverage, np.repeat(np.nan, bs - (intronTruthCoverage.size % bs))]).reshape(-1, bs)
            else:
                intronMappedCoverage100bp = intronMappedCoverage.reshape(-1, bs)
                intronTruthCoverage100bp = intronTruthCoverage.reshape(-1, bs)

            intronMappedCoverage100bp = np.nanmean(intronMappedCoverage100bp, axis=1)
            intronTruthCoverage100bp = np.nanmean(intronTruthCoverage100bp, axis=1)

        else :
            intronMappedCoverage100bp = None
            intronTruthCoverage100bp = None

        intronStats100bp = uniformityStats(intronMappedCoverage100bp, intronTruthCoverage100bp)

        exonMappedCoverage = None
        exonTruthCoverage = None

        for i in range(0, len(t.exons)):

            txid = t.exons.iloc[i]['transcript_id']
            chromosome = t.exons.iloc[i]['Chromosome']
            start = t.exons.iloc[i]['Start']
            end = t.exons.iloc[i]['End']

            mappedCoverage = np.zeros(end - start + 1, dtype=float)
            truthCoverage = np.zeros(end - start + 1, dtype=float)

            it = BlockedPerPositionIterator([
                LocationPileupIterator(bam_file, dict_chr2idx, chromosome, start, end),
                LocationPileupIterator(truth_file, dict_chr2idx, chromosome, start, end)])

            for loc, (mapped, truth) in it:

                truthPositionPileup = 0
                mappedPositionPileup = 0

                if truth:
                    for r in truth.pileups:
                        if not r.is_del and not r.is_refskip:
                            truthPositionPileup += 1

                if mapped:
                    for r in mapped.pileups:
                        if not r.is_del and not r.is_refskip:
                            mappedPositionPileup += 1

                truthCoverage[loc.start - start] = truthPositionPileup
                mappedCoverage[loc.start - start] = mappedPositionPileup

            if exonMappedCoverage is None:
                exonMappedCoverage = mappedCoverage
            else:
                exonMappedCoverage = np.concatenate((exonMappedCoverage, mappedCoverage))

            if exonTruthCoverage is None:
                exonTruthCoverage = truthCoverage
            else:
                exonTruthCoverage = np.concatenate((exonTruthCoverage, truthCoverage))

        exonStats = uniformityStats(exonMappedCoverage, exonTruthCoverage)

        bs = 10

        if not exonMappedCoverage is None and not exonTruthCoverage is None:

            if exonMappedCoverage.size % bs != 0:
                exonMappedCoverage10bp = np.concatenate([exonMappedCoverage, np.repeat(np.nan, bs - (exonMappedCoverage.size % bs))]).reshape(-1, bs)
                exonTruthCoverage10bp = np.concatenate([exonTruthCoverage, np.repeat(np.nan, bs - (exonTruthCoverage.size % bs))]).reshape(-1, bs)
            else:
                exonMappedCoverage10bp = exonMappedCoverage.reshape(-1, bs)
                exonTruthCoverage10bp = exonTruthCoverage.reshape(-1, bs)

            exonMappedCoverage10bp = np.nanmean(exonMappedCoverage10bp, axis = 1)
            exonTruthCoverage10bp = np.nanmean(exonTruthCoverage10bp, axis = 1)

        else :
            exonMappedCoverage10bp = None
            exonTruthCoverage10bp = None

        exonStats10bp = uniformityStats(exonMappedCoverage10bp, exonTruthCoverage10bp)

        bs = 100

        if not exonMappedCoverage is None and not exonTruthCoverage is None:

            if exonMappedCoverage.size % bs != 0:
                exonMappedCoverage100bp = np.concatenate([exonMappedCoverage, np.repeat(np.nan, bs - (exonMappedCoverage.size % bs))]).reshape(-1, bs)
                exonTruthCoverage100bp = np.concatenate([exonTruthCoverage, np.repeat(np.nan, bs - (exonTruthCoverage.size % bs))]).reshape(-1, bs)
            else:
                exonMappedCoverage100bp = exonMappedCoverage.reshape(-1, bs)
                exonTruthCoverage100bp = exonTruthCoverage.reshape(-1, bs)

            exonMappedCoverage100bp = np.nanmean(exonMappedCoverage100bp, axis=1)
            exonTruthCoverage100bp = np.nanmean(exonTruthCoverage100bp, axis=1)

        else :
            exonMappedCoverage100bp = None
            exonTruthCoverage100bp = None

        exonStats100bp = uniformityStats(exonMappedCoverage100bp, exonTruthCoverage100bp)

        print("\t".join([str(x) for x in [
            1 if is_converted_bam else 0,
            mapper,
            condition,
            tid,
            exonStats['truth'][0],
            exonStats['truth'][1],
            exonStats['truth'][2],
            exonStats['truth'][3],
            exonStats['truth'][4],
            exonStats['truth'][5],
            exonStats['observed'][0],
            exonStats['observed'][1],
            exonStats['observed'][2],
            exonStats['observed'][3],
            exonStats['observed'][4],
            exonStats['observed'][5],
            exonStats['global'][0],
            exonStats['global'][1],
            exonStats['global'][2],
            exonStats['global'][3],
            exonStats['global'][4],
            intronStats['truth'][0],
            intronStats['truth'][1],
            intronStats['truth'][2],
            intronStats['truth'][3],
            intronStats['truth'][4],
            intronStats['truth'][5],
            intronStats['observed'][0],
            intronStats['observed'][1],
            intronStats['observed'][2],
            intronStats['observed'][3],
            intronStats['observed'][4],
            intronStats['observed'][5],
            intronStats['global'][0],
            intronStats['global'][1],
            intronStats['global'][2],
            intronStats['global'][3],
            intronStats['global'][4]
        ]]), file=out_performance)

        print("\t".join([str(x) for x in [
            1 if is_converted_bam else 0,
            mapper,
            condition,
            tid,
            exonStats10bp['truth'][0],
            exonStats10bp['truth'][1],
            exonStats10bp['truth'][2],
            exonStats10bp['truth'][3],
            exonStats10bp['truth'][4],
            exonStats10bp['truth'][5],
            exonStats10bp['observed'][0],
            exonStats10bp['observed'][1],
            exonStats10bp['observed'][2],
            exonStats10bp['observed'][3],
            exonStats10bp['observed'][4],
            exonStats10bp['observed'][5],
            exonStats10bp['global'][0],
            exonStats10bp['global'][1],
            exonStats10bp['global'][2],
            exonStats10bp['global'][3],
            exonStats10bp['global'][4],
            intronStats10bp['truth'][0],
            intronStats10bp['truth'][1],
            intronStats10bp['truth'][2],
            intronStats10bp['truth'][3],
            intronStats10bp['truth'][4],
            intronStats10bp['truth'][5],
            intronStats10bp['observed'][0],
            intronStats10bp['observed'][1],
            intronStats10bp['observed'][2],
            intronStats10bp['observed'][3],
            intronStats10bp['observed'][4],
            intronStats10bp['observed'][5],
            intronStats10bp['global'][0],
            intronStats10bp['global'][1],
            intronStats10bp['global'][2],
            intronStats100bp['global'][3],
            intronStats100bp['global'][4]
        ]]), file=out_performance10bp)

        print("\t".join([str(x) for x in [
            1 if is_converted_bam else 0,
            mapper,
            condition,
            tid,
            exonStats100bp['truth'][0],
            exonStats100bp['truth'][1],
            exonStats100bp['truth'][2],
            exonStats100bp['truth'][3],
            exonStats100bp['truth'][4],
            exonStats100bp['truth'][5],
            exonStats100bp['observed'][0],
            exonStats100bp['observed'][1],
            exonStats100bp['observed'][2],
            exonStats100bp['observed'][3],
            exonStats100bp['observed'][4],
            exonStats100bp['observed'][5],
            exonStats100bp['global'][0],
            exonStats100bp['global'][1],
            exonStats100bp['global'][2],
            exonStats100bp['global'][3],
            exonStats100bp['global'][4],
            intronStats100bp['truth'][0],
            intronStats100bp['truth'][1],
            intronStats100bp['truth'][2],
            intronStats100bp['truth'][3],
            intronStats100bp['truth'][4],
            intronStats100bp['truth'][5],
            intronStats100bp['observed'][0],
            intronStats100bp['observed'][1],
            intronStats100bp['observed'][2],
            intronStats100bp['observed'][3],
            intronStats100bp['observed'][4],
            intronStats100bp['observed'][5],
            intronStats100bp['global'][0],
            intronStats100bp['global'][1],
            intronStats100bp['global'][2],
            intronStats100bp['global'][3],
            intronStats100bp['global'][4]
        ]]), file=out_performance100bp)

def initializeUniformityHeaders(out) :
    print("is_converted_bam\tmapper\tcondition\ttid", end="\t", file=out)
    print("true_exon_coverage_mean\ttrue_exon_coverage_stdev\ttrue_exon_coverage_waviness\ttrue_exon_coverage_gini\ttrue_exon_coverage_uniform_ks_test\ttrue_exon_coverage_uniform_chisq_test",end='\t', file=out)
    print("mapped_exon_coverage_mean\tmapped_exon_coverage_stdev\tmapped_exon_coverage_waviness\tmapped_exon_coverage_gini\tmapped_exon_coverage_uniform_ks_test\tmapped_exon_coverage_uniform_chisq_test",end='\t', file=out)
    print("exon_AUC\texon_RSME\texon_KL_deviation\texon_chisq\texon_KS", end='\t', file=out)
    print("true_intron_coverage_mean\ttrue_intron_coverage_stdev\ttrue_intron_coverage_waviness\ttrue_intron_coverage_gini\ttrue_intron_coverage_uniform_ks_test\ttrue_intron_coverage_uniform_chisq_test",end='\t', file=out)
    print("mapped_intron_coverage_mean\tmapped_intron_coverage_stdev\tmapped_intron_coverage_waviness\tmapped_intron_coverage_gini\tmapped_intron_coverage_uniform_ks_test\tmapped_intron_coverage_uniform_chisq_test", end='\t', file=out)
    print("intron_AUC\tintron_RSME\tintron_KL_deviation\tintron_chisq\tintron_KS", file=out)
    
    
    
#bam='/Volumes/groups/ameres/Niko/projects/Ameres/splicing/splice_sim/testruns/small4/small4/sim/bam_conv/STAR/small4.5min.STAR.TC.bam'
# bam='/Volumes/groups/ameres/Niko/projects/Ameres/splicing/splice_sim/testruns/small4/small4/sim/bam_conv/HISAT2_TLA/small4.5min.HISAT2_TLA.TC.bam'
     
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
    
    overlap_cutoff = config['overlap_cutoff'] if 'overlap_cutoff' in config else m.readlen * 0.8
    logging.info("Configured overlap cutoff: %i" % overlap_cutoff)
     
    fout = outdir + 'reads.tsv'
    fout2 = outdir + 'tid_performance.tsv'
    fout3 = outdir + 'splice_site_performance.tsv'
    fout4 = outdir + 'coverage_uniformity_performance.tsv'
    fout5 = outdir + 'coverage_uniformity_performance.10bp.tsv'
    fout6 = outdir + 'coverage_uniformity_performance.100bp.tsv'
    fout7 = outdir + 'mappability.tsv'
    fout8 = outdir + 'metadata.tsv'
    fout9 = outdir + 'mapq_hist.tsv'
    with open(fout, 'w') as out:
        with open(fout2, 'w') as out2:
            with open(fout3, 'w') as out3:
                with open(fout4, 'w') as out4:
                    with open(fout5, 'w') as out5:
                        with open(fout6, 'w') as out6:
                            with open(fout7, 'w') as out7:
                                with open(fout8, 'w') as out8:
                                    with open(fout9, 'w') as out9:
                                        print("mapped_coords\ttrue_coords\tclassification\ttid\tis_converted_bam\tmapper\tcondition_id\toverlap\ttrue_tid\ttrue_strand\ttrue_isoform\ttag\ttrue_chr\tn_true_seqerr\tn_tc_pos\tis_converted_read\toverlapping_tids\tread_name", file=out)
                                        print("is_converted_bam\tmapper\tcondition_id\tiso\ttid\tTP\tFP\tFN", file=out2)
                                        print("is_converted_bam\tmapper\tcondition\ttid\tintron_id\tdonor_rc_splicing_TP\tdonor_rc_splicing_FP\tdonor_rc_overlapping_TP\tdonor_rc_overlapping_FP\tdonor_rc_splicing_wrong\tacceptor_rc_splicing_TP\tacceptor_rc_splicing_FP\tacceptor_rc_overlapping_TP\tacceptor_rc_overlapping_FP\tacceptor_rc_splicing_wrong", file=out3)
                                        print("is_converted_bam\tmapper\tcondition\ttid\tintron_id\tchromosome\tstart\tend\tstrand\tacceptor_sj_mappability_TP\tacceptor_sj_mappability_FP\tacceptor_TP_reads\tacceptor_FP_reads\tmean_acceptorExonGenomeMappability\tacceptorExonBaseContent_A\tacceptorExonBaseContent_C\tacceptorExonBaseContent_T\tacceptorExonBaseContent_G\tmean_acceptorIntronGenomeMappability\tacceptorIntronBaseContent_A\tacceptorIntronBaseContent_C\tacceptorIntronBaseContent_T\tacceptorIntronBaseContent_G\tdonor_sj_mappability_TP\tdonor_sj_mappability_FP\tdonor_sj_TP_reads\tdonor_sj_FP_reads\tmean_donorExonGenomeMappability\tdonorExonBaseContent_A\tdonorExonBaseContent_C\tdonorExonBaseContent_T\tdonorExonBaseContent_G\tmean_donorIntronGenomeMappability\tdonorIntronBaseContent_A\tdonorIntronBaseContent_C\tdonorIntronBaseContent_T\tdonorIntronBaseContent_G",file=out7)
                                        print("bam\tis_converted_bam\tmapper\tcondition_id\tcondition_tp\tcondition_cr\tcondition_cov", file=out8)
                                        print("bam\tis_converted_bam\tmapper\tcondition_id\tmapq\tcount", file=out9)
        
                                        initializeUniformityHeaders(out4)
                                        initializeUniformityHeaders(out5)
                                        initializeUniformityHeaders(out6)
        
                                        for cond in m.conditions:
        
                                            bam_ori_truth = simdir + "bam_ori/TRUTH/" + config['dataset_name'] + "." + cond.id + ".simulated.bam"
                                            print('\t'.join([str(x) for x in [basename(bam_ori_truth), 0, 'NA', cond.id, cond.timepoint, cond.conversion_rate, cond.coverage ]]), file=out8)
                                            
                                            bam_conv_truth = simdir + "bam_conv/TRUTH/" + config['dataset_name'] + "." + cond.id + ".simulated+conversions.bam"
                                            print('\t'.join([str(x) for x in [basename(bam_conv_truth), 1, 'NA', cond.id, cond.timepoint, cond.conversion_rate, cond.coverage ]]), file=out8)
        
                                            for mapper in config['mappers'].keys():
                                                if mapper not in ["STAR", "HISAT-3N", "HISAT2_TLA", "MERANGS"]:
                                                    print("skipping unknown mapper config %s" % mapper)
                                                    continue
                                                bam_ori = simdir + "bam_ori/" + mapper + "/" + config['dataset_name'] + "." + cond.id + "." + mapper + ".nodup.bam"
                                                print('\t'.join([str(x) for x in [basename(bam_ori), 0, mapper, cond.id, cond.timepoint, cond.conversion_rate, cond.coverage ]]), file=out8)
                                                
                                                bam_conv = simdir + "bam_conv/" + mapper + "/" + config['dataset_name'] + "." + cond.id + ".conv." + mapper + ".nodup.bam"
                                                print('\t'.join([str(x) for x in [basename(bam_conv), 1, mapper, cond.id, cond.timepoint, cond.conversion_rate, cond.coverage ]]), file=out8)
                                                
                                                bam_out_dir_ori = outdir + "bam_ori/" + mapper + "/"
                                                if not os.path.exists(bam_out_dir_ori):
                                                    print("Creating dir " + bam_out_dir_ori)
                                                    os.makedirs(bam_out_dir_ori)
                                                bam_out_dir_conv = outdir + "bam_conv/" + mapper + "/"
                                                if not os.path.exists(bam_out_dir_conv):
                                                    print("Creating dir " + bam_out_dir_conv)
                                                    os.makedirs(bam_out_dir_conv)
                                                bam_ori_out = bam_out_dir_ori + config['dataset_name'] + "." + cond.id + "." + mapper + ".mismapped.bam"
                                                bam_ori_out_intron = bam_out_dir_ori + config['dataset_name'] + "." + cond.id + "." + mapper + ".intron.bam"
                                                mappability_tid_ori_out = bam_out_dir_ori + config['dataset_name'] + "." + cond.id + "." + mapper + ".tid_mappability.bedGraph"
                                                mappability_SJ_ori_out = bam_out_dir_ori + config['dataset_name'] + "." + cond.id + "." + mapper + ".SJ_mappability.bedGraph"
                                                bam_conv_out =  bam_out_dir_conv + config['dataset_name'] + "." + cond.id + "." + mapper + ".conv.mismapped.bam"
                                                bam_conv_out_intron = bam_out_dir_conv + config['dataset_name'] + "." + cond.id + "." + mapper + ".conv.intron.bam"
                                                mappability_tid_conv_out = bam_out_dir_conv + config['dataset_name'] + "." + cond.id + "." + mapper + ".conv.tid_mappability.bedGraph"
                                                mappability_SJ_conv_out = bam_out_dir_conv + config['dataset_name'] + "." + cond.id + "." + mapper + ".conv.SJ_mappability.bedGraph"
        
                                                ##evaluate_bam(bam_ori, bam_ori_out, False, m, mapper, cond, out, out2)
                                                ##evaluate_splice_sites(bam_ori, bam_ori_out_intron, False, m, mapper, cond.id, out3)
                                                # evaluate_coverage_uniformity(bam_ori, bam_ori_truth, False, m, mapper, cond.id, out4, out5, out6)
                                                calculate_splice_site_mappability(config, bam_ori, bam_ori_truth, False, m, mapper, cond.id, out7, mappability_SJ_ori_out, config["genome_fa"], config["readlen"])
                                                get_mapq_histogram(bam_ori, False, m, mapper, cond, out9)
                                                
                                                ##evaluate_bam(bam_conv, bam_conv_out, True, m, mapper, cond, out, out2)
                                                ##evaluate_splice_sites(bam_conv, bam_conv_out_intron, True, m, mapper, cond.id, out3)
                                                # evaluate_coverage_uniformity(bam_conv, bam_conv_truth, True, m, mapper, cond.id, out4, out5, out6)
                                                calculate_splice_site_mappability(config, bam_conv, bam_conv_truth, True, m, mapper, cond.id, out7, mappability_SJ_conv_out, config["genome_fa"], config["readlen"])
                                                get_mapq_histogram(bam_conv, True, m, mapper, cond, out9)
        
                                            # eval truth
                                            bam_out_dir_ori = outdir + "bam_ori/TRUTH/"
                                            if not os.path.exists(bam_out_dir_ori):
                                                print("Creating dir " + bam_out_dir_ori)
                                                os.makedirs(bam_out_dir_ori)
                                            bam_out_dir_conv = outdir + "bam_conv/TRUTH/"
                                            if not os.path.exists(bam_out_dir_conv):
                                                print("Creating dir " + bam_out_dir_conv)
                                                os.makedirs(bam_out_dir_conv)
        
                                            bam_ori_truth_intron = bam_out_dir_ori + config['dataset_name'] + "." + cond.id + ".TRUTH.intron.bam"
                                            bam_conv_truth_intron = bam_out_dir_conv + config['dataset_name'] + "." + cond.id + ".TRUTH.TC.intron.bam"
        
                                            evaluate_bam(bam_ori_truth, None, False, m, 'NA', cond, out, out2)
                                            evaluate_splice_sites(bam_ori_truth, bam_ori_truth_intron, False, m, 'NA', cond.id, out3)
                                            evaluate_bam(bam_conv_truth, None, True, m, 'NA', cond, out, out2)
                                            evaluate_splice_sites(bam_conv_truth, bam_conv_truth_intron, True, m, 'NA', cond.id, out3)

    bgzip(fout, delinFile=True, override=True)
    logging.info("All done in %s" % str(datetime.timedelta(seconds=time.time() - startTime)))
    print("All done in", datetime.timedelta(seconds=time.time() - startTime))
    
