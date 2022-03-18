'''
@author: niko.popitsch@imba.oeaw.ac.at
'''
from collections import *
import csv, datetime, time, logging, sys, os, json, pickle
import pysam
import numpy as np
from genomic_iterators import Location, get_chrom_dicts_from_bam, get_chrom_dicts
from genomic_iterators.bam_iterators import ReadIterator
from genomic_iterators.grouped_iterators import BlockedIterator, BlockLocationIterator, AnnotationOverlapIterator
from splice_sim.utils import *
from splice_sim.model import *
from functools import reduce

read_colors={
    'fn': '200,200,200',
    'fn+fp': '200,200,100',
    'fp+fn': '200,100,0',
    'fp': '200,0,0'
    }


def merge_bam_files(out_file, bam_files, sort_output=False, del_in_files=False):
    """ merge multiple BAM files and sort + index results """
    if bam_files is None or len(bam_files) == 0:
        print("no input BAM file provided")
        return None
    samfile = pysam.AlignmentFile(bam_files[0], "rb")  # @UndefinedVariable
    out = pysam.AlignmentFile(out_file + '.unsorted.bam', "wb", template=samfile)  # @UndefinedVariable
    for f in bam_files:
        samfile = None
        try:
            samfile = pysam.AlignmentFile(f, "rb")  # @UndefinedVariable
            for read in samfile.fetch(until_eof=True):
                out.write(read)  
        except Exception as e:
            print("error opening bam %s: %s" % (f, e))            
        finally:
            if samfile:
                samfile.close()
    out.close()
    if sort_output:
        try:
            pysam.sort("-o", out_file, out_file + '.unsorted.bam')  # @UndefinedVariable
            os.remove(out_file + '.unsorted.bam')
            if del_in_files:
                for f in bam_files + [b + '.bai' for b in bam_files]:
                    if os.path.exists(f):
                        os.remove(f)
        except Exception as e:
            print("error sorting this bam: %s" % e)
    else:
        os.rename(out_file + '.unsorted.bam', out_file)
        if del_in_files:
            for f in bam_files + [b + '.bai' for b in bam_files]:
                os.remove(f)
    # index
    try:
        pysam.index(out_file)  # @UndefinedVariable
    except Exception as e:
        print("error indexing bam: %s" % e)
    return out_file

def read_aligns_to_loc(loc, read):
    """ Tests whether a read aligns to the passed location """
    # chr is not checked
    for read_start,read_end in read.get_blocks():
        if loc.start <= read_end and read_start < loc.end:
            return True
    return False

def calc_overlap(a, b):
    """ overlap between two intervals """
    return max(0, min(a[1], b[1]) - max(a[0], b[0]))

def read_envelops_loc(iloc, read):
    """ Tests whether a spliced read covers the passed (intronic) location """
    return (read.reference_start < iloc.start) and (read.reference_end>=iloc.end) and ('N' in read.cigarstring)

def read_splices_intron(iloc, read, max_diff=5):
    """ Tests whether the passed read splices the passed intronic) location. 
    The maximum difference in coordinates between the intron and the splicing block (N section) can be configured.
    A value>0 allows for some fuzziness """
    os,oe=None,None
    for [start,end] in read.get_blocks():
        if os is not None:
            splice_block=(oe+1,start-1)
            diff=abs(splice_block[0]-iloc.start)+abs(splice_block[1]-iloc.end)
            if diff <= max_diff:
                return True
        os,oe=start,end
    return False

def get_class(read, features):
    """ Calculate the classification of the passed read wrt different classification groups """
    classification={
        'tx' : set(),
        'fx' : set(),
        'don':set(),
        'acc':set(),
        'spl':set()}
    for feat_loc,(tid,ftype,fid) in features: # NB: chrom is not checked, must be matched beforehand!
        strand_match = feat_loc.strand=='-' if read.is_reverse else feat_loc.strand=='+' 
        if not strand_match:
            continue
        # does read align into feature?
        if read_aligns_to_loc(feat_loc, read): # NB no need to check chr here
            classification['tx'].add(tid) 
            classification['fx'].add(fid) 
            if ftype=='intron':
                if read.reference_start < feat_loc.start and read.reference_end >= feat_loc.start:
                    classification['don' if feat_loc.strand=='+' else 'acc'].add(fid)
                elif read.reference_start < feat_loc.end and read.reference_end >= feat_loc.end:
                    classification['acc' if feat_loc.strand=='+' else 'don'].add(fid)
        # does read span over feature?
        if read_envelops_loc(feat_loc, read) and ftype=='intron' and read_splices_intron(feat_loc, read):
            classification['spl'].add(fid)
    return classification

def set_compare(truth, truth_ext, found):
    """ Compares two sets and returns TP, FP, FNs.
    truth contains the classification of the true read wrt. features of the true tx, 
    found contains the classification of the found read wrt. overlapping features and
    truth_ext contains the classification of the true read wrt. the union of features overlapping the found read and true tx features. """
    TP = truth & found
    FP = found.difference(truth_ext)
    FN = truth.difference(found)
    return TP, FP, FN

classtype2cat={
    'tx' : 'tx',
    'fx' : 'fx',
    'don': 'sj',
    'acc': 'sj',
    'spl': 'sj'}
empty_class={'tx' : set(),'fx' : set(),'don': set(),'acc': set(),'spl': set()}
 
def classify_read(found_read, found_overlapping_raw, tid2feat, performance, dict_samout, dict_chr2idx, min_mapq=20):
    # recreate true read
    true_tid, true_strand, true_isoform, read_tag, true_chr, true_start, true_cigar, n_seqerr, n_converted, is_converted_read = found_read.query_name.split('_')
    true_read = pysam.AlignedSegment()
    true_read.query_name=found_read.query_name
    true_read.cigarstring=true_cigar
    true_read.reference_start=int(true_start)
    true_read.is_reverse = true_strand=='-'
    cv1=1 if int(n_converted)>0 else 0
    cv2=1 if int(n_converted)>1 else 0
    se1=1 if int(n_seqerr)>0 else 0
    se2=1 if int(n_seqerr)>1 else 0
    # which features does the true read overlap with and whats the true iclass   
    true_features=set(tid2feat[true_tid]) 
    true_class=get_class(true_read, true_features)
    true_features_ext=set(tid2feat[true_tid] + found_overlapping_raw)
    true_class_ext=true_class if true_features_ext==true_features else get_class(true_read, true_features_ext)
    found_class={}
    found_class[0]=empty_class if found_read.is_unmapped else get_class(found_read, found_overlapping_raw) 
    found_class[1]=empty_class if ((found_read.is_unmapped) or (found_read.mapping_quality<min_mapq)) else found_class[0]
    # calculate TP,FP,FN (+/- mq filtering)
    for mq_fil in [0,1]:
        for class_type in true_class.keys():
            TP, FP, FN = set_compare(true_class[class_type], true_class_ext[class_type], found_class[mq_fil][class_type])
            for fid in TP:
                performance[class_type, mq_fil, fid, true_isoform, cv1, cv2, se1, se2, 'TP']+=1.0
            for fid in FN:
                performance[class_type, mq_fil, fid, true_isoform, cv1, cv2, se1, se2, 'FN']+=1.0
            for fid in FP:
                performance[class_type, mq_fil, fid, true_isoform, cv1, cv2, se1, se2, 'FP']+=1.0/len(FP)
                performance[class_type, mq_fil, fid, true_isoform, cv1, cv2, se1, se2, 'FP_raw']+=1.0
                
            # if class_type=='spl' and mq_fil==0 and found_read.query_name=='ENSMUST00000162558.7_+_mat_4-92_1_161036013_1=193N37=304N31=240N31=_0_0_1':
            #     print(true_class[class_type])  
            #     print(true_class_ext[class_type])  
            #     print(found_class[mq_fil][class_type])  
            #     print(TP, FP, FN)
            #     print(found_overlapping_raw)
            #     print(found_read.get_blocks())
            #     for feat_loc, (tid, ftype, fid) in found_overlapping_raw:
            #         if ftype=='intron':
            #             print(fid, read_envelops_loc(feat_loc, found_read), read_splices_intron(feat_loc, found_read), feat_loc)
                    
                    
                
            # write reads (no mq filtering!)  if mismapped
            if (mq_fil==0) and (dict_samout is not None) and (len(FN)+len(FP)>0):
                samout=dict_samout[classtype2cat[class_type]]
                # write FN
                if len(FN)>0:
                    is_correct_strand = ( found_read.is_reverse and true_strand=='-' ) or ((not found_read.is_reverse) and true_strand=='+')
                    true_read.query_sequence = found_read.query_sequence if is_correct_strand else reverse_complement(found_read.query_sequence)
                    true_read.query_qualities = found_read.query_qualities
                    true_read.reference_id = dict_chr2idx[true_chr]
                    true_read.tags = (("NM", n_seqerr + n_converted ),)
                    true_read.flag = 0 if true_strand=='+' else 16
                    true_read.set_tag(tag='YC', value=read_colors['fn+fp'] if len(FP)>0 else read_colors['fn'], value_type="Z")
                    samout.write(true_read)
                # write FP
                if len(FP)>0:
                    found_read.set_tag(tag='YC', value=read_colors['fp+fn'] if len(FN)>0 else read_colors['fp'], value_type="Z")
                    samout.write(found_read)
    return performance


def classify_region(target_region, config, bam_file, chrom2feat, tid2feat):
    """ Classify all reads in the passed target region """
    min_mapq=config['min_mapq'] if 'min_mapq' in config else 20
    write_bam=config['write_bam'] if 'write_bam' in config else True
    samin = pysam.AlignmentFile(bam_file, "rb")
    dict_chr2idx, dict_idx2chr, dict_chr2len=get_chrom_dicts_from_bam(bam_file)
    dict_samout = {
        'tx': pysam.AlignmentFile(out_file_prefix+'.tx.eval.%s.bam' % str(target_region), "wb", template=samin),
        'fx': pysam.AlignmentFile(out_file_prefix+'.fx.eval.%s.bam' % str(target_region), "wb", template=samin),
        'sj': pysam.AlignmentFile(out_file_prefix+'.sj.eval.%s.bam' % str(target_region), "wb", template=samin)
        } if write_bam else None
    performance=Counter()
    if target_region=='unmapped':
        for unmapped_read in samin.fetch(contig=None, until_eof=True):
            if not unmapped_read.is_unmapped:
                continue
            performance=classify_read(unmapped_read, [], tid2feat, performance, dict_samout, dict_chr2idx, min_mapq)        
    else:
        chrom = dict_idx2chr[target_region.chr_idx]
        if chrom in samin.references: # there are reads on this chrom
            return performance, dict_samout
            aits = [BlockLocationIterator(iter(sorted(chrom2feat[chrom] if chrom in chrom2feat else [], key=lambda x: x[0]) ))]
            rit = ReadIterator(samin, dict_chr2idx, reference=chrom, start=target_region.start, end=target_region.end, max_span=None, flag_filter=0) # max_span=m.max_ilen
            it = AnnotationOverlapIterator(rit, aits, check_read_alignment=False)
            for loc, (read, annos) in it:
                wr+=1
                overlapping_annos = annos[0] # select 1st entry as theree is only one ait
                performance=classify_read(read, overlapping_annos, tid2feat, performance, dict_samout, dict_chr2idx, min_mapq)
    # close SAM files    
    samin.close()   
    if dict_samout is not None:
        for cat,samout in dict_samout.items():
            samout.close() 
            dict_samout[cat]=out_file_prefix+'.%s.eval.%s.bam' % str(cat,target_region)
    return performance, dict_samout, wr

# evaluate
# --bam_file /Volumes/groups/ameres/Niko/projects/Ameres/splicing/splice_sim/testruns/big3_slamseq_nf/sim/final_bams/big3_slamseq_nf.cr0.1.STAR.final.bam
# --config /Volumes/groups/ameres/Niko/projects/Ameres/splicing/splice_sim/testruns/big3_slamseq_nf/eva/debug/config.json
# --model /Volumes/groups/ameres/Niko/projects/Ameres/splicing/splice_sim/testruns/big3_slamseq_nf/sim/reference_model/big3_slamseq_nf.model
# --outdir /Volumes/groups/ameres/Niko/projects/Ameres/splicing/splice_sim/testruns/big3_slamseq_nf/eva/debug/

def evaluate(config, m, bam_file, out_dir): 
    """ evaluate the passed BAM file """
    logging.info("Evaluating %s" % bam_file)
    assert os.path.exists(bam_file) and os.path.exists(bam_file+'.bai'), "Could not find bam file (or idx) for %s" % (bam_file)
    # parse cr/mapper from file name    
    assert '.cr' in bam_file, "Cannot parse bam file name %s" % (bam_file)
    tmp = bam_file[:-10] if bam_file.endswith('.final.bam') else bam_file[:-4] # file ends with .bam or .final.bam
    cr,mapper=tmp[tmp.find('.cr')+3:].rsplit('.', 1)
    out_file_prefix=out_dir+'/'+config['dataset_name']+'.cr'+cr+'.'+mapper
    cr=float(cr)
    # get config params
    dict_chr2idx, dict_idx2chr, dict_chr2len=get_chrom_dicts_from_bam(bam_file)
    # target regions
    if 'target_region' in config:
        target_regions=[Location.from_str(config['debug_region'], dict_chr2idx)]
    else:
        target_regions=[Location(dict_chr2idx[c], 1, c_len) for c, c_len in dict_chr2len.items()] + ['unmapped'] # unmapped reads    
    # collect= features for chrom
    chrom2feat=OrderedDict()
    tid2feat=OrderedDict()
    for tid,t in m.transcripts.items():
        c=t.chromosome
        if c not in chrom2feat:
            chrom2feat[c]=[]
        if tid not in tid2feat:
            tid2feat[tid]=[]
        exon_len=0
        for exon in t.exons:
            fid = tid + "_ex" + str(exon['exon_number'])
            loc=Location(dict_chr2idx[t.chromosome], exon['start'], exon['end'], strand=t.strand)
            chrom2feat[c].append((loc,(tid, 'exon', fid)))
            tid2feat[tid].append((loc,(tid, 'exon', fid)))
            exon_len+=len(loc)
        intron_len=0
        for intron in t.introns:
            loc=Location(dict_chr2idx[t.chromosome], intron['start'], intron['end'], strand=t.strand)
            fid = tid + "_in" + str(intron['exon_number']-1 if t.strand=='+' else intron['exon_number'])
            chrom2feat[c].append((loc,(tid, 'intron', fid)))
            tid2feat[tid].append((loc,(tid, 'intron', fid)))
            intron_len+=len(loc)
    # run parallel
    performance=Counter()
    dict_samout=None
    wr=0
    with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
        future_param = {executor.submit(classify_region, reg, bam_file, chrom2feat, tid2feat): (reg) for reg in target_regions }
        pbar = tqdm.tqdm(total=len(target_regions), desc="Annotated regions")
        for future in concurrent.futures.as_completed(future_param):
            para = future_param[future]
            region_performance, region_dict_samout, region_wr = future.result()
            wr+=region_wr
            performance.update(region_performance)
            if region_dict_samout is not None:
                if dict_samout is None:
                    dict_samout={'tx':[], 'fx':[], 'sj':[]}
                for k in dict_samout.keys():
                    dict_samout[k]+=region_dict_samout[k]
            pbar.set_description('%s: wrote %i reads' % (c, wr))
            pbar.update(1)
        pbar.close()             
    # write performance data
    with open(out_file_prefix+'.counts.tsv', 'w') as out_raw:
        with open(out_file_prefix+'.counts.mq%i.tsv' % min_mapq, 'w') as out_fil:
            header='\t'.join(str(x) for x in [
                    'mapper', 'conversion_rate', 'mq_fil', 'class_type', 'fid', 'true_isoform', 'cv1', 'cv2', 'se1', 'se2', 'classification', 'count'
                    ])
            print(header, file=out_raw)
            print(header, file=out_fil)
            for (class_type, mq_fil, fid, true_isoform, cv1, cv2, se1, se2, classification), count in performance.items():
                out=out_raw if mq_fil == 0 else out_fil
                print('\t'.join(str(x) for x in [
                    mapper, cr, mq_fil, class_type, fid, true_isoform, cv1, cv2, se1, se2, classification, count
                    ]), file=out)
    # merge, sort and index BAM files    
    if dict_samout is not None:
        for cat,bam_files in dict_samout.items():
            out_file = out_file_prefix+'.%s.eval.bam' % str(cat)
            merge_bam_files(out_file, bam_files, sort_output=True, del_in_files=True)
    print("All done. Considered %i reads" % wr)

def flatten(t):
    """ Flatten nested list """
    return [i for s in t for i in s]

def calc_mappability(genomeMappability, chrom, start, end):
    """ Calculate genomic mappability for giveen window """ 
    map = flatten([[float(x[3])]*(min(end+1,int(x[2]))-max(start,int(x[1]))) for x in genomeMappability.fetch('chr'+chrom, start, end, parser=pysam.asTuple())])
    map=map+[0]*(end-start+1-len(map)) # pad with missing zeros
    mea_map=np.mean(map) if len(map)>0 else 0 
    return mea_map


# extract_feature_metadata
# --config /Volumes/groups/ameres/Niko/projects/Ameres/splicing/splice_sim/testruns/big3_slamseq_nf/eva/debug/config.json
# --model /Volumes/groups/ameres/Niko/projects/Ameres/splicing/splice_sim/testruns/big3_slamseq_nf/sim/reference_model/big3_slamseq_nf.model
# --outdir /Volumes/groups/ameres/Niko/projects/Ameres/splicing/splice_sim/testruns/big3_slamseq_nf/eva/debug/
def extract_feature_metadata(config, m, out_dir):
    """ extract metadata fox tx/fx/sj """
    genomeMappability = pysam.TabixFile(config['genome_mappability'], mode="r")
    readlen=config['readlen']
    # calculate transcript overlap
    tid2info={}
    tiv = m.build_transcript_iv()
    with open(out_dir+'/'+config['dataset_name']+'.tx.metadata.tsv', 'w') as out_tx:
        with open(out_dir+'/'+config['dataset_name']+'.fx.metadata.tsv', 'w') as out_fx:
            with open(out_dir+'/'+config['dataset_name']+'.sj.metadata.tsv', 'w') as out_sj:
                print("""tid\tftype\trnk\tchromosome\tstart\tend\tstrand\t""" +
                      """A\tC\tT\tG\tmean_map\texon_len\tintron_len\tn_overlapping\toverlapping_tids""",
                                                  file=out_tx)
                print("""tid\tfid\tftype\trnk\tchromosome\tstart\tend\tstrand\t""" +
                      """A\tC\tT\tG\tmean_map""",
                                                  file=out_fx)
                print("""tid\tfid\tftype\trnk\tchromosome\tstart\tend\tstrand\t""" +
                      """don_ex_A\tdon_ex_C\tdon_ex_T\tdon_ex_G\tdon_in_A\tdon_in_C\tdon_in_T\tdon_in_G\tdon_win_map\t""" +
                      """acc_ex_A\tacc_ex_C\tacc_ex_T\tacc_ex_G\tacc_in_A\tacc_in_C\tacc_in_T\tacc_in_G\tacc_win_map""",
                                                  file=out_sj)
                for tid, t in m.transcripts.items():                    
                    # get RNA subseq from transcript seq. NOTE: may be shorter than 2xreadlen
                    rna_seq=reverse_complement(t.transcript_seq) if t.strand == "-" else t.transcript_seq
                    # calc mean mappability
                    tx_map = calc_mappability(genomeMappability, t.chromosome, t.start, t.end)
                    # get overlapping tids
                    overlapping=[x.data.tid for x in tiv[t.chromosome].overlap(t.start, t.end) if x.data.tid != tid]
                    l=len(rna_seq)
                    exon_len, intron_len=0,0   
                    # exons                 
                    for exon in t.exons:
                        rnk=exon['exon_number']
                        fid = tid + "_ex" + str(rnk)
                        ex_map = calc_mappability(genomeMappability, t.chromosome, exon['start'], exon['end'])
                        ex_rna = rna_seq[exon['start']-t.start:exon['end']-t.start+1]
                        exon_len+=exon['end']-exon['start']+1
                        print("\t".join([str(x) for x in [
                            tid,
                            fid,
                            'exon',
                            rnk,
                            t.chromosome,
                            exon['start'],
                            exon['end'],
                            exon['strand'],
                            ex_rna.count("A"),
                            ex_rna.count("C"),
                            ex_rna.count("T"),
                            ex_rna.count("G"),
                            ex_map
                        ]]), file=out_fx)
                    # introns
                    for intron in t.introns:
                        rnk=intron['exon_number']-1 if t.strand=='+' else intron['exon_number']
                        fid = tid + "_in" + str(rnk)
                        # tx
                        in_map = calc_mappability(genomeMappability, t.chromosome, intron['start'], intron['end'])
                        in_rna = rna_seq[intron['start']-t.start:intron['end']-t.start+1]
                        intron_len+=intron['end']-intron['start']+1
                        print("\t".join([str(x) for x in [
                            tid,
                            fid,                            
                            'intron',
                            rnk,
                            t.chromosome,
                            intron['start'],
                            intron['end'],
                            intron['strand'],
                            in_rna.count("A"),
                            in_rna.count("C"),
                            in_rna.count("T"),
                            in_rna.count("G"),
                            in_map
                        ]]), file=out_fx)
                        # SJ
                        # note that win can be shorter than readlen if we run out of sequence!
                        don_win_genomic = [intron['start'] - readlen, intron['start'] + readlen] if intron['strand'] == "+" else [intron['end'] - readlen, intron['end'] + readlen]
                        acc_win_genomic = [intron['start'] - readlen, intron['start'] + readlen] if intron['strand'] == "-" else [intron['end'] - readlen, intron['end'] + readlen]
                        don_win_map = calc_mappability(genomeMappability, t.chromosome, don_win_genomic[0], don_win_genomic[1])
                        acc_win_map = calc_mappability(genomeMappability, t.chromosome, acc_win_genomic[0], acc_win_genomic[1])
                        don_seq_rna_ex = rna_seq[max(0,t.end-intron['end']-readlen):min(t.end-intron['end'],l)] if intron['strand'] == "-" else rna_seq[max(0,intron['start']-readlen-t.start):min(l,intron['start']-t.start)]
                        don_seq_rna_in = rna_seq[max(0,t.end-intron['end']):min(t.end-intron['end']+readlen,l)] if intron['strand'] == "-" else rna_seq[max(0,intron['start']-t.start):min(l,intron['start']+readlen-t.start)]
                        acc_seq_rna_in = rna_seq[max(0,t.end-intron['start']-readlen+1):min(t.end-intron['start']+1,l)] if intron['strand'] == "-" else rna_seq[max(0,intron['end']-readlen-t.start+1):min(l,intron['end']-t.start+1)]
                        acc_seq_rna_ex = rna_seq[max(0,t.end-intron['start']+1):min(t.end-intron['start']+readlen+1,l)] if intron['strand'] == "-" else rna_seq[max(0,intron['end']-t.start+1):min(l,intron['end']+readlen-t.start+1)]
                        print("\t".join([str(x) for x in [
                            tid,
                            fid,
                            'intron',
                            rnk,
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
                            don_win_map,
                            # acceptor
                            acc_seq_rna_ex.count("A"),
                            acc_seq_rna_ex.count("C"),
                            acc_seq_rna_ex.count("T"),
                            acc_seq_rna_ex.count("G"),
                            acc_seq_rna_in.count("A"),
                            acc_seq_rna_in.count("C"),
                            acc_seq_rna_in.count("T"),
                            acc_seq_rna_in.count("G"),  
                            acc_win_map
                        ]]), file=out_sj)
                    # tx info
                    print("\t".join([str(x) for x in [
                                tid,
                                'tx',
                                len(t.exons),
                                t.chromosome,
                                t.start,
                                t.end,
                                t.strand,
                                rna_seq.count("A"),
                                rna_seq.count("C"),
                                rna_seq.count("T"),
                                rna_seq.count("G"),
                                tx_map,
                                exon_len,
                                intron_len,
                                len(overlapping),
                                ','.join(overlapping)
                            ]]), file=out_tx)
    print("All done.")

