import argparse
from collections import *
import csv, datetime, time, logging, sys, os, json
import pysam
from model import Model, Condition, Isoform, Transcript
import itertools
from config_creator import calculate_transcript_data
import random

# Table of reverse complement bases
COMP_TABLE = {
    "A": 'T', "C": 'G', "T": 'A', "G": 'C'
    }
def reverse_complement(seq):
    """ Calculate reverse complement DNA sequence """
    rev=[]
    for c in seq[::-1]:
        c=c.upper()
        rev+=[COMP_TABLE.get(c, 'N')]
    return ''.join(rev)
                        
def load_model_from_config(config_file):
    config = json.load(open(config_file), object_pairs_hook=OrderedDict)
    if 'transcripts' in config:
        logging.info("Reading transcript configuration from config file.")
    else:
        assert "transcript_data" in config, "Transcript data needs to be configured either in config file ('transcripts' section) or in an external file referenced via 'transcript_data'"
        tfile = config['transcript_data']
        if not os.path.isabs(tfile):
            tfile = config_dir+"/"+tfile
        if not os.path.isfile(tfile):
            logging.info("Creating external transcript config file %s" % tfile)
            config_dir = os.path.dirname(os.path.abspath(config_file))+'/'
            calculate_transcript_data(config, config_dir, outdir)
            
        logging.info("Reading transcript configuration from external config file %s" % tfile)
        tdata = json.load(open(tfile), object_pairs_hook=OrderedDict)
        config["transcripts"]=tdata
    m = Model( config )
    return m
                        
def transcript2genome_bam(transcript_bam_file, mod, out_bam):
    """ Reads a transcriptome-mapped BAM and converts to a genome BAM. 
    Transcript ids, strand, isoform and tags are parsed from the input read names and matched with the model """
    chromosomes = mod.genome.references
    dict_chr2idx = {k: v for v, k in enumerate(chromosomes)}
    dict_idx2chr = {v: k for v, k in enumerate(chromosomes)}
    dict_chr2len = {c: mod.genome.get_reference_length(c) for c in mod.genome.references}
    logging.info("Configured conditions: %s" % (", ".join(str(x) for x in mod.conditions)) )

    header = { 'HD': {'VN': '1.4'}, 
               'SQ': [{'LN': l, 'SN': c} for c,l in dict_chr2len.items()] }
    stats=Counter()
    with pysam.AlignmentFile(out_bam+'.tmp.bam', "wb", header=header) as bam_out:        
        sam = pysam.AlignmentFile(transcript_bam_file, "rb")
        for r in sam.fetch(until_eof=True):
            stats['input_reads']+=1
            # e.g. 'ENSMUST00000099151.5_+_mat_0-400'. tag ('0-400') consists of running number from 'abundance' and running number form ART
            tid,transcript_strand,iso_id,tag=r.query_name.split("_")
            read_strand = "-" if r.is_reverse else "+"
            stats['input_reads_'+read_strand]+=1
            if not transcript_strand == read_strand:
                stats['wrong_strand_reads']+=1
            else: # NOTE: reads that were mapped to opposite strand will be dropped later in post_filtering step!
                iso = mod.transcripts[tid].isoforms[iso_id]
                ablocks=iter(iso.aln_blocks)
                genome_cigatuples=[]
                seq=""
                qual=[]
                rpos=0
                gpos=iso.t.transcript.Start-1
                mismatches=[]
                n_deletions, n_insertions = 0,0
                # get read data
                if r.is_reverse:
                    r_pos = iso.aln_block_len - r.reference_end
                    r_cigartuples=reversed(r.cigartuples)
                    r_query_sequence=reverse_complement(r.query_sequence)
                    r_query_qualities=list(reversed(r.query_qualities))
                else:
                    r_pos = r.reference_start
                    r_cigartuples=r.cigartuples
                    r_query_sequence=r.query_sequence
                    r_query_qualities=r.query_qualities
                
                # consume r_pos bases and find 1st alignment block
                ablock=next(ablocks, None) # genome coords
                off=r_pos
                while off>=0:
                    blen=ablock[1]-ablock[0]+1
                    if gpos+off>=ablock[1]:
                        off-=blen
                        ablock=next(ablocks, None)
                        if ablock is None:
                            sys.exit("ERR1: no more ablocks")                    
                        gpos=ablock[0]-1
                    else:
                        gpos+=off
                        break
                start_pos=gpos
                for op,l in r_cigartuples:
                    if op in [0,7,8]: # M, =, X
                        inc_read=l
                        inc_pos=l
                    elif op ==1: # I
                        inc_read=l
                        inc_pos=0
                        n_insertions+=l
                    elif op ==2: # D
                        inc_read=0
                        inc_pos=l
                        n_deletions+=l
                    elif op in [4]: # S
                        inc_read=l
                        inc_pos=0
                    elif op in[3]: #  N
                        inc_read=0
                        inc_pos=l
                    elif op in [5,6]: # H, P
                        inc_read=0
                        inc_pos=0
                    else:
                        sys.exit("unsupported CIGAR op %i" % op)
                    diff=gpos+inc_pos-ablock[1]
                    while diff>=0:
                        if op==8:
                            #FIXME seq err streteches that are spliced are partially ignored
                            print('ignoring seq diff X')
                        # we need to insert a splice cigartuple
                        matchlen=l-diff
                        genome_cigatuples+=[(op, matchlen)]
                        gpos+=matchlen
                        inc_pos-=matchlen
                        if inc_read>0:
                            seq+=r_query_sequence[rpos:rpos+matchlen]
                            qual+=r_query_qualities[rpos:rpos+matchlen]
                            rpos+=matchlen
                            inc_read-=matchlen
                        # N-block
                        nblock=next(ablocks, None)
                        if nblock is None:
                            sys.exit("ERR: no more ablocks")
                        n_len = nblock[0]-ablock[1]-1
                        genome_cigatuples+=[(3, n_len)]
                        gpos+=n_len
                        # update l 
                        ablock=nblock
                        l-=matchlen
                        diff=gpos+inc_pos-ablock[1]

                    genome_cigatuples+=[(op, l)]
                    if op == 8 : # X - store mismatches
                        for i in range(l):
                            mismatches+=[gpos-start_pos+1] 
                    gpos+=inc_pos
                    if inc_read>0:
                        seq+=r_query_sequence[rpos:rpos+inc_read]
                        qual+=r_query_qualities[rpos:rpos+inc_read]
                        rpos+=inc_read

                
                
                read = pysam.AlignedSegment()
                read.query_sequence = seq
                read.query_qualities = qual
                read.reference_start = start_pos
                read.cigartuples=genome_cigatuples
                read.reference_id=dict_chr2idx[iso.t.transcript.Chromosome]
                read.mapping_quality = 60
                read.set_tag(tag="NM", value=len(mismatches)+n_insertions+n_deletions, value_type="i")
                # FIXME add MD tag, see https://samtools.github.io/hts-specs/SAMtags.pdf
                read.flag = 0 if iso.strand=='+' else 16
                # read name encodes: true_tid, true_strand, true_isoform, read_tag, true_chr, true_cigar, n_seqerr, n_converted
                read.query_name =  '_'.join([str(x) for x in [r.query_name, # true_tid, true_strand, true_isoform, tag
                                                    iso.t.transcript.Chromosome, # true_isoform
                                                    read.reference_start,
                                                    read.cigarstring, # cigarstring for calculating aligned blocks
                                                    len(mismatches)+n_insertions+n_deletions, # number of seq errors
                                                    0 # number of conversions
                                                    ]])
                
                # skip reads with too-long read names. Samtools maximum length of read name is 254.
                if len(read.query_name)>254:
                    stats['skipped_long_read_name']+=1
                    print('skipping read due to too long read_name %s' % read_name)
                    continue

                bam_out.write(read)
                stats['output_read']+=1
                stats['output_read_'+iso.strand]+=1        
   
    try:
        pysam.sort("-o", out_bam, out_bam+'.tmp.bam') # @UndefinedVariable
        os.remove(out_bam+'.tmp.bam')
        pysam.index(out_bam) # @UndefinedVariable
    except Exception as e:
        print("error sorting+indexing bam: %s" % e)
    print("All done.")
    print(stats)
    
    
def modify_bases(ref, alt, seq, is_reverse, conversion_rate):
    """ introduces TC / AG conversions
        returns the new sequence and the converted positions (always! 5'-3' as shown in genome browser)
        positions are 0-based """
    if is_reverse:
        ref=COMP_TABLE[ref]
        alt=COMP_TABLE[alt]
    convseq=""
    n_modified=0
    n_convertible=0
    for i,c in enumerate(seq):
        if c == ref:
            n_convertible+=1
            if random.uniform(0, 1) < conversion_rate:
                c=alt.lower()
                n_modified+=1
        convseq+=c
    return convseq, n_convertible, n_modified

def add_read_modifications(in_bam, out_bam, ref, alt, conversion_rate, tag_tc="xc"):
    """ modify single nucleotides in reads """
    sam = pysam.AlignmentFile(in_bam, "rb")
    with pysam.AlignmentFile(out_bam+'.tmp.bam', "wb", header=sam.header) as bam_out:
        for r in sam.fetch(until_eof=True):
            quals=r.query_qualities # save qualities as they will be deleted if query sequence is set
            r.query_sequence, n_convertible, n_modified=modify_bases(ref, alt, r.query_sequence, r.is_reverse, conversion_rate)
            r.set_tag(tag=tag_tc, value=n_modified, value_type="i")
            if n_modified>=0:
                # set blue read color intensity relative to number of tc conversions!
                intensity = min(255, 200.0*(1-n_modified/n_convertible))
                if r.is_reverse:
                    r.set_tag(tag='YC', value='%i,%i,255' % (intensity,intensity) )
                else:
                    r.set_tag(tag='YC', value='255,%i,%i' % (intensity,intensity) )
            true_tid, true_strand, true_isoform, read_tag, true_chr, true_start, true_cigar, n_seqerr, n_converted = r.query_name.split('_')
            n_converted=n_modified
            r.query_name='_'.join([str(x) for x in [true_tid, true_strand, true_isoform, read_tag, true_chr, true_start, true_cigar, n_seqerr, n_converted]])
            r.query_qualities=quals
            bam_out.write(r)
    try:
        pysam.sort("-o", out_bam, out_bam+'.tmp.bam') # @UndefinedVariable
        os.remove(out_bam+'.tmp.bam')
        pysam.index(out_bam) # @UndefinedVariable
    except Exception as e:
        print("error sorting+indexing bam: %s" % e)
    print("All done.")
    
    
def bam_to_fastq(in_bam, out_fastq):
    """ Convert a BAM to fastq """
    sam = pysam.AlignmentFile(in_bam, "rb")
    with open(out_fastq, 'w') as out_fq:
        for r in sam.fetch(until_eof=True):
            print(r)
            if r.is_reverse:
                r_query_sequence=reverse_complement(r.query_sequence)
                print(r.query_qualities)
                r_query_qualities=list(reversed(r.query_qualities))
            else:
                r_query_sequence=r.query_sequence
                r_query_qualities=r.query_qualities
            #write fq
            qstr = ''.join(map(lambda x: chr( x+33 ), r_query_qualities))
            print("@%s\n%s\n+\n%s" % ( r.query_name, r_query_sequence, qstr), file=out_fq )    
            
            
# DELME
def test_cigar_overlap():
    r1 = pysam.AlignedSegment()
    r1.cigarstring='100='
    r1.reference_start=45122380
    
    r2 = pysam.AlignedSegment()
    r2.cigarstring='15=1I28=1I12=1I11=1I18=1I11='
    r2.reference_start=45122380
    
    p1=set([x for _,x in r1.get_aligned_pairs() if x is not None])
    p2=set([x for _,x in r2.get_aligned_pairs() if x is not None])
    print(len(p1), len(p2), len(p1 & p2))

if __name__ == '__main__':
    
#     parser = argparse.ArgumentParser()
#     parser.add_argument("-i", "--in", type=str, required=True, dest="transcript_bam_file", help="Input transcript BAM file")
#     parser.add_argument("-c", "--config", type=str, required=True, dest="config_file", help="Input configuration file")
#     parser.add_argument("-of", "--out_fastq", type=str, required=True, dest="out_fastq",help="Output FASTQ file")
#     parser.add_argument("-ob", "--out_bam", type=str, required=True, dest="out_bam",help="Output BAM file")
#     args = parser.parse_args()
#     
#     m = load_model_from_config(args.config_file)
#     transcript2genome_bam(args.transcript_bam_file, m, args.out_fastq, args.out_bam)

#     add_read_modifications('/Users/niko.popitsch/Desktop/data/projects/Ameres/splicing/splice_sim/testruns/simple/simple/sim/tmp/simple.c_01.simulated.bam',
#                            '/Users/niko.popitsch/Desktop/data/projects/Ameres/splicing/splice_sim/testruns/simple/simple/sim/tmp/simple.c_01.simulated+tc.bam',
#                            'C', 'T',
#                            0.03)

    bam_to_fastq('/Users/niko.popitsch/Desktop/data/projects/Ameres/splicing/splice_sim/testruns/simple/simple/sim/tmp/simple.c_01.simulated+tc.bam',
                 '/Users/niko.popitsch/Desktop/data/projects/Ameres/splicing/splice_sim/testruns/simple/simple/sim/tmp/simple.c_01.simulated+tc.bam.fq')