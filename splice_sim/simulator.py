'''
@author: niko.popitsch@imba.oeaw.ac.at
'''
import csv, datetime, time, logging, sys, os, json
import pysam
from splice_sim.utils import *
from collections import Counter, OrderedDict
import random


default_isoform_colors = {
    ("mat","+"): "100,10,10",
    ("mat","-"): "10,10,100",
    ("old","+"): "200,200,200",
    ("old","-"): "200,200,200",
    ("pre","+"):"200,100,100",
    ("pre","-"):"100,100,200"
    }
    
def postfilter_bam( config, bam_in, out_dir):
    """ Filter secondary+supplementary reads and highlight isoforms """  
    # update isoform_colors from config
    isoform_colors=default_isoform_colors
    if 'isoform_colors' in config:
        for k in config['isoform_colors'].keys():
            isoform_colors[(k,'+')]=config['isoform_colors'][k][0]
            isoform_colors[(k,'-')]=config['isoform_colors'][k][1]
            
    samin = pysam.AlignmentFile(bam_in, "rb")
    bam_out = out_dir+'/'+os.path.splitext(bam_in)[0]+'.final.bam'
    mapped_out = pysam.AlignmentFile(bam_out, "wb", template=samin )
    n_reads=0
    f_reads=0
    for read in samin.fetch(until_eof=True): # NOTE important to add until_eof=True otherwise unmapped reads will be skipped.
        if read.is_secondary or read.is_supplementary:
            f_reads+=1
            continue
        # map isoform to color
        true_tid, true_strand, true_isoform, read_tag, true_chr, true_start, true_cigar, n_seqerr, n_converted, is_converted_read = read.query_name.split('_')   
        if (true_isoform, true_strand) in isoform_colors:
            read.set_tag(tag='YC', value=isoform_colors[true_isoform, true_strand], value_type="Z")
        mapped_out.write(read)
        n_reads+=1  
    samin.close()
    mapped_out.close()
    
    try:
        pysam.index(bam_out) # @UndefinedVariable
    except Exception as e:
        print("error indexing bam: %s" % e)
    return True, n_reads, f_reads

def transcript2genome_bam(m, transcript_bam_file, out_bam, overwrite=False):
    """ Reads a transcriptome-mapped BAM and converts to a genome BAM. 
    Transcript ids, strand, isoform and tags are parsed from the input read names and matched with the model """
    stats=Counter()
    if not overwrite and os.path.isfile(out_bam):
        print("Pre-existing file %s, will not recreate" % out_bam)
        return stats
    genome = pysam.FastaFile(m.config["genome_fa"])    
    chromosomes = genome.references
    dict_chr2idx = {k: v for v, k in enumerate(chromosomes)}
    dict_idx2chr = {v: k for v, k in enumerate(chromosomes)}
    dict_chr2len = {c: genome.get_reference_length(c) for c in genome.references}
    header = { 'HD': {'VN': '1.4'}, 
               'SQ': [{'LN': l, 'SN': c} for c,l in dict_chr2len.items()] }
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
                iso = m.transcripts[tid].isoforms[iso_id]
                is_converted_read = 1 if iso.is_labeled else 0 # does the read stem form a labeled or unlabeled isoform? 
                ablocks=iter(iso.aln_blocks)
                genome_cigatuples=[]
                seq=""
                qual=[]
                rpos=0
                gpos=iso.transcript.start-1
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
                    while diff>0:
                        if op==8 and l>1:
                            #FIXME seq err stretches that are spliced are partially ignored
                            print('ignoring seq diff stretch for ' + r.query_name+', '+r.cigarstring) # TODO we need to insert a splice cigartuple
                        matchlen=l-diff
                        if matchlen>0:
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
                            print("ERR: no more ablocks in iso " + str(iso) + ", read " + r.query_name+', '+r.cigarstring) # should not happen :-/
                            print(iso.aln_blocks)
                            break
                        else:
                            n_len = nblock[0]-ablock[1]-1
                            if n_len>0:
                                genome_cigatuples+=[(3, n_len)]
                            gpos+=n_len
                            # update l 
                            ablock=nblock
                            l-=matchlen
                            diff=gpos+inc_pos-ablock[1]
                    if l>0: 
                        genome_cigatuples+=[(op, l)]
                    if op == 8 : # X - store mismatches
                        for i in range(l):
                            mismatches+=[gpos-start_pos+1] 
                    gpos+=inc_pos
                    if inc_read>0:
                        seq+=r_query_sequence[rpos:rpos+inc_read]
                        qual+=r_query_qualities[rpos:rpos+inc_read]
                        rpos+=inc_read
                # create read
                read = pysam.AlignedSegment()
                read.query_sequence = seq
                read.query_qualities = qual
                read.reference_start = start_pos
                read.cigartuples=genome_cigatuples
                read.reference_id=dict_chr2idx[iso.transcript.chromosome]
                read.mapping_quality = 60
                read.set_tag(tag="NM", value=len(mismatches)+n_insertions+n_deletions, value_type="i")
                if is_converted_read==0:
                    r.set_tag(tag='YC', value='200,200,200') # old rna
                # FIXME add MD tag, see https://samtools.github.io/hts-specs/SAMtags.pdf
                read.flag = 0 if iso.strand=='+' else 16
                # read name encodes: true_tid, true_strand, true_isoform, read_tag, true_chr, true_cigar, n_seqerr, n_converted
                read.query_name =  '_'.join([str(x) for x in [r.query_name, # true_tid, true_strand, true_isoform, tag
                                                    iso.transcript.chromosome, # true_isoform
                                                    read.reference_start,
                                                    read.cigarstring, # cigarstring for calculating aligned blocks
                                                    len(mismatches)+n_insertions+n_deletions, # number of seq errors
                                                    0, # number of conversions
                                                    is_converted_read # read is chosen from the binomial distribution
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
    except Exception as e:
        print("error sorting: %s" % e)
    try:
        if os.path.exists(out_bam+'.bai'):
            os.remove(out_bam+'.bai')
        print("Indexing BAM ", out_bam)
        pysam.index(out_bam) # @UndefinedVariable
    except Exception as e:
        print("error indexing bam: %s" % e)    
    return(stats)
    
    
def modify_bases(chrom, start, ref, alt, seq, is_reverse, conversion_rate, conversion_mask_prob, tested_positions, masked_positions ):
    """ introduces nucleotide conversions.
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
            # is this a masked position?
            if conversion_mask_prob is not None:                
                pos=chrom+':'+str(int(start)+i)
                if pos in tested_positions and pos in masked_positions:
                    convseq+=c
                    continue # skip conversion
                tested_positions.add(pos)
                if random.uniform(0, 1) < conversion_mask_prob:
                    masked_positions.add(pos)
                    convseq+=c
                    continue # skip conversion
            if random.uniform(0, 1) < conversion_rate:
                c=alt.lower()
                n_modified+=1
        convseq+=c
    return convseq, n_convertible, n_modified, tested_positions, masked_positions

def add_read_modifications(in_bam, out_bam, ref, alt, conversion_rate, conversion_mask_prob, threads, tag_tc="xc", tag_isconvread="xm", overwrite=False):
    """ modify single nucleotides in reads. 
    """
    if not overwrite and os.path.isfile(out_bam):
        print("Pre-existing file %s, will not recreate" % out_bam)
        return
    sam = pysam.AlignmentFile(in_bam, "rb")
    masked_positions = set()
    tested_positions = set()
    with pysam.AlignmentFile(out_bam+'.tmp.bam', "wb", header=sam.header) as bam_out:
        for r in sam.fetch(until_eof=True):
            quals=r.query_qualities # save qualities as they will be deleted if query sequence is set
            true_tid, true_strand, true_isoform, read_tag, true_chr, true_start, true_cigar, n_seqerr, n_converted, is_converted_read  = r.query_name.split('_')
            is_converted_read=int(is_converted_read)
            if is_converted_read==1: # read chosen from binomial; modify its bases
                r.query_sequence, n_convertible, n_modified, tested_positions, masked_positions=modify_bases(true_chr, 
                                                                                                             true_start, ref, alt, r.query_sequence, 
                                                                                                             r.is_reverse, conversion_rate, 
                                                                                                             conversion_mask_prob, tested_positions, masked_positions)
            else:
                n_convertible=r.query_sequence.count(ref)
                n_modified=0
            r.set_tag(tag=tag_tc, value=n_modified, value_type="i")
            nm=r.get_tag('NM') if r.has_tag('NM') else 0
            r.set_tag(tag='NM', value=nm+1, value_type="i") # update NM tag
            r.set_tag(tag=tag_tc, value=n_modified, value_type="i")
            r.set_tag(tag=tag_isconvread, value=is_converted_read, value_type="i")
            if is_converted_read==0:
                r.set_tag(tag='YC', value='200,200,200') # old rna
            else:
                # set blue read color intensity relative to number of tc conversions!
                intensity = min(255, 200.0*(1-n_modified/n_convertible)) if n_convertible>0 else 255
                if r.is_reverse:
                    r.set_tag(tag='YC', value='%i,%i,255' % (intensity,intensity) )
                else:
                    r.set_tag(tag='YC', value='255,%i,%i' % (intensity,intensity) )

            n_converted=n_modified
            r.query_name='_'.join([str(x) for x in [true_tid, true_strand, true_isoform, read_tag, true_chr, true_start, true_cigar, n_seqerr, n_converted, is_converted_read]])
            r.query_qualities=quals
            bam_out.write(r)
        if conversion_mask_prob is not None:
            print("masked_positions %i/%i : " % ( str(len(masked_positions)),  str(len(tested_positions)) ) )
    try:
        pysam.sort("-@", str(threads), "-o", out_bam, out_bam+'.tmp.bam') # @UndefinedVariable
        os.remove(out_bam+'.tmp.bam')
        pysam.index(out_bam) # @UndefinedVariable
    except Exception as e:
        print("error sorting+indexing bam: %s" % e)
    print("All done.")
    
    
def bam_to_fastq(in_bam, out_fastq, overwrite=False):
    """ Convert a BAM to fastq """
    if not overwrite and os.path.isfile(out_fastq):
        print("Pre-existing file %s, will not recreate" % out_fastq)
        return
    sam = pysam.AlignmentFile(in_bam, "rb")
    with open(out_fastq, 'w') as out_fq:
        for r in sam.fetch(until_eof=True):
            if r.is_reverse:
                r_query_sequence=reverse_complement(r.query_sequence)
                r_query_qualities=list(reversed(r.query_qualities))
            else:
                r_query_sequence=r.query_sequence
                r_query_qualities=r.query_qualities
            #write fq
            qstr = ''.join(map(lambda x: chr( x+33 ), r_query_qualities))
            print("@%s\n%s\n+\n%s" % ( r.query_name, r_query_sequence, qstr), file=out_fq )
            


def create_genome_bam(m, art_sam_file, threads, out_dir):
    """ create genome bam files per condition and export to fastq """
    config=m.config

    # file names
    if not os.path.exists(out_dir):
        os.mkdirs(out_dir)
        
    # convert to genome BAM (truth file w/o read modifications)            
    f_bam_ori = out_dir+'/'+config["dataset_name"]+".cr0.truth.bam"
    f_fq_ori = out_dir+'/'+config["dataset_name"]+".cr0.truth.fq"
    stats=transcript2genome_bam(m, art_sam_file, f_bam_ori)
    print(stats)
    # convert to FASTQ
    bam_to_fastq(f_bam_ori, f_fq_ori)
        
    # add read modifications
    created_bams=[f_bam_ori]
    created_fqs=[f_fq_ori]
    for cr in m.condition.conversion_rates:
        f_bam_mod = out_dir+'/'+config["dataset_name"]+".cr"+str(cr)+".truth.bam"
        f_fq_mod = out_dir+'/'+config["dataset_name"]+".cr"+str(cr)+".truth.fq"
        add_read_modifications(f_bam_ori, 
                               f_bam_mod, 
                               m.condition.ref, 
                               m.condition.alt,
                               cr, 
                               m.conversion_mask_prob,
                               threads)            
            
        # convert to FASTQ
        bam_to_fastq(f_bam_mod,  f_fq_mod)
        created_bams+=[f_bam_mod]
        created_fqs+=[f_fq_mod]
    
    return created_bams, created_fqs