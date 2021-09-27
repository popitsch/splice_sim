'''
@author: niko.popitsch@imba.oeaw.ac.at
'''


from argparse import ArgumentParser, RawDescriptionHelpFormatter, \
    ArgumentTypeError
from collections import *
import csv, datetime, time, logging, sys, os, json
import pysam
import pyranges as pr
import pandas as pd
import numpy as np
from utils import bgzip, files_exist, pad_n, pipelineStep, reverse_complement, COMP_TABLE, removeFile,sambamba2bam
import random
import math
import gzip
import logging
from subprocess import check_output
from model import Model, Condition, Isoform, Transcript
from config_creator import calculate_transcript_data
from transcript2genome_bam import transcript2genome_bam
from pathlib import Path

def check_config(config):
    for section in ['dataset_name','mappers','genome_fa','gene_gff','conditions','isoform_mode','transcript_data']:
        assert section in config, "Missing configuration section %s" % section
def get_config(config, keys, default_value): 
    """ gets value from dict of dicts or default_value if path does not exist """
    d = config
    for k in keys:
        if k not in d:
            return default_value
        d = d[k]
    return d        
        
#============================================================================
# Mapping methods
#============================================================================
def get_version(exe):
    cmd = [exe, "--version"]
    return check_output(" ".join(cmd), shell=True)


def runSTAR(bam, ref, reads1=[], reads2=None,  gtf=None, force=True, run_flagstat=False, threads=1, chimeric=0, writeUnmapped=True, doSort=True, additionalParameters=[], STAR_EXE='STAR', SAMBAMBA_EXE='sambamba'):
    """ Run STAR, supports SE+PE and can merge multiple files (needs sambamba access). Can create samtools flagstats """
    success = True
    if files_exist(bam) and not force:
        logging.warn("BAM file " + bam + " already exists! use -f to recreate.")
        return True
    # check input params
    if len(reads1)==0 or (reads2 is not None and len(reads1)!=len(reads2)):
        logging.error("ERROR")
        return False
    iscompressed=reads1[0].endswith('.gz')
    if (type(reads1) is not list):
        reads1=[reads1]
        reads2=[reads2]
    # PAIRED or SINGLE END
    bam_chimeric = os.path.splitext(bam)[0] + ".chimeric.bam"
    tomerge=[]
    tomerge_chimeric=[]
    for idx, val in enumerate(reads1):
        pre = os.path.splitext(bam)[0]+"."+str(idx)+"."
        starout = pre+"Aligned.out.sam"
        sorted = pre+".bam"
        starout_chimeric = pre+"Chimeric.out.sam"
        sorted_chimeric = pre+".chimeric.bam"
        cmd=[STAR_EXE]
        if ( iscompressed ):
            cmd+=["--readFilesCommand", "zcat"]
        cmd+=["--outFileNamePrefix", pre]
        cmd+=["--runThreadN", str(threads)]
        cmd+=["--genomeDir", ref]
        cmd+=["--outSAMattributes", "NH HI AS nM MD"] # @see http://samtools.github.io/hts-specs/SAMtags.pdf
        if reads2 is not None:
            cmd+=["--readFilesIn", reads1[idx], reads2[idx]]
        else:
            cmd+=["--readFilesIn", reads1[idx]]
        if chimeric > 0:
            cmd+=["--chimSegmentMin", str(chimeric)]
            cmd+=["--chimOutType", "SeparateSAMold"]
        if writeUnmapped:
            #cmd+=["--outReadsUnmapped", "Fastx"]
            cmd+=["--outSAMunmapped", "Within"]
        if gtf:
            cmd+=["--sjdbGTFfile", gtf]
        cmd += additionalParameters
        logging.info(cmd)
        success = success and pipelineStep([reads1,reads2], starout, cmd, shell=True) 
        logging.info(starout)
        success = success and sambamba2bam(starout, sorted,     
                                          sort=doSort, 
                                          index=True, 
                                          delinFile=True,
                                          ncpu=threads, EXE=SAMBAMBA_EXE) 
        tomerge+=[sorted]
        # handle chimeric reads if configured
        if files_exist(starout_chimeric):
            success = success and sambamba2bam(starout_chimeric, sorted_chimeric, 
                                          sort=doSort, 
                                          index=True, 
                                          delinFile=True,
                                          ncpu=threads, EXE=SAMBAMBA_EXE) 
            tomerge_chimeric+=[sorted_chimeric]
    #print(tomerge)
    # merge files
    if len(tomerge)==1:
        os.rename(tomerge[0], bam)
        os.rename(tomerge[0]+".bai", bam+".bai")
    else:
        logging.info("Merging from multiple input files")
        success = success and sambambamerge(bam, tomerge, ncpu=threads, EXE=SAMBAMBA_EXE)
    if len(tomerge_chimeric)==0:
        logging.info("NOTE: no chimeric read output configured!")
    else:
        if len(tomerge_chimeric)==1:
            os.rename(tomerge_chimeric[0], bam_chimeric)
            os.rename(tomerge_chimeric[0]+".bai", bam_chimeric+".bai")
        else:
            logging.info("Merging from multiple input files")
            success = success and sambambamerge(bam_chimeric, tomerge_chimeric, ncpu=threads, EXE=SAMBAMBA_EXE)
        flagstat( bam_chimeric )
    if run_flagstat:            
        flagstat( bam )        
    return success


def runHISAT_3N(bam, fq, idx, base_ref='T', base_alt='C', known_splicesites=None, force=True, run_flagstat=False, threads=1, doSort=True, additionalParameters=[], HISAT_3N_EXE='hisat-3n', SAMBAMBA_EXE='sambamba'):
    """ Run HISAT-3N, see http://daehwankimlab.github.io/hisat2/hisat-3n/ """
    success = True
    if files_exist(bam) and not force:
        logging.warn("BAM file " + bam + " already exists! use -f to recreate.")
        return True
    
    todel=[]
    # temporary unzip as hisat2 cannot read gzipped inputs
    if fq.endswith(".gz"):
        fq1=os.path.splitext(fq)[0]
        cmd=[ "gunzip", "-c", fq]
        logging.info(cmd)
        success = success and pipelineStep([fq], [fq1], cmd, shell=True, stdout=fq1)
        todel+=[fq1]
        fq=fq1
    # run hisat-3n
    sam = bam+".sam"
    cmd=[ HISAT_3N_EXE,
         "--base-change", "%s,%s"%(base_ref, base_alt),
         "--reference", ref,
         "--index", idx,
         "-U", fq,
         "--threads", str(threads),
         "-S", sam
        ]
    if known_splicesites is not None:
        cmd+=[ "--known-splicesite-infile", known_splicesites ]
    logging.info(cmd)
    success = success and pipelineStep([fq], sam, cmd, shell=True)
    
    # sort + index
    success = success and sambamba2bam(sam, bam, 
                              sort=doSort, 
                              index=True, 
                              delinFile=True,
                              ncpu=threads, 
                              EXE=SAMBAMBA_EXE) 
    
    removeFile(todel)
    if run_flagstat:            
        flagstat( bam )        
    return success


def runMERANGS(bam, fq, idx, known_splicesites=None, force=True, run_flagstat=False, threads=1, additionalParameters=[], MERANG_EXE='meRanGs', STAR_EXE='STAR', FASTQTOBAM_EXE='fastqtobam'):
    """ Run meRanGs, see http://daehwankimlab.github.io/hisat2/hisat-3n/ """
    success = True
    if files_exist(bam) and not force:
        logging.warn("BAM file " + bam + " already exists! use -f to recreate.")
        return True
    # get base name for FQ
    fq_name=Path(fq).stem
    if fq_name.endswith('fq'):
        fq_name=Path(fq_name).stem
    udir=str(Path(bam).parent.absolute())+"/meRanGs_unaligned_reads" # dir for storing unaligned reads 
    
    # get base name for sam
    sam_name=Path(bam).stem
    
    todel=[]
    cmd=[ MERANG_EXE, "align",
         "-f", fq,
         "-t", str(threads),
         "-S", sam_name+'.sam',
         "-un", "-ud", udir,  # Report unaligned reads
         #"-MM", "-star_outFilterMultimapNmax", "20"  # Save multimappers?
         "-id", idx,
         "-mbp", # m‐Bias plot (for QC only)
         "-starcmd", STAR_EXE, # STAR cmd to use
         "-star_outSAMattributes", "NH HI AS nM MD"  # @see http://samtools.github.io/hts-specs/SAMtags.pdf
         ] + additionalParameters
    if known_splicesites is not None:
        cmd+=[ "-star_sjdbGTFfile", known_splicesites ]
    logging.info(cmd)
    success = success and pipelineStep([fq], sam, cmd, shell=True)
        
    # converted unaligned FASTQ to BAM
    ureadsfq=udir+'/'+fq_name+'_unmapped.fq.gz'
    cmd=[ FASTQTOBAM_EXE, "gz=1", ureadsfq]
    logging.info(cmd)
    success = success and pipelineStep([ureadsfq], sam, cmd, shell=True, stdout=udir+'/'+fq_name+'_unmapped.bam')
    
    # merge aligend and unaligned bams
    cmd=[ 'samtools', 'merge', bam, sam_name+'_sorted.bam', udir+'/'+fq_name+'_unmapped.bam']
    logging.info(cmd)
    success = success and pipelineStep([sam_name+'_sorted.bam', udir+'/'+fq_name+'_unmapped.bam'], bam, cmd, shell=True)
    if run_flagstat:            
        flagstat( bam )        
    return success

@DeprecationWarning
def runHISAT2_TLA(bam, ref, fq, idx1, idx2, base_ref='T', base_alt='C', known_splicesites=None, force=True, run_flagstat=False, threads=1, doSort=True, additionalParameters=[], HISAT2_EXE='hisat2', SAMBAMBA_EXE='sambamba'):
    """ Run HISAT2_TLA. Depreated """
    success = True
    if files_exist(bam) and not force:
        logging.warn("BAM file " + bam + " already exists! use -f to recreate.")
        return True
    
    todel=[]
    # temporary unzip as hisat2 cannot read gzipped inputs
    if fq.endswith(".gz"):
        fq1=os.path.splitext(fq)[0]
        cmd=[ "gunzip", "-c", fq]
        success = success and pipelineStep([fq], [fq1], cmd, shell=True, stdout=fq1)
        todel+=[fq1]
        fq=fq1
    
    # run hisat2_tla
    sam = bam+".sam"
    cmd=[ HISAT2_EXE,
         "--TLA",
         "--base-change", "%s%s"%(base_ref, base_alt),
         "--reference", ref,
         "--index1", idx1,
         "--index2", idx2,
         "-U", fq,
         "--threads", str(threads),
         "-S", sam
        ]
    if known_splicesites is not None:
        cmd+=[ "--known-splicesite-infile", known_splicesites ]
    logging.info(cmd)
    success = success and pipelineStep([fq], sam, cmd, shell=True)
    
    # sort + index
    success = success and sambamba2bam(sam, bam, 
                              sort=doSort, 
                              index=True, 
                              delinFile=True,
                              ncpu=threads, 
                              EXE=SAMBAMBA_EXE) 
    
    removeFile(todel)
    if run_flagstat:            
        flagstat( bam )        
    return success

# 
# create tdf file
#
def create_tdf(bam, tdf, chrom_sizes, igvtools_cmd='igvtools', logfile="igvtools.log"):
    """ Create TDF files """
    cmd=[igvtools_cmd, "count", "-w", "1", "-z", "7", "--strands", "read", bam, tdf, chrom_sizes]
    return pipelineStep(bam, tdf, cmd, shell=True, stdout=logfile)

def postfilter_bam( bam_in, bam_out):
    """ Filter secondary+supplementary reads """  
    samin = pysam.AlignmentFile(bam_in, "rb")
    mapped_out = pysam.AlignmentFile(bam_out, "wb", template=samin )
    n_reads=0
    f_reads=0
    for read in samin.fetch(until_eof=True): # NOTE important to add until_eof=True otherwise unmapped reads will be skipped.
        if read.is_secondary or read.is_supplementary:
            f_reads+=1
            continue
        mapped_out.write(read)
        n_reads+=1  
    samin.close()
    mapped_out.close()
    
    try:
        pysam.index(bam_out) # @UndefinedVariable
    except Exception as e:
        print("error indexing bam: %s" % e)
    return True, n_reads, f_reads

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
                    while diff>0:
                        if op==8 and l>1:
                            #FIXME seq err stretc1Ghes that are spliced are partially ignored
                            print('ignoring seq diff stretch for ' + r.query_name+', '+r.cigarstring) # TODO we need to insert a splice cigartuple
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
                            print("ERR: no more ablocks in iso " + str(iso) + ", read " + r.query_name+', '+r.cigarstring) # should not happen :-/
                            print(iso.aln_blocks)
                            break
                        else:
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
                # create read
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
    return(stats)
    
    
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
                intensity = min(255, 200.0*(1-n_modified/n_convertible)) if n_convertible>0 else 255
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
            if r.is_reverse:
                r_query_sequence=reverse_complement(r.query_sequence)
                r_query_qualities=list(reversed(r.query_qualities))
            else:
                r_query_sequence=r.query_sequence
                r_query_qualities=r.query_qualities
            #write fq
            qstr = ''.join(map(lambda x: chr( x+33 ), r_query_qualities))
            print("@%s\n%s\n+\n%s" % ( r.query_name, r_query_sequence, qstr), file=out_fq )
            

def simulate_dataset(config, config_dir, outdir, overwrite=False):
    """ Simulates a dataset. """
    startTime = time.time()
    logging.info("Simulating dataset %s" % config['dataset_name'])
    
    # output dirs
    if not os.path.exists(outdir):
        print("Creating dir " + outdir)
        os.makedirs(outdir)
    tmpdir = outdir + "/tmp/"
    if not os.path.exists(tmpdir):
        os.makedirs(tmpdir)
    
    # write ORI files?
    write_uncoverted=config["write_uncoverted"] if "write_uncoverted" in config else True
    
    # executables
    art_cmd =      config["executables"]["art_cmd"] if "executables" in config and "art_cmd" in config["executables"] else 'art_illumina'
    sambamba_cmd = config["executables"]["sambamba_cmd"] if "executables" in config and "sambamba_cmd" in config["executables"] else 'sambamba'
    igvtools_cmd = config["executables"]["igvtools_cmd"] if "executables" in config and "igvtools_cmd" in config["executables"] else 'igvtools'
    threads=config["threads"] if "threads" in config else 1
        
    if overwrite:
        logging.info("NOTE: existing result files will be overwritten!")
    
    # read transcript data from external file if not in config
    if 'transcripts' in config:
        logging.info("Reading transcript configuration from config file.")
    else:
        assert "transcript_data" in config, "Transcript data needs to be configured either in config file ('transcripts' section) or in an external file referenced via 'transcript_data'"
        tfile = config['transcript_data']
        if not os.path.isabs(tfile):
            tfile = config_dir+"/"+tfile
        if not files_exist(tfile):
            logging.info("Creating external transcript config file %s" % tfile)
            calculate_transcript_data(config, config_dir, outdir)
            
        logging.info("Reading transcript configuration from external config file %s" % tfile)
        tdata = json.load(open(tfile), object_pairs_hook=OrderedDict)
        config["transcripts"]=tdata
    
    
    # instantiate model
    m = Model( config )
    logging.info("Configured conditions: %s" % (", ".join(str(x) for x in m.conditions)) )
    
    # write considered transcripts to GFF
    m.write_gff(outdir)
    
    # add number of transcripts to stats
    stats=Counter()
    stats["transcripts"]= len(m.transcripts)
    
    # now write one fasta file per condition with transcript/isoform entries that are proportional to the configured abundance
    logging.info("Calculating isoform data")
    for cond in m.conditions:
        fout = tmpdir + config['dataset_name'] + "." + cond.id + ".fa"
        if overwrite or (not files_exist(fout) and not files_exist(fout+".gz")):
            buf=[]
            with open(fout, 'w') as out:
                for t in m.transcripts.values():
                    sequences = t.get_sequences()
                    #print("Writing to %s" % (fout) )            
                    for iso in t.isoforms.keys():
                        for i in range(sequences[cond][iso][1]):
                            buf+=[">%s_%s_%s_%i\n" %(t.tid, t.transcript.Strand, iso, i)]
                            buf+=[pad_n(sequences[cond][iso][0], config['readlen'])+"\n"]
                            if len(buf)>10000:
                                out.writelines(buf)
                                buf=[]
                if len(buf)>0:
                    out.writelines(buf)
        else:
            logging.info("Will not re-create existing file %s" % (fout))
    
    # run art simulator 
    logging.info("Running read simulator")
    artlog=tmpdir + config['dataset_name'] + ".art.log" 
    stats=Counter()
    for cond in m.conditions:
        bamdir_truth = outdir + "bam_ori/TRUTH/"
        bamdir_truth_conv = outdir + "bam_conv/TRUTH/"
        if not os.path.exists(bamdir_truth):
            os.makedirs(bamdir_truth)
        if not os.path.exists(bamdir_truth_conv):
            os.makedirs(bamdir_truth_conv)
                    
        f = tmpdir + config['dataset_name'] + "." + cond.id + ".fa"
        # art output
        art_out_prefix = tmpdir + config['dataset_name'] + "." + cond.id + ".art_output"
        f_fq_art = art_out_prefix + ".fq"
        f_sam_art =  art_out_prefix + ".sam"
        # output BAM/fastq files
        f_bam = bamdir_truth + config['dataset_name'] + "." + cond.id + ".simulated.bam"
        f_bam_conv = bamdir_truth_conv + config['dataset_name'] + "." + cond.id + ".simulated+conversions.bam"
        f_fq = tmpdir + config['dataset_name'] + "." + cond.id + ".simulated.fq"
        f_fq_conv = tmpdir + config['dataset_name'] + "." + cond.id + ".simulated+conversions.fq"
        additional_art_params = get_config(config, ['simulator','additional_art_params'],"") # e.g., --insRate (default: 0.00009)
        if overwrite or not files_exist([f_bam]):
            # NOTE: use 2x coverage as ~50% of reads will be simulated for wrong strand and will be dropped in postprocessing
            cmd=[art_cmd, "-ss", "HS25", "-i", f, "-l", str(config["readlen"]), "-f", str(cond.coverage * 2), "-na", "--samout", "-o", art_out_prefix, additional_art_params ]
            if "random_seed" in config:
                cmd+=["-rs", str(config["random_seed"])]
            success = pipelineStep(f, [f_fq_art, f_sam_art], cmd, shell=True, stdout=artlog, append=True)    
            if not success:
                logging.error('error simulating %s' % (f_fq_art))
                
            # convert to genome BAM (truth file w/o read modifications)            
            stats[cond.id]=transcript2genome_bam(f_sam_art, m, f_bam)
            
            # convert to FASTQ
            bam_to_fastq(f_bam, f_fq)
            bgzip(f_fq, threads=threads)

            # add read modifications
            add_read_modifications(f_bam, f_bam_conv, cond.ref, cond.alt, cond.conversion_rate)            
                    
            # convert to FASTQ
            bam_to_fastq(f_bam_conv,  f_fq_conv)
            bgzip(f_fq_conv, threads=threads)
                
            if success:
                if not get_config(config, ['simulator','keep_sam'], False):
                    removeFile([f_fq_art, f_sam_art])
    
        else:
            logging.warn("Will not re-create existing file %s" % (f_bam))
    
    print(stats)
            
    # write ROIs 
    logging.info("Writing ROI bed")
    f_roi = outdir + "roi.bed"
    with open(f_roi, 'w') as out:
        for t in m.transcripts.values():
            print(t.to_bed(), file=out)
    
    
    
    # map reads if mappers configured
    bams={}
    if 'mappers' in config:
        logging.info("Mapping reads")
        # print version info
        for mapper in config['mappers'].keys():
            EXE=None
            if mapper == "STAR":
                EXE=config['mappers'][mapper]['star_cmd'] if 'star_cmd' in config['mappers'][mapper] else 'STAR'
            elif mapper == "HISAT2_TLA":
                EXE=config['mappers'][mapper]['hisat2_cmd'] if 'hisat2_cmd' in config['mappers'][mapper] else 'hisat2'
            elif mapper == "MERANGS":
                EXE=config['mappers'][mapper]['merangs_cmd'] if 'merangs_cmd' in config['mappers'][mapper] else 'meRanGs'
            elif mapper == "HISAT-3N":
                EXE=config['mappers'][mapper]['hisat-3n_cmd'] if 'hisat-3n_cmd' in config['mappers'][mapper] else 'hisat-3n'
            if EXE:
                logging.info("%s version: %s" % ( mapper, get_version(EXE)) )

        for cond in m.conditions:
            for mapper in config['mappers'].keys():
                bamdir_all = outdir + "bam_ori/" + mapper + "/"
                bamdir_conv  = outdir + "bam_conv/" + mapper + "/"
                if not os.path.exists(bamdir_all):
                    os.makedirs(bamdir_all)
                if not os.path.exists(bamdir_conv):
                    os.makedirs(bamdir_conv)
                    
                f_fq =      tmpdir + config['dataset_name'] + "." + cond.id + ".simulated.fq.gz"
                f_fq_conv = tmpdir + config['dataset_name'] + "." + cond.id + ".simulated+conversions.fq.gz"
                
                f_bam =      tmpdir + config['dataset_name'] + "." + cond.id + "." +    mapper +".bam"
                f_bam_conv = tmpdir + config['dataset_name'] + "." + cond.id + ".conv."+mapper+".bam"
                
                final_bam      = bamdir_all +   config['dataset_name'] + "." + cond.id +      "."+mapper+".nodup.bam"
                final_bam_conv = bamdir_conv  + config['dataset_name'] + "." + cond.id + ".conv."+mapper+".nodup.bam"
                
                logging.info("Mapping %s reads with %s" % (cond.id, mapper))
                if mapper == "STAR":
                    # run STAR mapper 
                    STAR_EXE=config['mappers'][mapper]['star_cmd'] if 'star_cmd' in config['mappers'][mapper] else 'STAR'
                    star_splice_gtf=config['mappers'][mapper]['star_splice_gtf'] if 'star_splice_gtf' in config['mappers'][mapper] else None
                    star_genome_idx=config['mappers'][mapper]['star_genome_idx']
                    if write_uncoverted:
                        if overwrite or not files_exist(f_bam):
                            runSTAR(f_bam, 
                                    star_genome_idx, 
                                    reads1=[f_fq], 
                                    gtf=star_splice_gtf,
                                    threads=threads, 
                                    STAR_EXE=STAR_EXE,
                                    SAMBAMBA_EXE=sambamba_cmd,
                                    force=overwrite )
                            postfilter_bam(f_bam, final_bam)
                        else:
                            logging.warn("Will not re-create existing file %s" % (f_bam))
                    if overwrite or not files_exist(f_bam_conv):
                        runSTAR(f_bam_conv, 
                            star_genome_idx, 
                            reads1=[f_fq_conv], 
                            gtf=star_splice_gtf,
                            threads=threads, 
                            STAR_EXE=STAR_EXE,
                            SAMBAMBA_EXE=sambamba_cmd,
                            force=overwrite )  
                        postfilter_bam(f_bam_conv, final_bam_conv) 
                    else:
                        logging.warn("Will not re-create existing file %s" % (f_bam_conv))       
                elif mapper == "HISAT-3N":
                    HISAT_3N_EXE=config['mappers'][mapper]['hisat-3n_cmd'] if 'hisat-3n_cmd' in config['mappers'][mapper] else 'hisat-3n'
                    hisat_3n_idx=config['mappers'][mapper]['hisat-3n_idx']
                    hisat_3n_kss=config['mappers'][mapper]['hisat-3n_kss'] if 'hisat-3n_kss' in config['mappers'][mapper] else None
                    if write_uncoverted:
                        if overwrite or not files_exist(f_bam):
                            runHISAT_3N(f_bam, 
                                    fq,
                                    hisat_3n_idx,
                                    base_ref=cond.ref,
                                    base_alt=cond.alt,
                                    known_splicesites=hisat_3n_kss,
                                    threads=threads, 
                                    HISAT_3N_EXE=HISAT_3N_EXE,
                                    SAMBAMBA_EXE=sambamba_cmd,
                                    force=overwrite )
                            postfilter_bam(f_bam, final_bam)
                        else:
                            logging.warn("Will not re-create existing file %s" % (f_bam))
                    if overwrite or not files_exist(f_bam_conv):
                        runHISAT_3N(f_bam_conv, 
                            f_fq_conv,
                            hisat_3n_idx,
                            base_ref=cond.ref,
                            base_alt=cond.alt,
                            known_splicesites=hisat_3n_kss,
                            threads=threads, 
                            HISAT_3N_EXE=HISAT_3N_EXE,
                            SAMBAMBA_EXE=sambamba_cmd,
                            force=overwrite )
                        postfilter_bam(f_bam_conv, final_bam_conv) 
                    else:
                        logging.warn("Will not re-create existing file %s" % (f_bam_conv))                
                elif mapper == "HISAT2_TLA": #@DeprecationWarning
                    HISAT2_EXE=config['mappers'][mapper]['hisat2_cmd'] if 'hisat2_cmd' in config['mappers'][mapper] else 'hisat2'
                    hisat2_idx1=config['mappers'][mapper]['hisat2_idx1']
                    hisat2_idx2=config['mappers'][mapper]['hisat2_idx2']
                    hisat2_kss=config['mappers'][mapper]['hisat2_kss'] if 'hisat2_kss' in config['mappers'][mapper] else None
                    if write_uncoverted:
                        if overwrite or not files_exist(f_bam):
                            runHISAT2_TLA(f_bam, 
                                    config["genome_fa"], 
                                    fq=f_fq,
                                    idx1=hisat2_idx1,
                                    idx2=hisat2_idx2,
                                    base_ref=cond.ref,
                                    base_alt=cond.alt,
                                    known_splicesites=hisat2_kss,
                                    threads=threads, 
                                    HISAT2_EXE=HISAT2_EXE,
                                    SAMBAMBA_EXE=sambamba_cmd,
                                    force=overwrite )
                            postfilter_bam(f_bam, final_bam)
                        else:
                            logging.warn("Will not re-create existing file %s" % (f_bam))
                    if overwrite or not files_exist(f_bam_conv):
                        runHISAT2_TLA(f_bam_conv, 
                            config["genome_fa"], 
                            fq=f_fq_conv,
                            idx1=hisat2_idx1,
                            idx2=hisat2_idx2,
                            known_splicesites=hisat2_kss,
                            threads=threads, 
                            HISAT2_EXE=HISAT2_EXE,
                            SAMBAMBA_EXE=sambamba_cmd,
                            force=overwrite )
                        postfilter_bam(f_bam_conv, final_bam_conv) 
                    else:
                        logging.warn("Will not re-create existing file %s" % (f_bam_conv))
                elif mapper == "MERANGS":
                    # run meRanGs mapper 
                    MERANGS_EXE=config['mappers'][mapper]['merangs_cmd'] if 'merangs_cmd' in config['mappers'][mapper] else 'meRanGs'
                    STAR_EXE=config['mappers'][mapper]['star_cmd'] if 'star_cmd' in config['mappers'][mapper] else 'STAR'
                    merangs_splice_gtf=config['mappers'][mapper]['merangs_splice_gtf'] if 'merang_splice_gtf' in config['mappers'][mapper] else None
                    merangs_genome_idx=config['mappers'][mapper]['merangs_genome_idx']
                    if write_uncoverted:
                        if overwrite or not files_exist(f_bam):
                            runMERANGS(f_bam, 
                                    merangs_genome_idx, 
                                    f_fq, 
                                    gtf=merangs_splice_gtf,
                                    threads=threads, 
                                    MERANGS_EXE=MERANGS_EXE,
                                    STAR_EXE=STAR_EXE,
                                    SAMBAMBA_EXE=sambamba_cmd,
                                    force=overwrite )
                            postfilter_bam(f_bam, final_bam)
                        else:
                            logging.warn("Will not re-create existing file %s" % (f_bam))
                    if overwrite or not files_exist(f_bam_conv):
                        runMERANGS(f_bam_conv, 
                                    merangs_genome_idx, 
                                    f_fq_conv, 
                                    gtf=merangs_splice_gtf,
                                    threads=threads, 
                                    MERANGS_EXE=MERANGS_EXE,
                                    STAR_EXE=STAR_EXE,
                                    SAMBAMBA_EXE=sambamba_cmd,
                                    force=overwrite )
                        postfilter_bam(f_bam_conv, final_bam_conv) 
                    else:
                        logging.warn("Will not re-create existing file %s" % (f_bam_conv))      
                else:
                    logging.warn("Unknown mapper %s configured" % (mapper))
                    

                if files_exist([final_bam]):
                    bams[cond.id + "."+mapper]=os.path.abspath(final_bam)
                if files_exist([final_bam_conv]):
                    bams[cond.id + "."+mapper]=os.path.abspath(final_bam_conv)
                    
    # create TDF files
    if 'create_tdf' in config and config['create_tdf']:
        logging.info("Creating TDF files for %i bams" % (len(bams)))
        igvtools_log=tmpdir + config['dataset_name'] + ".igvtools.log" 
        for b in list(bams.values()):
            if overwrite or not files_exist(b+".tdf"):
                create_tdf(b, b+".tdf", m.chrom_sizes, igvtools_cmd=igvtools_cmd, logfile=igvtools_log)
    
    logging.info("All done in %s" % str(datetime.timedelta(seconds=time.time()-startTime)))
    print("All done in", datetime.timedelta(seconds=time.time()-startTime))