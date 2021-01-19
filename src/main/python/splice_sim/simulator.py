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
from utils import *
import random
import math
import gzip
import logging
from subprocess import check_output
from model import Model, Condition, Isoform, Transcript
from config_creator import calculate_transcript_data

def check_config(config):
    for section in ['dataset_name','mappers','genome_fa','gene_gff','conditions','isoform_mode','transcript_data']:
        assert section in config, "Missing configuration section %s" % section
        
        
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


def runHISAT2_TLA(bam, ref, fq, idx1, idx2, known_splicesites=None, force=True, run_flagstat=False, threads=1, doSort=True, additionalParameters=[], HISAT2_EXE='hisat2', SAMBAMBA_EXE='sambamba'):
    """ Run HISAT2_TLA """
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
         "--base-change", "TC",
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
    cmd=[igvtools_cmd, "count", "-w", "1", "-z", "7", "--strands", "read", bam, tdf, chrom_sizes]
    return pipelineStep(bam, tdf, cmd, shell=True, stdout=logfile)

#============================================================================
#    Classes
#============================================================================

class Stat():
    """ For collecting mapper-specific stats """
    def __init__(self, id, value, cond=None, mapper=None ):
        self.id=id
        self.cond=cond
        self.mapper=mapper
        self.value=value
        
    def get_header(self):
        ret = ("id\tcond\tmapper\tvalue")
        return (ret)
    
    def __str__(self):
        ret = ("%s\t%s\t%s\t%s" % (self.id if self.id is not None else "NA", 
                                   self.cond if self.cond is not None else "NA",
                                   self.mapper if self.mapper is not None else "NA",
                                   str(self.value)) )
        return (ret)  
    


def postfilter_bam( bam_in, bam_out, tag_tc="xc", tag_mp="xp" ):
    """ set proper read names and bam tags """  
    if tag_tc is None:
        tag_tc="xc"
    if tag_mp is None:
        tag_mp="xp"    
    samin = pysam.AlignmentFile(bam_in, "rb")
    mapped_out = pysam.AlignmentFile(bam_out, "wb", template=samin )
    n_reads=0
    f_reads=0
    for read in samin.fetch(until_eof=True): # NOTE important to add until_eof=True otherwise unmapped reads will be skipped.
        if read.is_secondary or read.is_supplementary:
            f_reads+=1
            continue
        #print(read.is_reverse)
        # parse from read name
        read_name = read.query_name
        # example    : ENSMUST00000100497.10_-_mat_2-151_5_1-2,3-4,5-6_142905567,142905452
        # example (tc): ENSMUST00000100497.10_-_mat_3-253_5_1-2,3-4,5-6_142905629_72,58,9
        is_tc_read = read_name.count('_')==7
        if not is_tc_read:
            read_name+='_NA' # no tc tag!
            
        true_tid,true_strand,true_isoform,tag,true_chrom,true_read_cigar,true_seqerr,tc_pos = read_name.split("_")
        #true_seqerr=true_seqerr.split(',') if true_seqerr != 'NA' else None        
        # calc absolute t/c conversions
        valid_tc=0
        if tc_pos != 'NA':
            tc_pos = [str(x) for x in tc_pos.split(',') if x is not None and x !=''] # t/c positions, relative coords of simulated read
            valid_tc=len(tc_pos)
            tc_pos = ','.join(tc_pos) if len(tc_pos)>0 else 'NA' # drop invalid (softclipped) t/c conversions
        is_correct_strand = ( read.is_reverse and true_strand=='-' ) or ((not read.is_reverse) and true_strand=='+') 
        
        read.query_name = '_'.join([true_tid,true_strand,true_isoform,tag,true_chrom,true_read_cigar,true_seqerr,tc_pos])

        if is_tc_read:
            if valid_tc==0: # no valid TC conversion 
                read.set_tag(tag=tag_tc, value=0, value_type="i")
                # set read color to light grey!
                read.set_tag(tag='YC', value='200,200,200' if is_correct_strand else '255,0,0', value_type="Z")
            else:
                read.set_tag(tag=tag_tc, value=valid_tc, value_type="i")
                # set blue read color intensity relative to number of tc conversions!
                intensity = 255/min(valid_tc, 25) if valid_tc > 0 else 255
                read.set_tag(tag='YC', value='%i,%i,255' % (intensity,intensity) if is_correct_strand else '255,0,0', value_type="Z")
        mapped_out.write(read)  
        n_reads=n_reads+1
    samin.close()
    mapped_out.close()
    
    pysam.index(bam_out)
    return True, n_reads, f_reads


def create_cigatuples(rtuples):
    """ convert aligned block tuples to pysam cigartuples """
    cigar_tuples=[]
    last=None
    for t in rtuples:
        if last  is not None:
            cigar_tuples+=[(3, t[0]-last-1)] # N-block
        cigar_tuples+=[(0, t[1]-t[0]+1)] # M-block
        last=t[1]
    return cigar_tuples

def fastq_to_bam(fastq_file, m, bam_file, tag_tc='xc', SAMBAMBA_EXE='sambamba', threads=1):
    """ Convert simulated FASTQ files to BAM """
    if tag_tc is None:
        tag_tc="xc"
    chromosomes = m.genome.references
    dict_chr2idx = {k: v for v, k in enumerate(chromosomes)}
    dict_idx2chr = {v: k for v, k in enumerate(chromosomes)}
    dict_chr2len = {c: m.genome.get_reference_length(c) for c in m.genome.references}
    header = { 'HD': {'VN': '1.4'}, 
               'SQ': [{'LN': l, 'SN': c} for c,l in dict_chr2len.items()] }
    with pysam.AlignmentFile(bam_file+'.tmp.bam', "wb", header=header) as bam_out:
        with gzip.open(fastq_file, 'rt') as file_in:
            for l1 in file_in:
                l2=str(next(file_in)).rstrip()
                l3=str(next(file_in)).rstrip()
                l4=str(next(file_in)).rstrip()
                read_name=l1.rstrip()[1:]
                is_tc_read = read_name.count('_')==7
                if not is_tc_read:
                    read_name+='_NA' # no tc tag!
                true_tid,true_strand,true_isoform,tag,true_chrom,true_read_cigar,true_seqerr,tc_pos = read_name.split("_")
                n_true_seqerr = len(true_seqerr.split(',')) if true_seqerr != 'NA' else 0
                n_tc_pos = len(tc_pos.split(',')) if tc_pos != 'NA' else 0
                read = pysam.AlignedSegment()
                read.query_name = read_name
                read.query_sequence = l2 if true_strand=='+' else reverse_complement(l2)
                read.query_qualities = pysam.qualitystring_to_array(l4)
                true_tuples = [tuple([int(y) for y in x.split('-')]) for x in true_read_cigar.split(',')]
                read.reference_start = true_tuples[0][0]-1
                read.cigartuples=create_cigatuples(true_tuples)
                read.reference_id = dict_chr2idx[true_chrom]
                read.mapping_quality = 60
                read.tags = (("NM", n_true_seqerr + n_tc_pos ),(tag_tc, n_tc_pos))
                read.flag = 0 if true_strand=='+' else 16
                bam_out.write(read)
    bam_out.close()
    sambambasortbam(bam_file+'.tmp.bam', 
                    bam_file, 
                    index=True, 
                    override=True, 
                    delinFile=True,
                    ncpu=threads, 
                    EXE=SAMBAMBA_EXE)
    

# config = json.load(open('/Volumes/groups/ameres/Niko/projects/Ameres/splicing/splice_sim/testruns/small4/config.json'), object_pairs_hook=OrderedDict)
# config["transcripts"] = json.load(open('/Volumes/groups/ameres/Niko/projects/Ameres/splicing/splice_sim/testruns/small4/' + config['transcript_data']), object_pairs_hook=OrderedDict)
# for k in ['gene_gff', 'genome_fa']:
#     config[k] = '/Volumes' + config[k] 
# m = Model(config)
# fq_file = '/Volumes/groups/ameres/Niko/projects/Ameres/splicing/splice_sim/testruns/small4/small4/sim/tmp/small4.5min.TC.fq.gz'
# bam_file= '/Volumes/groups/ameres/Niko/projects/Ameres/splicing/splice_sim/testruns/small4/delme.bam'
# fastq_to_bam(fq_file, m, bam_file)



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
    stats=[Stat("transcripts", len(m.transcripts))]
    
    # now write one fasta file per transcript/cond
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
    simulated_read_stats=Counter()
    incosistent_data=0
    for cond in m.conditions:
        f = tmpdir + config['dataset_name'] + "." + cond.id + ".fa"
        # art output
        art_out_prefix = tmpdir + config['dataset_name'] + "." + cond.id + ".bothstrands"
        f_fq_both = art_out_prefix + ".fq"
        f_sam_both =  art_out_prefix + ".sam"
        # filtered final output
        f_fq =  tmpdir + config['dataset_name'] + "." + cond.id + ".fq"
        if overwrite or not files_exist([f_fq+".gz" ]):
            # NOTE: use 2x coverage as ~50% of reads will be simulated for wrong strand and will be dropped in postprocessing
            cmd=[art_cmd, "-ss", "HS25", "-i", f, "-l", str(config["readlen"]), "-f", str(cond.coverage * 2), "-na", "--samout", "-o", art_out_prefix ]
            if "random_seed" in config:
                cmd+=["-rs", str(config["random_seed"])]
            success = pipelineStep(f, [f_fq_both, f_sam_both], cmd, shell=True, stdout=artlog, append=True)    
            if not success:
                logging.error('error simulating %s' % (f_fq_both))
    
            # filter reads that map to wrong strand and extend read names
            # get alignment positions from art aln file. Note that seq-read errors are encoded in the 
            # CIGAR string (e.g., '16=1X4=1X78=') means that bases 17 and 22 are mismatches wrt. reference.
            # FIXME: reads with INDELs are currently skipped! 
            # deletion example: 80=1D20=
            # insertion example: 45=1X47=1I6=
            # TODO: 1) provide proper cigarblocks for INDEL reads
            #       2) provide seqerr positions for INDELs 
            with open(f_fq, 'w') as out_fq:
                sam = pysam.AlignmentFile(f_sam_both, "rb")
                for r in sam.fetch(until_eof=True):
                    # e.g. 'ENSMUST00000099151.5_+_mat_0-400'. tag ('0-400') consists of running number from 'abundance' and running number form ART
                    tid,transcript_strand,iso_id,tag=r.query_name.split("_")
                    read_strand = "-" if r.is_reverse else "+"
                    if transcript_strand == read_strand: # NOTE: reads that were mapped to opposite strand will be dropped later in post_filtering step!
                        iso = m.transcripts[tid].isoforms[iso_id]
                        start_abs, bid_start = iso.rel2abs_pos(r.reference_start)
                        end_abs, bid_end = iso.rel2abs_pos(r.reference_end-1)
                        if read_strand == '-': # swap start/end coords
                            start_abs, bid_start, end_abs, bid_end = end_abs, bid_end, start_abs, bid_start 
                            
                        
                        block_tuples, seqerr_pos = parse_art_cigar(r)
                        read_cigartuples=[]
                        read_cigartuples_len=0
                        read_spliced=False
                        for (x,y) in block_tuples:
                            start_abs, bid_start = iso.rel2abs_pos(x)
                            end_abs, bid_end = iso.rel2abs_pos(y)
                            if not read_spliced and bid_start != bid_end:
                                read_spliced=True
                            rc, rl = iso.calc_cigartuples(start_abs, end_abs)
                            read_cigartuples+=rc
                            read_cigartuples_len+=rl
                        read_cigar  = ','.join(["%i-%i" % (a,b) for (a,b) in read_cigartuples]) # to parse cigar_tuples = [tuple(x.split('-')) for x in read_cigar.split(',')]
                        read_name = "%s_%s_%s_%s" % (r.query_name,                                                             
                                                        iso.t.transcript.Chromosome, 
                                                        read_cigar, 
                                                        (",".join(str(iso.rel2abs_pos(seqerr_pos)[0]) for p in rel_pos) )  if len(rel_pos)>0 else 'NA' )
                        qstr = ''.join(map(lambda x: chr( x+33 ), r.query_qualities))
                        # double check whether seq length and calculates CIGAR length match and if not skip read.
                        # FIXME: why does this happen?
                        if ( read_cigartuples_len != len(r.query_sequence)):
                            print("Data inconsistent for read %s with cigar %s. Ignoring..." % (read_name, r.cigarstring))
                            incosistent_data+=1
                        else:
                            simulated_read_stats[cond.id,tid]+=1
                            print("@%s\n%s\n+\n%s" % ( read_name, r.query_sequence, qstr), file=out_fq )
            if success:
                removeFile([f_fq_both, f_sam_both])
    
        else:
            logging.warn("Will not re-create existing file %s" % (f_fq+".gz"))
    if incosistent_data > 0:
        logging.warn("Dropped %i simulated reads due to inconsistent data" % incosistent_data)
    
    # write simulated reead stats
    f_srs = outdir + 'simulated_read_stats.tsv'
    with open(f_srs, 'w') as out:
        print("cond\ttranscript_id\ttlen\tn_reads\tcov", file=out)
        for k in simulated_read_stats.keys():
            condid,transcript_id = k
            tlen=m.transcripts[tid].len
            n_reads=simulated_read_stats[k]
            cov=n_reads*config['readlen']/tlen
            print("%s\t%s\t%i\t%i\t%f.2" % (condid, transcript_id, tlen, n_reads, cov), file=out)
            
    # concat files per condition, introduce T/C conversions and bgzip
    for cond in m.conditions:
        f_all = tmpdir + config['dataset_name'] + "." + cond.id + ".fq"
        f_tc  = tmpdir + config['dataset_name'] + "." + cond.id + ".TC.fq"
        hist_tc=[]
        if overwrite or not files_exist(f_tc+".gz"):
            with open(f_tc, 'w') as out:
                with open(f_all, 'r') as fin:
                    buf=[]
                    lines = []
                    for line in fin:
                        lines.append(line.rstrip())
                        if len(lines) == 4:
                            names=lines[0].split("_")
                            seq, conv=add_tc(lines[1], names[1], cond.conversion_rate)
                            hist_tc+=[len(conv)]
                            buf+=["_".join(names)+'_'+",".join(str(x) for x in conv)+"\n"]
                            buf+=[seq+"\n"]
                            buf+=[lines[2]+"\n"]
                            buf+=[lines[3]+"\n"]
                            if len(buf)>10000:
                                out.writelines(buf)
                                buf=[]
                            lines = []
                    if len(buf)>0:
                        out.writelines(buf)
            if not write_uncoverted:
                removeFile(f_all)
            else:
                bgzip(f_all, override=True, delinFile=True, threads=threads)
            bgzip(f_tc, override=True, delinFile=True, threads=threads)
        else:
            logging.warn("Will not re-create existing file %s" % (f_tc+".gz"))
        stats+=[Stat("mean_tc", np.mean(hist_tc),cond=cond)]
        stats+=[Stat("median_tc", np.median(hist_tc),cond=cond)]
        stats+=[Stat("hist_tc", ",".join(str(x) for x in Counter(hist_tc).keys()) + ":" + ",".join(str(x) for x in Counter(hist_tc).values()),cond=cond)]
            
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
            if mapper == "STAR":
                EXE=config['mappers'][mapper]['star_cmd'] if 'star_cmd' in config['mappers'][mapper] else 'STAR'
            elif mapper == "HISAT2_TLA":
                EXE=config['mappers'][mapper]['hisat2_cmd'] if 'hisat2_cmd' in config['mappers'][mapper] else 'hisat2'
            logging.info("%s version: %s" % ( mapper, get_version(EXE)) )

        for cond in m.conditions:
            for mapper in config['mappers'].keys():
                f_all = tmpdir + config['dataset_name'] + "." + cond.id + ".fq.gz"
                b_all = tmpdir + config['dataset_name'] + "." + cond.id + "."+mapper+".bam"
                f_tc  = tmpdir + config['dataset_name'] + "." + cond.id + ".TC.fq.gz"
                b_tc  = tmpdir + config['dataset_name'] + "." + cond.id + "."+mapper+".TC.bam"
                bamdir_all = outdir + "bam_ori/" + mapper + "/"
                bamdir_tc  = outdir + "bam_tc/" + mapper + "/"
                if not os.path.exists(bamdir_all):
                    os.makedirs(bamdir_all)
                if not os.path.exists(bamdir_tc):
                    os.makedirs(bamdir_tc)
                final_all  = bamdir_all + config['dataset_name'] + "." + cond.id + "."+mapper+".bam"
                final_tc   = bamdir_tc  + config['dataset_name'] + "." + cond.id + "."+mapper+".TC.bam"
                tag_tc = config['mappers'][mapper]['tag_tc'] if 'tag_tc' in config['mappers'][mapper] else None
                tag_mp = config['mappers'][mapper]['tag_mp'] if 'tag_mp' in config['mappers'][mapper] else None
                
                logging.info("Mapping %s reads with %s" % (cond.id, mapper))
                if mapper == "STAR":
                    # run STAR mapper 
                    STAR_EXE=config['mappers'][mapper]['star_cmd'] if 'star_cmd' in config['mappers'][mapper] else 'STAR'
                    star_splice_gtf=config['mappers'][mapper]['star_splice_gtf'] if 'star_splice_gtf' in config['mappers'][mapper] else None
                    star_genome_idx=config['mappers'][mapper]['star_genome_idx']
                    if write_uncoverted:
                        if overwrite or not files_exist(final_all):
                            runSTAR(b_all, 
                                    star_genome_idx, 
                                    reads1=[f_all], 
                                    gtf=star_splice_gtf,
                                    threads=threads, 
                                    STAR_EXE=STAR_EXE,
                                    SAMBAMBA_EXE=sambamba_cmd,
                                    force=overwrite )
                        else:
                            logging.warn("Will not re-create existing file %s" % (b_all))
                    if overwrite or not files_exist(final_tc):
                        runSTAR(b_tc, 
                            star_genome_idx, 
                            reads1=[f_tc], 
                            gtf=star_splice_gtf,
                            threads=threads, 
                            STAR_EXE=STAR_EXE,
                            SAMBAMBA_EXE=sambamba_cmd,
                            force=overwrite )   
                    else:
                        logging.warn("Will not re-create existing file %s" % (b_tc))            
                elif mapper == "HISAT2_TLA":
                    HISAT2_EXE=config['mappers'][mapper]['hisat2_cmd'] if 'hisat2_cmd' in config['mappers'][mapper] else 'hisat2'
                    hisat2_idx1=config['mappers'][mapper]['hisat2_idx1']
                    hisat2_idx2=config['mappers'][mapper]['hisat2_idx2']
                    hisat2_kss=config['mappers'][mapper]['hisat2_kss'] if 'hisat2_kss' in config['mappers'][mapper] else None
                    if write_uncoverted:
                        if overwrite or not files_exist(final_all):
                            runHISAT2_TLA(b_all, 
                                    config["genome_fa"], 
                                    fq=f_all,
                                    idx1=hisat2_idx1,
                                    idx2=hisat2_idx2,
                                    known_splicesites=hisat2_kss,
                                    threads=threads, 
                                    HISAT2_EXE=HISAT2_EXE,
                                    SAMBAMBA_EXE=sambamba_cmd,
                                    force=overwrite )
                        else:
                            logging.warn("Will not re-create existing file %s" % (b_all))
                    if overwrite or not files_exist(final_tc):
                        runHISAT2_TLA(b_tc, 
                            config["genome_fa"], 
                            fq=f_tc,
                            idx1=hisat2_idx1,
                            idx2=hisat2_idx2,
                            known_splicesites=hisat2_kss,
                            threads=threads, 
                            HISAT2_EXE=HISAT2_EXE,
                            SAMBAMBA_EXE=sambamba_cmd,
                            force=overwrite )
                    else:
                        logging.warn("Will not re-create existing file %s" % (b_tc))
                else:
                    logging.warn("Unknown mapper %s configured" % (mapper))
                    
                # Now filter all reads that mapped to the wrong strand as ART has no strand support and write all TC mutation
                if files_exist([b_all]):
                    success,n_reads, f_reads = postfilter_bam( b_all, final_all, tag_tc, tag_mp)
                    if success:
                        removeFile([b_all, b_all+".bai"]) # remove unfiltered/original BAM
                        stats+=[Stat("all_reads", n_reads,cond=cond.id, mapper=mapper)]
                        stats+=[Stat("all_filtered_secondary_alignments", f_reads,cond=cond.id, mapper=mapper)]
    
                if files_exist([b_tc]):
                    success,n_reads,f_reads = postfilter_bam( b_tc, final_tc, tag_tc, tag_mp)
                    if success:
                        removeFile([b_tc, b_tc+".bai"]) # remove unfiltered/original BAM
                        stats+=[Stat("tc_reads", n_reads,cond=cond.id, mapper=mapper)]
                        stats+=[Stat("tc_filtered_secondary_alignments", f_reads,cond=cond.id, mapper=mapper)]
    
                # store bams for tdf creation
                if files_exist([final_all]):
                    bams[cond.id + "."+mapper]=os.path.abspath(final_all)
                if files_exist([final_tc]):
                    bams[cond.id + "."+mapper]=os.path.abspath(final_tc)
                    
            # create truth BAMs
            bamdir_all = outdir + "bam_ori/TRUTH/"
            if not os.path.exists(bamdir_all):
                os.makedirs(bamdir_all)
            bamdir_tc  = outdir + "bam_tc/TRUTH/"
            if not os.path.exists(bamdir_tc):
                os.makedirs(bamdir_tc)
            final_all_truth  = bamdir_all + config['dataset_name'] + "." + cond.id + ".TRUTH.bam"
            final_tc_truth   = bamdir_tc  + config['dataset_name'] + "." + cond.id + ".TRUTH.TC.bam"
            if not files_exist([final_all_truth, final_tc_truth]):
                fastq_to_bam(f_all, m, final_all_truth, tag_tc=tag_tc, threads=threads)
                fastq_to_bam(f_tc, m, final_tc_truth, tag_tc=tag_tc, threads=threads)
    
    
    # write stats. FIXME: stats are not complete if pipeline was restarted with 
    logging.info("Writing stats")
    f_roi = outdir + "stats.tsv"
    with open(f_roi, 'w') as out:
        print(Stat("",1).get_header(), file=out)
        for s in stats:
            print(s, file=out)
            
    
    f_roi = outdir + "stats.tsv"
    with open(f_roi, 'w') as out:
        print(Stat("",1).get_header(), file=out)
        for s in stats:
            print(s, file=out)
    
#     # write slamstr config
#     if 'slamstr_config_template' in config and 'mappers' in config:
#         logging.info("Writing slamstr config files")
#         with open(config['slamstr_config_template'], 'r') as f:
#             x = f.read().splitlines()
#         for mapper in config['mappers'].keys():
#             modi=['tc','ori'] if write_uncoverted else ['tc']
#             for mode in modi:
#                 id = "slamstr."+mapper+"."+mode
#                 f_config = outdir + id+".json"
#                 with open(f_config, 'w') as out:
#                     tokens={}
#                     tokens["ID"]=id
#                     tokens["genome_fa"]=os.path.abspath(config['genome_fa'])
#                     tokens["gff"]=os.path.abspath(f_anno)
#                     tokens["datasets"]="\n"
#                     tokens["timepoints"]="\n"
#                     for idx, cond in enumerate(m.conditions):
#                         tokens["datasets"]+='\t"' + cond.id + '": "' + bams[cond.id + "."+mapper] + '"' + ("," if idx < len(m.conditions)-1 else "") + "\n"
#                         tokens["timepoints"]+='\t"' + cond.id + '": "' + str(cond.timepoint)+ '"' + ("," if idx < len(m.conditions)-1 else "")+ "\n"
#                     for l in x:
#                         print(replace_tokens(l, tokens), file=out)
                        
    # create TDF files
    if 'create_tdf' in config and config['create_tdf']:
        logging.info("Creating TDF files for %i bams" % (len(bams)))
        igvtools_log=tmpdir + config['dataset_name'] + ".igvtools.log" 
        for b in list(bams.values()):
            if overwrite or not files_exist(b+".tdf"):
                create_tdf(b, b+".tdf", m.chrom_sizes, igvtools_cmd=igvtools_cmd, logfile=igvtools_log)
    
    logging.info("All done in %s" % str(datetime.timedelta(seconds=time.time()-startTime)))
    print("All done in", datetime.timedelta(seconds=time.time()-startTime))