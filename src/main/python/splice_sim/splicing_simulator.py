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
from model import Condition, Isoform, Transcript

#============================================================================
# splice sim readme
#
# NOTE: splice_sim will not simulate reads for annotations that are shorter than the readlength (as this is 
# not supported by ART.
#============================================================================
VERSION = "0.1"
logo = '''
===================================
Slamstr splice_sim v%s
===================================
''' % VERSION


#============================================================================
# Mapping methods
#============================================================================


#
# Run STAR, supports SE+PE and can merge multiple files (needs sambamba access). Can create samtools flagstats
#
def runSTAR(bam, ref, reads1=[], reads2=None,  gtf=None, force=True, run_flagstat=False, threads=1, chimeric=0, writeUnmapped=False, doSort=True, additionalParameters=[], STAR_EXE='STAR'):
    success = True
    if files_exist(bam) and not force:
        print("BAM file " + bam + " already exists! use -f to recreate.")
        return True
    # check input params
    if len(reads1)==0 or (reads2 is not None and len(reads1)!=len(reads2)):
        print("ERROR")
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
            cmd+=["--outReadsUnmapped", "Fastx"]
        if gtf:
            cmd+=["--sjdbGTFfile", gtf]
        cmd += additionalParameters
        print(cmd)
        success = success and pipelineStep([reads1,reads2], starout, cmd, shell=True) 
        print(starout)
        success = success and sambamba2bam(starout, sorted,     
                                          sort=doSort, 
                                          index=True, 
                                          delinFile=True,
                                          ncpu=threads, EXE=sambamba_cmd) 
        tomerge+=[sorted]
        # handle chimeric reads if configured
        if files_exist(starout_chimeric):
            success = success and sambamba2bam(starout_chimeric, sorted_chimeric, 
                                          sort=doSort, 
                                          index=True, 
                                          delinFile=True,
                                          ncpu=threads, EXE=sambamba_cmd) 
            tomerge_chimeric+=[sorted_chimeric]
    print(tomerge)
    # merge files
    if len(tomerge)==1:
        os.rename(tomerge[0], bam)
        os.rename(tomerge[0]+".bai", bam+".bai")
    else:
        logging.info("Merging from multiple input files")
        success = success and sambambamerge(bam, tomerge, ncpu=threads, EXE=sambamba_cmd)
    if len(tomerge_chimeric)==0:
        logging.info("NOTE: no chimeric read output configured!")
    else:
        if len(tomerge_chimeric)==1:
            os.rename(tomerge_chimeric[0], bam_chimeric)
            os.rename(tomerge_chimeric[0]+".bai", bam_chimeric+".bai")
        else:
            logging.info("Merging from multiple input files")
            success = success and sambambamerge(bam_chimeric, tomerge_chimeric, ncpu=threads, EXE=sambamba_cmd)
        flagstat( bam_chimeric )
    if run_flagstat:            
        flagstat( bam )        
    return success


#
# Run HISAT2_TLA
#
def runHISAT2_TLA(bam, ref, fq, idx1, idx2, known_splicesites=None, force=True, run_flagstat=False, threads=1, doSort=True, additionalParameters=[], HISAT2_EXE='hisat2'):
    success = True
    if files_exist(bam) and not force:
        print("BAM file " + bam + " already exists! use -f to recreate.")
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
    print(cmd)
    success = success and pipelineStep([fq], sam, cmd, shell=True)
    
    # sort + index
    success = success and sambamba2bam(sam, bam, 
                              sort=doSort, 
                              index=True, 
                              delinFile=True,
                              ncpu=threads) 
    
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
    


#============================================================================
def postfilter_bam( bam_in, bam_out, tag_tc=None, tag_mp=None):
    samin = pysam.AlignmentFile(bam_in, "rb")
    samout = pysam.AlignmentFile(bam_out, "wb", template=samin, )
    if tag_tc is None:
        tag_tc="xc" 
    if tag_mp is None:
        tag_mp="xp" 
    n_reads=0
    for read in samin.fetch():
        #print(read.is_reverse)
        comp=read.query_name.split("_")
        is_correct_strand = ( read.is_reverse and comp[1]=='-' ) or ((not read.is_reverse) and comp[1]=='+') 
        read.query_name = "_".join(comp[0:4]) # cut of _tc section from read name.
        if len(comp)>4:
            tc_pos = comp[4][3:].split(",")
            if not tc_pos[0]: # no TC conversion
                read.set_tag(tag=tag_tc, value=0, value_type="i")
                # set read color to light grey!
                read.set_tag(tag='YC', value='200,200,200' if is_correct_strand else '255,0,0', value_type="Z")
            else:
                slamdunk_type=2 if read.is_reverse else 16 # see paper supplement
                valid_tc=0
                mp_str=[]
                for x in tc_pos:
                    rpos = int(x)
                    gpos = readidx2genpos_rel(read, rpos)
                    if gpos is not None:
                        valid_tc+=1
                        mp_str+=[str(slamdunk_type)+":"+str(rpos)+":"+str(gpos)]
                read.set_tag(tag=tag_tc, value=valid_tc, value_type="i")
                # set blue read color intensity relative to number of tc conversions!
                intensity = 255/min(valid_tc, 25) if valid_tc > 0 else 255
                read.set_tag(tag='YC', value='%i,%i,255' % (intensity,intensity) if is_correct_strand else '255,0,0', value_type="Z")

                if valid_tc>0:
                    read.set_tag(tag=tag_mp, value=",".join(mp_str), value_type="Z")
        samout.write(read)  
        n_reads=n_reads+1

    samin.close()
    samout.close()
    pysam.index(bam_out)
    return True, n_reads

if __name__ == '__main__':
    #============================================================================
    usage = '''                           
    
      Copyright 2019 Niko Popitsch. All rights reserved.
      
      Licensed under the Apache License 2.0
      http://www.apache.org/licenses/LICENSE-2.0
      
      Distributed on an "AS IS" basis without warranties
      or conditions of any kind, either express or implied.
    
    USAGE
    '''
    
    
    parser = ArgumentParser(description=usage, formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument("-c","--config", type=existing_file, required=True, dest="confF", help="JSON config file")
    parser.add_argument("-f","--force", required=False, action='store_true', dest="force", help="If set, existing results files will be overwritten")
    parser.add_argument("-o","--out", type=str, required=False, dest="outdir", metavar="outdir", help="Output folder")
    args = parser.parse_args()   
    startTime = time.time()
    
    print(logo)
    
    # output dir (current dir if none provided)
    outdir = args.outdir if args.outdir else os.getcwd()
    if not outdir.endswith("/"):
        outdir += "/"
    if not os.path.exists(outdir):
        print("Creating dir " + outdir)
        os.makedirs(outdir)
    tmpdir = outdir + "/tmp/"
    if not os.path.exists(tmpdir):
        os.makedirs(tmpdir)
    
    # load + check config
    #config = json.load(open('/Users/niko.popitsch/eclipse-workspace/slamstr/src/test/splicing_simulator/testconfig.json'), object_pairs_hook=OrderedDict)
    config = json.load(open(args.confF), object_pairs_hook=OrderedDict)
    if "random_seed" in config:
        random.seed(config["random_seed"])
    write_uncoverted=config["write_uncoverted"] if "write_uncoverted" in config else True
    
    # executables
    art_cmd =      config["executable"]["art_cmd"] if "executable" in config and "art_cmd" in config["executable"] else 'art_illumina'
    sambamba_cmd = config["executable"]["sambamba_cmd"] if "executable" in config and "sambamba_cmd" in config["executable"] else 'sambamba'
    igvtools_cmd = config["executable"]["igvtools_cmd"] if "executable" in config and "igvtools_cmd" in config["executable"] else 'igvtools'
    threads=config["threads"] if "threads" in config else 1
    
    
    # get conditions and configuration
    conditions=[]
    for id in config["conditions"].keys():
        conditions+=[Condition(id, config["conditions"][id][0], config["conditions"][id][1], config["conditions"][id][2] )]
    print("Configured conditions: %s" % (", ".join(str(x) for x in conditions)) )
    
    if args.force:
        print("NOTE: existing result files will be overwritten!")
    
    # read transcript data from external file if not in config
    if 'transcripts' not in config:
        assert "transcript_data" in config, "Transcript data needs to be configured either in config file ('transcripts' section) or in an external file referenced via 'transcript_data'"
        tfile = config['transcript_data']
        if not os.path.isabs(tfile):
            tfile = os.path.dirname(os.path.abspath(args.confF))+"/"+tfile
        tdata = json.load(open(tfile), object_pairs_hook=OrderedDict)
        config["transcripts"]=tdata
    
    
    
    
    # load + filter gene gff
    print("Loading gene GFF")
    gff = pr.read_gff3(config["gene_gff"])
    gff_df = gff.df
    gff_df.Start = gff_df.Start + 1 # correct for pyranges bug?
    df = gff_df[gff_df['transcript_id'].isin(list(config['transcripts'].keys()))] # reduce to contain only configured transcripts
    
    print("Writing filtered gene GFF")
    f_anno=outdir+"gene_anno.gff3"
    d=pd.read_csv(config["gene_gff"],delimiter='\t',encoding='utf-8')
    if not files_exist(f_anno+".gz"):
        with open(f_anno, 'w') as out:
            for index, row in d.iterrows():
                keep=False
                for k in row[8].split(';'):
                    if k.startswith('transcript_id='):
                        keep=k[len('transcript_id='):] in config['transcripts'].keys()
                if keep:
                    print('\t'.join(str(x) for x in row), file=out)
        bgzip(f_anno, override=True, delinFile=True, index=True, threads=threads)
    f_anno=f_anno+".gz"
          
    # load genome
    print("Loading genome")
    genome = pysam.FastaFile(config["genome_fa"])
    chrom_sizes = config["chrom_sizes"] if 'chrom_sizes' in config else config["genome_fa"]+".chrom.sizes"
    
    # instantiate transcripts
    transcripts=OrderedDict()
    max_ilen = config["max_ilen"] if 'max_ilen' in config else None
    stats=[]
    for tid in list(config['transcripts'].keys()):
        tid_df = df[df['transcript_id']==tid] # extract data for this transcript only
        if tid_df.empty:
            print("No annotation data found for configured tid %s, skipping..." % (tid))
        else:
            t = Transcript(tid, tid_df, genome, conditions, max_ilen=max_ilen) 
            if t.is_valid:
                transcripts[tid] = t
    stats+=[Stat("transcripts", len(transcripts))]
    
    # now write one fasta file per transcript/cond
    print("Calculating isoform data")
    for cond in conditions:
        fout = tmpdir + config['dataset_name'] + "." + cond.id + ".fa"
        if args.force or (not files_exist(fout) and not files_exist(fout+".gz")):
            buf=[]
            with open(fout, 'w') as out:
                for t in transcripts.values():
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
            print("Will not re-create existing file %s" % (fout))
    
    # run art simulator 
    print("Running read simulator")
    artlog=tmpdir + config['dataset_name'] + ".art.log" 
    valid_reads={}
    for cond in conditions:
        valid_reads[cond]=set()
        f = tmpdir + config['dataset_name'] + "." + cond.id + ".fa"
        # art output
        art_out_prefix = tmpdir + config['dataset_name'] + "." + cond.id + ".bothstrands"
        f_fq_both = art_out_prefix + ".fq"
        f_sam_both =  art_out_prefix + ".sam"
        # filtered final output
        f_truth = tmpdir + config['dataset_name'] + "." + cond.id + ".truth_tsv"
        f_fq =  tmpdir + config['dataset_name'] + "." + cond.id + ".fq"
        if args.force or not files_exist([f_fq+".gz", f_truth+".gz" ]):
            # NOTE: use 2x coverage as ~50% of reads will be simulated for wrong strand and will be dropped in postprocessing
            cmd=[art_cmd, "-ss", "HS25", "-i", f, "-l", str(config["readlen"]), "-f", str(cond.coverage * 2), "-na", "--samout", "-o", art_out_prefix ]
            if "random_seed" in config:
                cmd+=["-rs", str(config["random_seed"])]
            success = pipelineStep(f, [f_fq_both, f_sam_both], cmd, shell=True, stdout=artlog, append=True)    
            if not success:
                print('error simulating %s' % (f_fq_both))
    
            # filter reads that map to wrong strand and write truth file
            # get alignment positions from art aln file. Note that seq-read errors are encoded in the 
            # CIGAR string (e.g., '16=1X4=1X78=') means that bases 17 and 22 are mismatches wrt. reference.
            # Q: why are there no D/I entries in the cigarstrings? How to find seqerr INDELs?
            with open(f_truth, 'w') as out_truth:
                with open(f_fq, 'w') as out_fq:
                    print("read_name\tstart_rel\tend_rel\tseq_err_pos_rel\tchr_abs\tstart_abs\tend_abs\tread_spliced\tseq_err_pos_abs", file=out_truth)
                    sam = pysam.AlignmentFile(f_sam_both, "rb")
                    read_stats={}
                    for r in sam.fetch():
                        # e.g. 'ENSMUST00000099151.5_+_mat_0-400'. tag ('0-400') consists of running number from 'abundance' and running number form ART
                        tid,transcript_strand,iso_id,tag=r.query_name.split("_")
                        read_strand = "-" if r.is_reverse else "+"
                        read_stats['all_'+read_strand]=read_stats['all_'+read_strand]+1 if 'all_'+read_strand in read_stats else 1
                        if transcript_strand == read_strand: # NOTE: reads that were mapped to opposite strand will be dropped later in post_filtering step!
                            valid_reads[cond].add(r.query_name)
                            iso = transcripts[tid].isoforms[iso_id]
                            start_abs, bid_start = iso.rel2abs_pos(r.reference_start)
                            end_abs, bid_end = iso.rel2abs_pos(r.reference_end-1)
                            if read_strand == '-': # swap start/end coords
                                start_abs, bid_start, end_abs, bid_end = end_abs, bid_end, start_abs, bid_start 
                            read_spliced = 1 if bid_start != bid_end else 0
                            rel_pos = cigar_to_rel_pos(r)
                            print("%s\t%i\t%i\t%s\t%s\t%i\t%i\t%i\t%s" % (r.query_name, 
                                                                          r.reference_start,
                                                                          r.reference_end, 
                                                                          (",".join(str(p) for p in rel_pos) ) if len(rel_pos)>0 else 'NA',
                                                                          iso.t.transcript.Chromosome,
                                                                          start_abs,
                                                                          end_abs, 
                                                                          read_spliced, 
                                                                          (",".join(str(iso.rel2abs_pos(p)[0]) for p in rel_pos) )  if len(rel_pos)>0 else 'NA'  
                                                                          ), file=out_truth )
                            qstr = ''.join(map(lambda x: chr( x+33 ), r.query_qualities))
                            print("@%s\n%s\n+\n%s" % ( r.query_name, r.query_sequence, qstr), file=out_fq )
                        else:
                            read_stats['dropped_'+read_strand]=read_stats['dropped_'+read_strand]+1 if 'dropped_'+read_strand in read_stats else 1
            if success:
                bgzip(f_truth, override=True, delinFile=True, threads=threads)
                removeFile([f_fq_both, f_sam_both])
            print(read_stats)
    
        else:
            print("Will not re-create existing file %s" % (f_fq+".gz"))
    
                   
    # concat files per condition, introduce T/C conversions and bgzip
    for cond in conditions:
        f_all = tmpdir + config['dataset_name'] + "." + cond.id + ".fq"
        f_tc  = tmpdir + config['dataset_name'] + "." + cond.id + ".TC.fq"
        hist_tc=[]
        if args.force or not files_exist(f_tc+".gz"):
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
                            buf+=["_".join(names)+"_tc:"+",".join(str(x) for x in conv)+"\n"]
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
            print("Will not re-create existing file %s" % (f_tc+".gz"))
        stats+=[Stat("mean_tc", np.mean(hist_tc),cond=cond)]
        stats+=[Stat("median_tc", np.median(hist_tc),cond=cond)]
        stats+=[Stat("hist_tc", ",".join(str(x) for x in Counter(hist_tc).keys()) + ":" + ",".join(str(x) for x in Counter(hist_tc).values()),cond=cond)]
            
    # write ROIs 
    print("Writing ROI bed")
    f_roi = outdir + "roi.bed"
    with open(f_roi, 'w') as out:
        for t in transcripts.values():
            print(t.to_bed(), file=out)
    
    
    # map reads if mappers configured
    bams={}
    if 'mappers' in config:
        print("Mapping reads")
        for cond in conditions:
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
                
                print("Mapping %s reads with %s" % (cond.id, mapper))
                if mapper == "STAR":
                    # run STAR mapper 
    
                    STAR_EXE=config['mappers'][mapper]['star_cmd'] if 'star_cmd' in config['mappers'][mapper] else 'STAR'
                    star_splice_gtf=config['mappers'][mapper]['splice_gtf'] if 'splice_gtf' in config['mappers'][mapper] else None
                    star_genome_idx=config['mappers'][mapper]['star_genome_idx']
                    if write_uncoverted:
                        if args.force or not files_exist(final_all):
                            runSTAR(b_all, 
                                    star_genome_idx, 
                                    reads1=[f_all], 
                                    gtf=star_splice_gtf,
                                    threads=threads, 
                                    STAR_EXE=STAR_EXE,
                                    force=args.force )
                        else:
                            print("Will not re-create existing file %s" % (b_all))
                    if args.force or not files_exist(final_tc):
                        runSTAR(b_tc, 
                            star_genome_idx, 
                            reads1=[f_tc], 
                            gtf=star_splice_gtf,
                            threads=threads, 
                            STAR_EXE=STAR_EXE,
                            force=args.force )   
                    else:
                        print("Will not re-create existing file %s" % (b_tc))            
                elif mapper == "HISAT2_TLA":
                    HISAT2_EXE=config['mappers'][mapper]['hisat2_cmd'] if 'hisat2_cmd' in config['mappers'][mapper] else 'hisat2'
                    hisat2_idx1=config['mappers'][mapper]['hisat2_idx1']
                    hisat2_idx2=config['mappers'][mapper]['hisat2_idx2']
                    hisat2_kss=config['mappers'][mapper]['hisat2_kss'] if 'hisat2_kss' in config['mappers'][mapper] else None
                    if write_uncoverted:
                        if args.force or not files_exist(final_all):
                            runHISAT2_TLA(b_all, 
                                    config["genome_fa"], 
                                    fq=f_all,
                                    idx1=hisat2_idx1,
                                    idx2=hisat2_idx2,
                                    known_splicesites=hisat2_kss,
                                    threads=threads, 
                                    HISAT2_EXE=HISAT2_EXE,
                                    force=args.force )
                        else:
                            print("Will not re-create existing file %s" % (b_all))
                    if args.force or not files_exist(final_tc):
                        runHISAT2_TLA(b_tc, 
                            config["genome_fa"], 
                            fq=f_tc,
                            idx1=hisat2_idx1,
                            idx2=hisat2_idx2,
                            known_splicesites=hisat2_kss,
                            threads=threads, 
                            HISAT2_EXE=HISAT2_EXE,
                            force=args.force )
                    else:
                        print("Will not re-create existing file %s" % (b_tc))
                else:
                    print("Unknown mapper %s configured" % (mapper))
                    
                # Now filter all reads that mapped to the wrong strand as ART has no strand support and write all TC mutation
                if files_exist([b_all]):
                    success,n_reads = postfilter_bam( b_all, final_all, tag_tc, tag_mp)
                    if success:
                        removeFile([b_all, b_all+".bai"])
                        stats+=[Stat("all_reads", n_reads,cond=cond, mapper=mapper)]
    
                if files_exist([b_tc]):
                    success,n_reads = postfilter_bam( b_tc, final_tc, tag_tc, tag_mp)
                    if success:
                        removeFile([b_tc, b_tc+".bai"])
                        stats+=[Stat("tc_reads", n_reads,cond=cond, mapper=mapper)]
    
                # store bams for tdf creation
                if files_exist([final_all]):
                    bams[cond.id + "."+mapper]=os.path.abspath(final_all)
                if files_exist([final_tc]):
                    bams[cond.id + "."+mapper]=os.path.abspath(final_tc)
    
    
    # write stats. FIXME: stats are not complete if pipeline was restarted with 
    print("Writing stats")
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
    
    # write slamstr config
    if 'slamstr_config_template' in config and 'mappers' in config:
        print("Writing slamstr config files")
        with open(config['slamstr_config_template'], 'r') as f:
            x = f.read().splitlines()
        for mapper in config['mappers'].keys():
            modi=['tc','ori'] if write_uncoverted else ['tc']
            for mode in modi:
                id = "slamstr."+mapper+"."+mode
                f_config = outdir + id+".json"
                with open(f_config, 'w') as out:
                    tokens={}
                    tokens["ID"]=id
                    tokens["genome_fa"]=os.path.abspath(config['genome_fa'])
                    tokens["gff"]=os.path.abspath(f_anno)
                    tokens["datasets"]="\n"
                    tokens["timepoints"]="\n"
                    for idx, cond in enumerate(conditions):
                        tokens["datasets"]+='\t"' + cond.id + '": "' + bams[cond.id + "."+mapper] + '"' + ("," if idx < len(conditions)-1 else "") + "\n"
                        tokens["timepoints"]+='\t"' + cond.id + '": "' + str(cond.timepoint)+ '"' + ("," if idx < len(conditions)-1 else "")+ "\n"
                    for l in x:
                        print(replace_tokens(l, tokens), file=out)
                        
    # create TDF files
    if 'create_tdf' in config and config['create_tdf']:
        print("Creating TDF files for %i bams" % (len(bams)))
        igvtools_log=tmpdir + config['dataset_name'] + ".igvtools.log" 
        for b in list(bams.values()):
            if args.force or not files_exist(b+".tdf"):
                create_tdf(b, b+".tdf", chrom_sizes, igvtools_cmd=igvtools_cmd, logfile=igvtools_log)
    
    print("All done in", datetime.timedelta(seconds=time.time()-startTime))