'''
Utililty methods

@author: niko
'''


if __name__ == '__main__':
    pass

from argparse import ArgumentTypeError
import hashlib
import logging
import os
import random, string
import signal
from subprocess import *
import subprocess
import sys
from zipfile import ZipFile
import gzip
from shutil import copyfile

#============================================================================
# constants
#============================================================================
# @see https://pysam.readthedocs.io/en/latest/api.html#pysam.AlignedSegment.cigar
BAM_CMATCH=0
BAM_CINS=1
BAM_CDEL=2
BAM_CREF_SKIP=3
BAM_CSOFT_CLIP=4
BAM_CHARD_CLIP=5
BAM_CPAD=6
BAM_CEQUAL=7
BAM_CDIFF=8
BAM_CBACK=9

# Table of reverse complement bases
COMP_TABLE = {
    "A": 'T', "C": 'G', "T": 'A', "G": 'C'
    }
# Colors for read classifications
colors = {"exonexonfalse" : "228,26,28",
          "exonexontrue" : "55,126,184",
          "exonintron" : "77,175,74",
          "garbage" : "255,127,0",
          "intron" : "200,200,200"
}
def reverse_complement(seq):
    """ Calculate reverse complement DNA sequence """
    rev=[]
    for c in seq[::-1]:
        c=c.upper()
        rev+=[COMP_TABLE.get(c, 'N')]
    return ''.join(rev)

def to_region(dat):
    """ Create region string (chr:start-end) """
    return(dat.Chromosome + ":" + str(dat.Start) + "-" + str(dat.End))
def overlaps(a,b):
    """ x1 <= y2 && y1 <= x2 """
    return a[0]<=b[1] and a[1]>=b[0]
def calc_cov_a(a,b):
    """ coverage of interval a """
    cov = min(a[1],b[1]) - max(a[0],b[0]) + 1
    if cov<0:
        cov=0
    return cov
def calc_coverage(true_chr, read_chr, true_tuples, read_tuples):
    """ Calculate the coverage (number of bases) of a mapped read with the passed aligned blocks (r.get_blocks()) wrt. the passed true tuples. """
    if true_chr != read_chr:
        return 0
    # note: intervals are non-overlapping and there are usually very few (1-2) intervals per set
    cov=0
    for a in true_tuples:
        for b in read_tuples:
            #print(a,b,calc_cov_a(a,b))
            cov+=calc_cov_a(a,b)
    return cov
def pad_n(seq, minlen):
    """ Add up-/downstream padding with N's to ensure a given minimum length of the passed sequence """ 
    ret = seq
    if ( len(ret)<minlen):
        pad0="N" * int((minlen-len(ret))/2)
        pad1="N" * int(minlen-(len(ret)+len(pad0))) 
        ret=pad0+ret+pad1
    return (ret)
def add_tc(seq, strand, conversion_rate):
    """ introduces TC / AG conversions
        returns the new sequence and the converted positions (always! 5'-3' as shown in genome browser)
        positions are 0-based """
    ref="T"
    alt="C"
    conv=[]
    convseq=""
    for i,c in enumerate(seq):
        if c == ref and random.uniform(0, 1) < conversion_rate:
            c=alt.lower()
            conv+=[i if strand=="+" else len(seq)-i-1] # converted positions
        convseq+=c
    return convseq, conv  

def replace_tokens(s, tok):
    """ replace all tokens in a strings """
    for t in tok.keys():
        s=s.replace("@"+t+"@", tok[t])
    return(s)


def readidx2genpos_rel( read, idx ):
    """ Converts read positions to genomic positions (1-based) by taking all cigar operations into account. 
        The returned coords are relative to the alignment start, ie. absolute coords require to add read.reference_start """
    if len(read.get_aligned_pairs(matches_only=True))<idx:
        #print("invalid index %i, possibly softclipped? " % ( idx) )
        #print(read)
        return None
    return (read.get_aligned_pairs(matches_only=True)[idx-1][1]-read.reference_start+1)

def readidx2genpos_abs( read, idx ):
    """ Converts read positions to genomic positions (1-based) by taking all cigar operations into account. """
    pairs = read.get_aligned_pairs(matches_only=True)
    if idx<0 or len(pairs)<=idx or len(pairs[idx-1])==0:
        return None
    return (pairs[idx][1]+1)
 
def cigar_to_rel_pos(read):
    """ converts a (art_illumina) cigar string to a set of relative ref-mismatch (=seqerr) positions (0-based) """
    pos=[]
    off = 0
    for op, len in read.cigartuples:
        if (op == BAM_CMATCH) or (op == BAM_CEQUAL): # M or =
           off+=len
        elif op == BAM_CDIFF: # X
          for i in range(0, len):
              off+=1
              p=readidx2genpos_abs(read, off)
              if p is not None:
                  pos+=[p-1]    
    return pos

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

def parse_art_cigar(read):
    """ converts a (art_illumina) cigar string to a set of relative ref-mismatch (=seqerr) positions (0-based) """
    seqerr_pos=[]
    block_tuples=[]
    off = 0
    read_off=0
    start=0
    for op, l in read.cigartuples:
        if (op == BAM_CMATCH) or (op == BAM_CEQUAL): # M or =
           off+=l
           read_off+=l
        elif op == BAM_CDIFF: # X
            seqerr_pos+=[readidx2genpos_abs(read, read_off)]
            if seqerr_pos[-1] is None:
                logging.warn("Invalid seqerr position calculated for %s / %i. Ignoring!" % (read.query_name, read_off))
                del seqerr_pos[-1]
            off+=l
            read_off+=l
        elif op == BAM_CINS: # I
            # TODO seqerr reporting does not work for insertions
            #seqerr_pos+=[readidx2genpos_abs(read, off-1)]
            #if seqerr_pos[-1] is None:
            #    logging.warn("Invalid seqerr position calculated for %s / %i. Ignoring!" % (read.query_name, off))
            #    del seqerr_pos[-1]
            read_off+=l
        elif op == BAM_CDEL: # D
            block_tuples+=[(start, off-1)]
            # TODO  seqerrr reporting does not work for deletions
            #seqerr_pos+=[readidx2genpos_abs(read, off)]
            #if seqerr_pos[-1] is None:
            #    logging.warn("Invalid seqerr position calculated for %s / %i. Ignoring!" % (read.query_name, off))
            #    del seqerr_pos[-1]
            off+=l
            start=off
        else:
            print("UNSUPPORTED CIGAR operator %i" % op)
    block_tuples+=[(start, off-1)]   
    if read.is_reverse:
        block_tuples = [(y,x) for (x,y) in reversed(block_tuples)]
    return block_tuples, seqerr_pos


# FIXME: there is a problem with called subprocessed that do not terminate! This will hang this code!
# example-cmd:
# 'bowtie2', '-t', '10', '-x', '/project/ngs-work/meta/reference/genomes/hg19_human/bowtie2_path/hg19', '-1', '../FASTQ/Probe7_CAGATC_L004_R1_001.fastq.gz', '-2', '../FASTQ/Probe7_CAGATC_L004_R2_001.fastq.gz', '-S', '../FASTQ/Probe7_CAGATC_L004_R1_001.fastq.gz.bt2.sam'
# This will result in a non-terminating PERL job. IN turn, the python program also stalls. 
# run a commandline task
def runTask(cmd, shell=False):
    logging.info(cmd)
    if shell:
        out = check_output(" ".join(cmd), shell=True, stderr=subprocess.STDOUT)
    else: 
        out = check_output(cmd, stderr=subprocess.STDOUT)
    return out
        
def files_exist(files):
    if (type(files) is list) :
        for f in files:
            if f is None:
                return False
            if not os.path.exists(f):
                return False
    else:
        if files is None:
            return False
        if not os.path.exists(files):
            return False
    return True

# use in argument parser
def existing_file(files):
    if files_exist(files):
        return files
    if (type(files) is list) :
        raise ArgumentTypeError("Not all files exist ["+",".join(files)+"]")
    else:
        raise ArgumentTypeError("Not all files exist ["+(files)+"]")
    
# remove a (list of) file(s) (if it/they exists)
def removeFile(files):
    if (type(files) is list) :
        for f in files:
            if f is not None and os.path.exists(f):
                os.remove(f)
    else:
        if os.path.exists(files):
            os.remove(files)

# check whether a file exists and exit if not        
def checkFile(files):
    if (type(files) is list) :
        for f in files:
            checkFile(f)
    else :
        if not files is None and not os.path.exists(files):
            print("Error: file", files, "was not found! Exiting...")
            sys.exit(1)       
                
def pipelineStep(inputfile, outFile, cmd, shell=False, stdout=None, append=False, logfile=None):
    try:
        if inputfile is not None:
            if (type(inputfile) is list):
                for f in inputfile:
                    checkFile(f)
            else:
                checkFile(inputfile) 
        
        if stdout is not None:
            if shell is False:
                raise ArgumentTypeError("When using the parameter stdout, the shell parameter must be set to True!")
            else: 
                if append:
                    cmd.append(">>")
                else:
                    cmd.append(">")
                cmd.append(stdout)                                                                                                                                                                                                                                                                 

                
        # Overwrite output file?!
        if outFile is not None:
            if (type(outFile) is list) :
                for f in outFile:
                    if files_exist(f):
                        logging.warn("Overwriting file %s", f)
                else:
                    if files_exist(outFile):
                        logging.warn("Overwriting file %s", outFile)

        out = runTask(cmd, shell)         
        if logfile is None:
            logging.info(out)
        else:
            with open(logfile, "a") as log:
                log.write(out)
        
        if outFile is not None:
            checkFile(outFile)
        return True
    except CalledProcessError as e:
        logging.error(e.output)
        logging.error("ERROR %s - removing outputfile %s", e, outFile)
        if outFile is not None:
            if (type(outFile) is list) :
                for f in outFile:
                    removeFile([f]) 
            else:
                removeFile([outFile]) 
        return False        

# sample <number> of reads randomly from inFile (.fa/.fq) and writes to outFile. <number> can either be a absolute number (e.g. 1 000 000) or a percentage (e.g. 0.1 for 10 %)
def sampleFromFastX(inFile, outFile, number, seed=100):     
    return pipelineStep(inFile, outFile, ["seqtk", "sample", "-s", str(seed), inFile, str(number)], shell=True, stdout=outFile)


# Replaces the file extension of inFile to with <newExtension> and adds a suffix
# Example replaceExtension("reads.fq", ".sam", suffix="_ngm") => reads_ngm.sam
def replaceExtension(inFile, newExtension, suffix=""):    
    return os.path.splitext(inFile)[0] + suffix + newExtension        

def replaceExtension2(inFile, oldExtension, newExtension, suffix=""):    
    if inFile.endswith(oldExtension):
        inFile = inFile[:-len(oldExtension)]
    return inFile + suffix + newExtension        


# returns the filename without path 
def getFilename(f):    
    return os.path.basename(f)
# returns the filename without path and extension
def getFilenameNoExt(f):    
    return os.path.splitext(os.path.basename(f))[0]
        
# Converts a bam file to fq using bam2fastq
def bam2fq(inFile, outFile, delinFile=False, override=False):
    success = True
    
    if((not files_exist(outFile[::-1].replace("#"[::-1], "_1"[::-1], 1)[::-1]) and not files_exist(outFile[::-1].replace("#", "", 1)[::-1])) or override):
        success = pipelineStep(inFile, outFile, ["bam2fastq", "-f --aligned --unaligned --no-filtered", "-o", outFile, inFile], shell=True)
        if(success and delinFile):
            removeFile(inFile)
    else:
        logging.info("bam2fq skipped. File " + outFile[::-1].replace("#"[::-1], "_1"[::-1], 1)[::-1] + " already exists.")
    return success

def runFASTQC(inFile, outFolder, threads=1, additionalParameters="", override=False):
    success = True
    if (not files_exist(outFolder)):
        os.makedirs(outFolder)
    outFile = os.path.join(outFolder, os.path.basename(inFile).split(".")[0] + "_fastqc.zip")
    if(not files_exist(outFile) or override):
        if (type(inFile) is list) :
            success = pipelineStep(inFile, None, ["fastqc", "-t", str(threads), "-o", outFolder, additionalParameters] + inFile, shell=True)
        else:
            success = pipelineStep(inFile, None, ["fastqc", "-t", str(threads), "-o", outFolder, additionalParameters, inFile], shell=True)
        
    return success

def indexBam(inFileBam, override=False):
    success = True
    idxFile = inFileBam + ".bai"
    if(not files_exist(idxFile) or override):
        success = success and pipelineStep(inFileBam, idxFile, ["samtools", "index", inFileBam], shell=True)
    return success

def tabix(inFileGz, override=False, additionalParameters=[]):
    success = True
    idxFile = inFileGz + ".tbi"
    if(not files_exist(idxFile) or override):
        success = success and pipelineStep(inFileGz, idxFile, ["tabix"] + additionalParameters +[inFileGz], shell=True)
    return success

def indexVcf(inFileVcf, override=False):
    success = True
    idxFile = inFileVcf + ".tbi"
    if(not files_exist(idxFile) or override):
        success = success and pipelineStep(inFileVcf, idxFile, ["tabix", "-p", "vcf", inFileVcf], shell=True)
    return success
    
def sam2bam(inFile, outFile, index=True, sort=True, delinFile=False, override=False, onlyUnique=False, onlyProperPaired=False, onlyPrimary=False, filterMQ=0, maxmem=None):
    if(onlyUnique and filterMQ == 0):
        filterMQ = 1;
        
    success = True    
    if(not files_exist(outFile) or override):        
        cmd = ["samtools", "view", "-Sb", inFile, "-o", outFile]
        if filterMQ > 0:
            cmd += ["-q", str(filterMQ)]
        if onlyProperPaired:
            cmd += ["-f", "2"]
        if onlyPrimary:
            cmd += ["-F", "256"]
        success = success and pipelineStep(inFile, outFile, cmd, shell=True)
        
        if(sort):         
            tmp = outFile + "_tmp"
            os.rename(outFile, tmp)   
            cmd = ["samtools", "sort"]
            if not maxmem is None:
                cmd += ["-m", maxmem ]
            cmd+=["-o", outFile]    
            cmd += [tmp]             
            success = success and pipelineStep(tmp, outFile, cmd, shell=True)
            if(success):
                removeFile(tmp)
        if(index):
            success = success and indexBam(outFile)
            
        if(success and delinFile):
            removeFile(inFile)
    else:
        logging.info("Skipping sam2bam. " + outFile + " already exists.");
    return success

def sambamba2bam(inFile, outFile, index=True, sort=True, delinFile=False, override=False, ncpu=1, maxmem=None, sortByName=False, EXE="sambamba"):

    success = True    
    if(not files_exist(outFile) or override):        
        cmd = [EXE, "view"]
        if ( inFile.endswith("sam")):
            cmd+=["-S"] 
        cmd+=["-f", "bam", "-t", str(ncpu), inFile, "-o", outFile]
        success = success and pipelineStep(inFile, outFile, cmd, shell=True)
        if(sort):    
            sortedBam = replaceExtension(outFile, ".sorted.bam")  
            cmd = [EXE, "sort", "-t", str(ncpu)]
            if sortByName:
                cmd+= ["-n"]
            if not maxmem is None:
                cmd += ["-m", maxmem ]    
            cmd += [outFile, sortedBam]             
            success = success and pipelineStep(outFile, sortedBam, cmd, shell=True)
            if(success):
                os.rename(sortedBam, outFile)   
                removeFile(replaceExtension(outFile, ".sorted.bam.bai"))
        if(index):
            success = success and indexBam(outFile)
            
        if(success and delinFile):
            removeFile(inFile)
    else:
        logging.info("Skipping sambamba2bam. " + outFile + " already exists.");
    return success

def sambambamerge(outFile, inFiles=[], ncpu=1, override=False,delinFile=True, EXE="sambamba"):
    success = True    
    if(not files_exist(outFile) or override):        
        cmd = [EXE, "merge", "-t", str(ncpu), outFile]
        for bam in inFiles:
            cmd+=[bam]
        success = success and pipelineStep(inFiles, outFile, cmd, shell=True)
        if(success and delinFile):
            removeFile(inFiles)
            for f in inFiles:
                removeFile(f+".bai")
    else:
        logging.info("Skipping sambambamerge. " + outFile + " already exists.");
    return success

def sortbam(inFile, outFile, index=False, override=False, delinFile=False, maxmem=None):
    success = True    
    if(not files_exist(outFile) or override):     
        cmd = ["samtools", "sort", "-o", outFile]
        if not maxmem is None:
            cmd += ["-m", maxmem ]
        cmd += [ inFile ]             
        success = success and pipelineStep(inFile, outFile, cmd, shell=True)
        if(index):
            success = success and indexBam(outFile)          
        if(success and delinFile):
            removeFile(inFile)
    else:
        logging.info("Skipping sortbam. " + outFile + " already exists.");
    return success

def sambambasortbam(inFile, outFile, index=False, override=False, delinFile=False, maxmem=None, ncpu=1, EXE="sambamba"):
    success = True    
    if(not files_exist(outFile) or override):     
        cmd = [EXE, "sort", "-t", str(ncpu)]
        cmd += ["-o", outFile]
        if not maxmem is None:
            cmd += ["-m", maxmem ]
        cmd += [ inFile]             
        success = success and pipelineStep(inFile, outFile, cmd, shell=True)
        if(index):
            success = success and indexBam(outFile)          
        if(success and delinFile):
            removeFile(inFile)
    else:
        logging.info("Skipping sortbam. " + outFile + " already exists.");
    return success

def gzipFile(inFile):
    success = True    
    success = success and pipelineStep(inFile, inFile + ".gz", ["gzip", inFile], shell=True)
    return success

def bgzip(inFile, outFile=None, index=False, override=False, delinFile=False, threads=None):
    if outFile == None:
        outFile = inFile+".gz"
    success = True    
    if(not files_exist(outFile) or override):      
        cmd=["bgzip", "-c", inFile]
        if threads:
            cmd+=["--threads", str(threads)]
        success = success and pipelineStep(inFile, outFile, cmd, shell=True, stdout=outFile)
        if(index):
            if outFile.endswith(".vcf"):
                success = success and indexVcf(outFile)
            else:
                success = success and tabix(outFile)
        if(success and delinFile):
            removeFile(inFile)
    else:
        logging.info("Skipping bgzip. " + outFile + " already exists.");
    return success

def flagstat(bam):
    success = True
    flagstat = bam + ".flagstat"
    if not files_exist(flagstat):
        success = success and pipelineStep(None, flagstat, [ "samtools", "flagstat", bam], shell=True, stdout=flagstat)
    return success

def bamstats(bam):
    success = True
    bamstat = bam + ".bamstat"
    if not files_exist(bamstat):
        success = success and pipelineStep(None, bamstat, [ "bamstats.py", "-o", bamstat, bam], shell=True)
    return success

def runPrePy(inFile, ref, outFile=None, sort=False, passOnly=False, additionalParameters=[]):
    success = True   
    if sort:
        sortedF=os.path.splitext(inFile)[0]+".sorted.vcf.gz"
        sortVcf(inFile, sortedF)
        inFile = sortedF
    if outFile == None:
        outFile = os.path.splitext(inFile)[0]+".norm.vcf.gz"
    cmd = ["pre.py", "-r", ref]
    if additionalParameters:
        cmd+=additionalParameters
    if passOnly:
        cmd+=["--pass-only"]   
    cmd+=[ inFile, outFile ] 
    success = success and pipelineStep(inFile, outFile, cmd, shell=True )


def sortVcf(inFile, outFile, override=False):
    if outFile == None:
        outFile = os.path.splitext(inFile)[0]+".sorted.vcf.gz"
    success = True    
    if(not files_exist(outFile) or override):     
        presort = outFile + ".pre"
        success = success and pipelineStep(inFile, presort, ["vcf-sort", "-c", inFile], shell=True, stdout=presort)
        sortedF = outFile + ".sorted"
        success = success and pipelineStep(presort, sortedF, ["ngs-tools", "VCFTools", "sortCanonical", "-i", presort, "-o", sortedF], shell=True)
        if outFile.endswith(".gz"):
            bgzip(sortedF, outFile, index=True, delinFile=True)
        else:
            os.rename(sortedF, outFile)       
        if success:
            removeFile(presort)
    else:
        logging.info("Skipping sortVCF. " + outFile + " already exists.");
    return success


def sortBed(inFile, outFile, override=False, fnc=False, can=False, delinFile=False):
    if outFile == None:
        outFile = os.path.splitext(inFile)[0]+".sorted.bed.gz"
    success = True    
    if(not files_exist(outFile) or override):     
        presort = outFile + ".pre"
        success = success and pipelineStep(inFile, presort, ["sort", "-k", "1,1", "-k2,2n", inFile], shell=True, stdout=presort)
        sortedF = outFile + ".sorted"
        cmd = ["ngs-tools", "BEDTools", "sort", "-i", presort, "-o", sortedF]
        if fnc:
            cmd+=["-fnc"]
        if can:
            cmd+=["-can"]
        success = success and pipelineStep(presort, sortedF, cmd, shell=True)
        if outFile.endswith(".gz"):
            bgzip(sortedF, outFile, index=True, delinFile=True)
        else:
            os.rename(sortedF, outFile)       
        if success:
            removeFile(presort)
            if delinFile:
                removeFile(inFile)
    else:
        logging.info("Skipping sortBed. " + outFile + " already exists.");
    return success


def multiqc(dirs, outDir):
    success = True
    cmd = ["multiqc", "-o", outDir]
    if (type(dirs) is list) :
        cmd+=dirs
    else:
        cmd+=[dirs]
    success = success and pipelineStep(None, os.path.join(outDir, "multiqc_report.html"),cmd, shell=True)
    return success

def sambambaflagstat(bam, ncpu=1, EXE="sambamba"):
    success = True
    flagstat = bam + ".flagstat"
    if not files_exist(flagstat):
        success = success and pipelineStep(None, flagstat, [ EXE, "flagstat", "-t", str(ncpu), bam], shell=True, stdout=flagstat)
    return success

def mapqHist(bam):
    success = True
    mapqhist = bam + ".mapqhist"
    if not files_exist(mapqhist):
        success = success and pipelineStep(bam, mapqhist, ["ngs-tools", "QualityDistribution", "printMAPQHistogram", "-r", bam, "-o", mapqhist ])
    return success

def rtg_vcfstats(vcf):
    success = True
    vcfstat = vcf + ".rtgstat"
    if not files_exist(vcfstat):
        success = success and pipelineStep(None, vcfstat, [ "rtg", "vcfstats", vcf], shell=True, stdout=vcfstat)
    return success

def hashsum(string):
    ret=""
    if (type(string) is list) :
        for s in string:
            ret += hashsum(s)
    else:
        m = hashlib.md5()
        m.update(string)
        ret+=m.hexdigest()
    return ret

def randstr(length):
    return ''.join(random.choice(string.ascii_lowercase) for i in range(length))

def touch(fname, times=None):
    with open(fname, 'a'):
        os.utime(fname, times)
        
        
def runParallel(cmdFile, jobs=4):
    success = True
    cmd = ["parallel", "--eta", "-j", str(jobs), "--load", "80%", "--noswap", "-a", cmdFile]
    success = success and pipelineStep(cmdFile, None, cmd, shell=True)
    return success

def addOrReplaceReadGroups(inbam, outbam, RGID="0", RGLB="lib1", RGPL="illumina", RGPU="unit1", RGSM="20", delinFile=False):
    success = True
    cmd = ["picard-tools", "AddOrReplaceReadGroups", 
           "I="+inbam, 
           "O="+outbam, 
           "RGID="+RGID,
           "RGLB="+RGLB,
           "RGPL=" + RGPL,
           "RGPU=" + RGPU,
           "RGSM="+ RGSM]
    if(success and delinFile):
            removeFile(inbam)
    success = success and pipelineStep(inbam, outbam, cmd, shell=True)
    return success


def filterBAM(bam, out, ncpu=1, minMapQ=20):
    success = True
    cmd = ["samtools", "view", "-@", str(ncpu), "-q", str(minMapQ), "-F", "4", "-b", bam]
    success = success and pipelineStep(bam, out, cmd, shell=True, stdout=out)
    success = success and indexBam(out)
    success = success and sambambaflagstat(out, ncpu=ncpu)
    success = success and bamstats( out )
    return success



def wig2bw(wig, out, chrSize):
    success = True
    cmd = ["wigToBigWig", wig, chrSize, out, "-clip"]
    success = success and pipelineStep([wig, chrSize], out, cmd, shell=True)
    return success

#
# if strandBy is set, additionalParameters is ignored.
#
def bam2tdf(inbam, outtdf, chrsizes, scale=1.0, strandBy=0, additionalParameters=[]):
    if ( strandBy == 1): # STRANDS_BY_READ
        additionalParameters=["-countFlags", "1"]
    if ( strandBy == 2): # STRANDS_BY_FIRST_IN_PAIR
        additionalParameters=["-countFlags", "2"]
    if ( strandBy == 4): # STRANDS_BY_SECOND_IN_PAIR
        additionalParameters=["-countFlags", "4"]
    success = True
    cmd = ["genomic-tools", "bamtools", "extractScaledTDF"]
    cmd +=["-i", inbam] 
    cmd +=["-o", outtdf] 
    cmd +=["-c", chrsizes] 
    cmd +=["-s", str(scale)] 
    cmd +=["-i", inbam] 
    cmd += additionalParameters
    success = success and pipelineStep(inbam, outtdf, cmd, shell=True)
    return success


#
# if strandBy is set, additionalParameters is ignored.
#
def bam2bed(inbam, outbed, additionalParameters=[]):
    success = True
    cmd = ["bedtools", "bamtobed"]
    cmd +=["-i", inbam] 
    cmd += additionalParameters
    tmp = outbed+".tmp"
    success = success and pipelineStep(inbam, tmp, cmd, shell=True, stdout=tmp)    
    sortBed(tmp, outbed, delinFile=True)
    return success


# @see https://stackoverflow.com/questions/15008758/parsing-boolean-values-with-argparse
def str2bool(v):
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise ArgumentTypeError('Boolean value expected.')
    
# runs wget on a (list of) URI(s). @return the created file(s)    
def wget(uri, redo=True, additionalParameters=[]):
    success = True
    ret=[]
    if (type(uri) is list) :
        for u in uri:
            ret+=wget(u, additionalParameters)
    else:
        fn=uri[uri.rfind('/')+1:]
        if not files_exist( fn ) or redo:
            cmd = ["wget", "-q", uri] + additionalParameters
            success = success and pipelineStep(None, fn, cmd, shell=True)
        if success: 
            ret+=[fn]       
    return ret

def fileNumLines(fname):
    if ( fname.endswith(".gz")):
        with gzip.open(fname, 'rb') as f:
            for i, l in enumerate(f):
                pass
        return i + 1
    else:
        with open(fname) as f:
            for i, l in enumerate(f):
                pass
        return i + 1

def numFastqReads(fname):
    fl = fileNumLines(fname)
    return fl / 4.0


def createIgvGenome(outdir, id, fasta, bed=None, alias=None):
    if not outdir.endswith("/"):
        outdir += "/"
    if not os.path.exists(outdir):
            print("Creating dir " + outdir)
            os.makedirs(outdir)
    
    # local bed file?
    localbed = None
    if bed is not None:
        localbed = outdir + id + ".genes.bed"
        copyfile(bed, localbed)        
    # local alias file?
    localalias = None
    if alias is not None:
        localalias = outdir + id + ".alias"
        copyfile(alias, localalias)    
       
    propfile = outdir + "property.txt"
    with open(propfile, 'w') as f:
        f.write("fasta=true\n")
        f.write("fastaDirectory=false\n")
        f.write("ordered=true\n")
        f.write("id="+id+"\n")
        f.write("name="+id+"\n")
        if localbed is not None:
            f.write("geneFile="+localbed+"\n")
        if localalias is not None:
            f.write("chrAliasFile="+localalias+"\n")
        f.write("sequenceLocation="+fasta+"\n")
        
    genomefile = outdir + id + ".genome"
    with ZipFile(genomefile, 'w') as f:
        f.write(propfile, os.path.basename(propfile))
        removeFile(propfile)
        if localbed is not None:
            f.write(localbed, os.path.basename(localbed))
            removeFile(localbed)
        if localalias is not None:
            f.write(localalias, os.path.basename(localalias))
            removeFile(localalias)

def partitionDict(dictionary, chunks=2):

    partition = [dict() for idx in range(chunks)]
    idx = 0
    for k, v in dictionary.items():
        partition[idx][k] = v
        if idx < chunks - 1:  # indexes start at 0
            idx += 1
        else:
            idx = 0
    return partition