import pysam
import os, sys

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

def parse_info(info):
    """ parse GFF3 info section """
    return {k:v for k,v in [a.split('=') for a in info.split(';') if '=' in a]}

def reverse_complement(seq):
    """ Calculate reverse complement DNA sequence """
    rev=[]
    for c in seq[::-1]:
        c=c.upper()
        rev+=[COMP_TABLE.get(c, 'N')]
    return ''.join(rev)

def pad_n(seq, minlen):
    """ Add up-/downstream padding with N's to ensure a given minimum length of the passed sequence """ 
    ret = seq
    if ( len(ret)<minlen):
        pad0="N" * int((minlen-len(ret))/2)
        pad1="N" * int(minlen-(len(ret)+len(pad0))) 
        ret=pad0+ret+pad1
    return (ret)

def bgzip_and_tabix(in_file, out_file=None, create_index=True, del_uncompressed=True, seq_col=0, start_col=1, end_col=1, line_skip=0, zerobased=False):
    """ Will BGZIP the passed file and creates a tabix index with the given params if create_index is True"""
    if out_file is None:
        out_file=in_file+'.gz'
    pysam.tabix_compress(in_file, out_file, force=True) # @UndefinedVariable
    if create_index:
        pysam.tabix_index(out_file, force=True, seq_col=seq_col, start_col=start_col, end_col=end_col, meta_char='#', line_skip=line_skip, zerobased=zerobased) # @UndefinedVariable
    if del_uncompressed:
        os.remove(in_file)
        
def localize_config(config):
    """ Recursively iterates json and adds /Volumes prefix to /group paths"""
    for k, v in config.items():
        if isinstance(v, str) and v.startswith('/groups'):
            config[k]='/Volumes'+ config[k]
        elif isinstance(v, dict):
            config[k]=localize_config(v)
        elif isinstance(v, list):
            config[k]=['/Volumes'+x if str(x).startswith('/groups') else x for x in v]
    return config
