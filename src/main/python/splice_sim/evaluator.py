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


def evaluate_bam(bam_file, category, mapper, condition, out):
    samin = pysam.AlignmentFile(bam_file, "rb")
    for read in samin.fetch():
        read_name=read.query_name
        true_tid,true_strand,true_isoform,tag,true_chr,true_read_cigar,true_seqerr,tc_pos = read_name.split("_")
        read_chr = 'NA' if read.is_unmapped else read.reference_name
        read_coord='NA' if read.is_unmapped else '%s:%i-%i' % (read.reference_name, read.reference_start, read.reference_end )
        true_tuples=[tuple([int(y) for y in x.split('-')]) for x in true_read_cigar.split(',')]
        read_tuples=read.get_blocks()
        overlap=calc_coverage(true_chr, read_chr, true_tuples, read_tuples)
        print("%s\t%s\t%s\t%s\t%s\t%i" % (read_coord, read_name,category,mapper,condition,overlap), file=out)
        

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
            tfile = config_dir+"/"+tfile
        assert files_exist(tfile), "transcript_data file not found: %s" % tdata
        logging.info("Reading transcript configuration from external config file %s" % tfile)
        tdata = json.load(open(tfile), object_pairs_hook=OrderedDict)
        config["transcripts"]=tdata
        
    # instantiate model
    m = Model( config )
    logging.info("Configured conditions: %s" % (", ".join(str(x) for x in m.conditions)) )

    fout = outdir+'reads.tsv'
    with open(fout, 'w') as out:
        print("coords\tread_name\tcategory\tmapper\tcondition\toverlap", file=out)
        for cond in m.conditions:
            for mapper in config['mappers'].keys():
                bamdir_all = simdir + "bam_ori/" + mapper + "/"
                bamdir_tc  = simdir + "bam_tc/" + mapper + "/"
                final_all  = bamdir_all + config['dataset_name'] + "." + cond.id + "."+mapper+".bam"
                final_tc   = bamdir_tc  + config['dataset_name'] + "." + cond.id + "."+mapper+".TC.bam"
                evaluate_bam(final_all,'all',mapper,cond.id, out)    
                evaluate_bam(final_tc, 'tc', mapper,cond.id, out)    

    logging.info("All done in %s" % str(datetime.timedelta(seconds=time.time()-startTime)))
    print("All done in", datetime.timedelta(seconds=time.time()-startTime))