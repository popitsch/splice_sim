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



     
def evaluate_dataset(conf_file, config_dir, overwrite=False):
    """ Evaluates a dataset. """
    startTime = time.time()

    # output dirs
    outdir = outdir + '/' + config['dataset_name'] + '/eva/'
    if not os.path.exists(outdir):
        print("Creating dir " + outdir)
        os.makedirs(outdir)
    tmpdir = outdir + "/tmp/"
    if not os.path.exists(tmpdir):
        os.makedirs(tmpdir)

    print("Logging to %s" % outdir+'/splice_sim.log')
    logging.basicConfig(filename=outdir+'/splicing_simulator.log', level=logging.DEBUG)
    logging.info("Evaluating dataset %s" % config['dataset_name'])

    
    logging.info("All done in %s" % str(datetime.timedelta(seconds=time.time()-startTime)))
    print("All done in", datetime.timedelta(seconds=time.time()-startTime))