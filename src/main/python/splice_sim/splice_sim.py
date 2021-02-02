#
#
# Copyright (C) 2022 XXX.  All rights reserved.
#
# This file is part of splice_sim
# 
# See the file LICENSE for redistribution information.
#
# @author: niko.popitsch
"""
    Main file that does all commandline handling
"""
from argparse import ArgumentParser, RawDescriptionHelpFormatter
import os, sys, json, logging, random
from utils import existing_file
from collections import OrderedDict
import os, sys
from simulator import simulate_dataset
from evaluator import evaluate_dataset

VERSION = "0.1"
LOGO = """
Splice_sim v%s
""" % VERSION

#============================================================================
# splice sim readme
#
# NOTE: splice_sim will not simulate reads for annotations that are shorter than the readlength (as this is 
# not supported by ART.
#============================================================================

def check_config(config):
    for section in ['dataset_name','mappers','genome_fa','gene_gff','conditions','isoform_mode']:
        assert section in config, "Missing configuration section %s" % section
    if 'transcripts' not in config:
        assert "transcript_data" in config, "Transcript data needs to be configured either in config file ('transcripts' section) or in an external file referenced via 'transcript_data'"
        
usage = '''                           

  Copyright (C) 2021 XXX.  All rights reserved.

  Distributed on an "AS IS" basis without warranties
  or conditions of any kind, either express or implied.

USAGE
'''
if __name__ == '__main__':
    MODULES = ["simulate", "evaluate"]
    #============================================================================    
    if len(sys.argv) <= 1 or sys.argv[1] in ['-h', '--help']:
        print("usage: splice_sim.py [-h] " + ",".join(MODULES))
        sys.exit(1)
    mod = sys.argv[1]  
    if mod not in MODULES:
        print("Invalid module %s selected. Please use one of " + ",".join(MODULES))
        sys.exit(1)
        
    parser = {}
    parser["simulate"] = ArgumentParser(description=usage, formatter_class=RawDescriptionHelpFormatter)
    parser["simulate"].add_argument("-c", "--config", type=existing_file, required=True, dest="config_file", metavar="config_file", help="JSON config file")
    parser["simulate"].add_argument("-o", "--outdir", type=str, required=False, dest="outdir", metavar="outdir", help="output directory (default is current dir)")
    parser["simulate"].add_argument("--overwrite", required=False, action="store_true", default=False, dest="overwrite", help="If set, existing files will be overwritten")

    parser["evaluate"] = ArgumentParser(description=usage, formatter_class=RawDescriptionHelpFormatter)
    parser["evaluate"].add_argument("-c", "--config", type=existing_file, required=True, dest="config_file", metavar="config_file", help="JSON config file")
    parser["evaluate"].add_argument("-o", "--outdir", type=str, required=False, dest="outdir", metavar="outdir", help="output directory (default is current dir)")
    parser["evaluate"].add_argument("--overwrite", required=False, action="store_true", default=False, dest="overwrite", help="If set, existing files will be overwritten")
   
    print(LOGO)
    print("module: " + mod)
    args = parser[mod].parse_args(sys.argv[2:])
    #============================================================================
    config = json.load(open(args.config_file), object_pairs_hook=OrderedDict)
    check_config(config)
    cdir = os.path.dirname(os.path.abspath(args.config_file))+'/'
    print("config dir", cdir)
    
    # output dir (current dir if none provided)
    outdir = (args.outdir if args.outdir else os.getcwd()) + '/' + config['dataset_name'] +'/'
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    print("Logging to %s" % outdir+'splice_sim.log')
    logging.basicConfig(filename=outdir+'splice_sim.log', level=logging.DEBUG)    
    with open(outdir+'/splicing_simulator.effective_conf.json', 'w') as out:
        print(json.dumps(config, indent=4, sort_keys=True), file=out)
    if "random_seed" in config:
        random.seed(config["random_seed"])
    logging.info(LOGO)

    if mod == "simulate":
        simulate_dataset(config, cdir, outdir + 'sim/', args.overwrite)
    if mod == "evaluate":
        evaluate_dataset(config, cdir, outdir + 'sim/', outdir + 'eva/', args.overwrite)
