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
from collections import OrderedDict
from splice_sim.model import Model
from splice_sim.simulator import create_genome_bam

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
    for section in ['dataset_name','mappers','genome_fa','gene_gff','isoform_mode']:
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
         
    parser = {}
    
    parser["build_model"] = ArgumentParser(description=usage, formatter_class=RawDescriptionHelpFormatter)
    parser["build_model"].add_argument("-c", "--config", type=str, required=True, dest="config_file", metavar="config_file", help="JSON config file")
    parser["build_model"].add_argument("-o", "--outdir", type=str, required=False, dest="out_dir", metavar="out_dir", help="output directory (default is current dir)")

    parser["create_genome_bam"] = ArgumentParser(description=usage, formatter_class=RawDescriptionHelpFormatter)
    parser["create_genome_bam"].add_argument("-m", "--model", type=str, required=True, dest="model_file", metavar="config_file", help="model file")
    parser["create_genome_bam"].add_argument("-a", "--art_sam", type=str, required=True, dest="art_sam_file", metavar="config_file", help="ART sam file")
    parser["create_genome_bam"].add_argument("-t", "--threads", type=int, required=False, default=1,  dest="threads", metavar="threads", help="threads")
    parser["create_genome_bam"].add_argument("-o", "--outdir", type=str, required=False, dest="out_dir", metavar="out_dir", help="output directory (default is current dir)")

    
    #============================================================================    
    if len(sys.argv) <= 1 or sys.argv[1] in ['-h', '--help']:
        print("usage: splice_sim.py [-h] " + ",".join(parser.keys()))
        sys.exit(1)
    mod = sys.argv[1]  
    if mod not in parser.keys():
        print("Invalid module '%s' selected. Please use one of %s" % (mod, ",".join(parser.keys())))
        sys.exit(1)

    args = parser[mod].parse_args(sys.argv[2:])
    #============================================================================

    # output dir (current dir if none provided)
    out_dir = os.path.abspath(args.out_dir if args.out_dir else os.getcwd())+'/'
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
        
    # logging    
    print("Logging to %s" % out_dir+'splice_sim.log')
    logging.basicConfig(filename=out_dir+'splice_sim.log', level=logging.DEBUG)    
    
    logging.info(LOGO)


    if mod == "build_model":
        # load and check onfig
        config=json.load(open(args.config_file), object_pairs_hook=OrderedDict)
        check_config(config)   
        Model.build_model(args.config_file, out_dir)

    if mod == "create_genome_bam":
        # load model
        m = Model.load_from_file(args.model_file)
        config = m.config
        # set random seed
        if "random_seed" in config:
            random.seed(config["random_seed"])
            logging.info("setting random seed to %i" % config["random_seed"])
        create_genome_bam(m, args.art_sam_file, args.threads, out_dir)