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
from splice_sim.simulator import create_genome_bam, postfilter_bam
from splice_sim.evaluator import *

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
    parser["build_model"].add_argument("-o", "--outdir", type=str, required=False, dest="outdir", metavar="outdir", help="output directory (default is current dir)")

    parser["create_genome_bam"] = ArgumentParser(description=usage, formatter_class=RawDescriptionHelpFormatter)
    parser["create_genome_bam"].add_argument("-m", "--model", type=str, required=True, dest="model_file", metavar="model_file", help="model file")
    parser["create_genome_bam"].add_argument("-a", "--art_sam", type=str, required=True, dest="art_sam_file", metavar="config_file", help="ART sam file")
    parser["create_genome_bam"].add_argument("-t", "--threads", type=int, required=False, default=1,  dest="threads", metavar="threads", help="threads")
    parser["create_genome_bam"].add_argument("-o", "--outdir", type=str, required=False, dest="outdir", metavar="outdir", help="output directory (default is current dir)")

    parser["postfilter_bam"] = ArgumentParser(description=usage, formatter_class=RawDescriptionHelpFormatter)
    parser["postfilter_bam"].add_argument("-c", "--config", type=str, required=True, dest="config_file", metavar="config_file", help="JSON config file")
    parser["postfilter_bam"].add_argument("-b", "--bam", type=str, required=True, dest="bam_file", metavar="bam_file", help="input bam file")
    parser["postfilter_bam"].add_argument("-o", "--outdir", type=str, required=True, dest="outdir", metavar="outdir", help="output dir")

    parser["evaluate_bam_performance"] = ArgumentParser(description=usage, formatter_class=RawDescriptionHelpFormatter)
    parser["evaluate_bam_performance"].add_argument("-c", "--config", type=str, required=True, dest="config_file", metavar="config_file", help="JSON config file")
    parser["evaluate_bam_performance"].add_argument("-m", "--model", type=str, required=True, dest="model_file", metavar="model_file", help="model file")
    parser["evaluate_bam_performance"].add_argument("-b", "--bam_file", type=str, required=True, dest="bam_file", metavar="bam_file", help="input bam file")
    parser["evaluate_bam_performance"].add_argument("-o", "--outdir", type=str, required=True, dest="outdir", metavar="outdir", help="output dir")

    parser["evaluate_splice_site_performance"] = ArgumentParser(description=usage, formatter_class=RawDescriptionHelpFormatter)
    parser["evaluate_splice_site_performance"].add_argument("-c", "--config", type=str, required=True, dest="config_file", metavar="config_file", help="JSON config file")
    parser["evaluate_splice_site_performance"].add_argument("-m", "--model", type=str, required=True, dest="model_file", metavar="model_file", help="model file")
    parser["evaluate_splice_site_performance"].add_argument("-b", "--bam_file", type=str, required=True, dest="bam_file", metavar="bam_file", help="input bam file")
    parser["evaluate_splice_site_performance"].add_argument("-t", "--truth_bam_dir", type=str, required=True, dest="truth_bam_dir", metavar="truth_bam_dir", help="truth bam dir")
    parser["evaluate_splice_site_performance"].add_argument("-o", "--outdir", type=str, required=True, dest="outdir", metavar="outdir", help="output dir")

    parser["calculate_splice_site_mappability"] = ArgumentParser(description=usage, formatter_class=RawDescriptionHelpFormatter)
    parser["calculate_splice_site_mappability"].add_argument("-c", "--config", type=str, required=True, dest="config_file", metavar="config_file", help="JSON config file")
    parser["calculate_splice_site_mappability"].add_argument("-m", "--model", type=str, required=True, dest="model_file", metavar="model_file", help="model file")
    parser["calculate_splice_site_mappability"].add_argument("-b", "--bam_file", type=str, required=True, dest="bam_file", metavar="bam_file", help="input bam file")
    parser["calculate_splice_site_mappability"].add_argument("-t", "--truth_bam_dir", type=str, required=True, dest="truth_bam_dir", metavar="truth_bam_dir", help="truth bam dir")
    parser["calculate_splice_site_mappability"].add_argument("-o", "--outdir", type=str, required=True, dest="outdir", metavar="outdir", help="output dir")

    parser["calc_feature_overlap"] = ArgumentParser(description=usage, formatter_class=RawDescriptionHelpFormatter)
    parser["calc_feature_overlap"].add_argument("-c", "--config", type=str, required=True, dest="config_file", metavar="config_file", help="JSON config file")
    parser["calc_feature_overlap"].add_argument("-m", "--model", type=str, required=True, dest="model_file", metavar="model_file", help="model file")
    parser["calc_feature_overlap"].add_argument("-b", "--bam_file", type=str, required=True, dest="bam_file", metavar="bam_file", help="input bam file")
    parser["calc_feature_overlap"].add_argument("-o", "--outdir", type=str, required=True, dest="outdir", metavar="outdir", help="output dir")

    parser["extract_splice_site_features"] = ArgumentParser(description=usage, formatter_class=RawDescriptionHelpFormatter)
    parser["extract_splice_site_features"].add_argument("-c", "--config", type=str, required=True, dest="config_file", metavar="config_file", help="JSON config file")
    parser["extract_splice_site_features"].add_argument("-m", "--model", type=str, required=True, dest="model_file", metavar="model_file", help="model file")
    parser["extract_splice_site_features"].add_argument("-o", "--outdir", type=str, required=True, dest="outdir", metavar="outdir", help="output dir")

    parser["extract_transcript_features"] = ArgumentParser(description=usage, formatter_class=RawDescriptionHelpFormatter)
    parser["extract_transcript_features"].add_argument("-c", "--config", type=str, required=True, dest="config_file", metavar="config_file", help="JSON config file")
    parser["extract_transcript_features"].add_argument("-m", "--model", type=str, required=True, dest="model_file", metavar="model_file", help="model file")
    parser["extract_transcript_features"].add_argument("-o", "--outdir", type=str, required=True, dest="outdir", metavar="outdir", help="output dir")

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
    outdir = os.path.abspath(args.outdir if args.outdir else os.getcwd())+'/'
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    # logging
    print("Logging to %s" % outdir+'splice_sim.log')
    logging.basicConfig(filename=outdir+'splice_sim.log', level=logging.DEBUG)

    logging.info(LOGO)


    if mod == "build_model":
        # load and check onfig
        config=json.load(open(args.config_file), object_pairs_hook=OrderedDict)
        check_config(config)
        Model.build_model(args.config_file, outdir)

    if mod == "create_genome_bam":
        # load model
        m = Model.load_from_file(args.model_file)
        config = m.config
        # set random seed
        if "random_seed" in config:
            random.seed(config["random_seed"])
            logging.info("setting random seed to %i" % config["random_seed"])
        create_genome_bam(m, args.art_sam_file, args.threads, outdir)

    if mod == "postfilter_bam":
        postfilter_bam(args.config_file, args.bam_file, args.outdir)

    if mod == "evaluate_bam_performance":
        # load config to be able to react to config changes after model was built!
        config=json.load(open(args.config_file), object_pairs_hook=OrderedDict)
        # load model
        m = Model.load_from_file(args.model_file)
        # set random seed
        if "random_seed" in config:
            random.seed(config["random_seed"])
            logging.info("setting random seed to %i" % config["random_seed"])
        evaluate_bam_performance(config, m, args.bam_file, outdir)

    if mod == "evaluate_splice_site_performance":
        # load config to be able to react to config changes after model was built!
        config=json.load(open(args.config_file), object_pairs_hook=OrderedDict)
        # load model
        m = Model.load_from_file(args.model_file)
        # set random seed
        if "random_seed" in config:
            random.seed(config["random_seed"])
            logging.info("setting random seed to %i" % config["random_seed"])
        evaluate_splice_site_performance(config, m, args.bam_file, args.truth_bam_dir, outdir)

    if mod == "calculate_splice_site_mappability":
        # load config to be able to react to config changes after model was built!
        config=json.load(open(args.config_file), object_pairs_hook=OrderedDict)
        # load model
        m = Model.load_from_file(args.model_file)
        # set random seed
        if "random_seed" in config:
            random.seed(config["random_seed"])
            logging.info("setting random seed to %i" % config["random_seed"])
        calculate_splice_site_mappability(config, m, args.bam_file, args.truth_bam_dir, outdir)

    if mod == "calc_feature_overlap":
        # load config to be able to react to config changes after model was built!
        config=json.load(open(args.config_file), object_pairs_hook=OrderedDict)
        # load model
        m = Model.load_from_file(args.model_file)
        # set random seed
        if "random_seed" in config:
            random.seed(config["random_seed"])
            logging.info("setting random seed to %i" % config["random_seed"])
        calc_feature_overlap(config, m, args.bam_file, outdir)

    if mod == "extract_splice_site_features":
        # load config to be able to react to config changes after model was built!
        config=json.load(open(args.config_file), object_pairs_hook=OrderedDict)
        # load model
        m = Model.load_from_file(args.model_file)
        # set random seed
        if "random_seed" in config:
            random.seed(config["random_seed"])
            logging.info("setting random seed to %i" % config["random_seed"])
        extract_splice_site_features(config, m, outdir)

    if mod == "extract_transcript_features":
        # load config to be able to react to config changes after model was built!
        config=json.load(open(args.config_file), object_pairs_hook=OrderedDict)
        # load model
        m = Model.load_from_file(args.model_file)
        # set random seed
        if "random_seed" in config:
            random.seed(config["random_seed"])
            logging.info("setting random seed to %i" % config["random_seed"])
        extract_transcript_features(config, m, outdir)
