#!/usr/bin/env python

from argparse import ArgumentParser, RawDescriptionHelpFormatter
import os, sys
import gzip
import pysam
import re
import pandas as pd
import numpy as np

#
# Extract TC mismatches over slamdunk reads and compare with simulated mismatches
#
usage = '''                           

  Copyright (C) 2021 XXX.  All rights reserved.

  Distributed on an "AS IS" basis without warranties
  or conditions of any kind, either express or implied.

USAGE
'''
if __name__ == '__main__':
    parser = ArgumentParser(description=usage, formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument("-i", "--bam", type=str, required=True, dest="bam_file",
                        metavar="bam_file", help="BAM file")
    parser.add_argument("-b", "--bed", type=str, required=True, dest="bed_file",
                        metavar="bed_file", help="BED file")
    parser.add_argument("-o", "--outfile", type=str, required=False, dest="outfile", metavar="outfile",
                        help="output file")
    args = parser.parse_args()

    opener = gzip.open if ".gz" in args.bed_file else open

    if args.outfile:
        fout = open(args.outfile, 'w')
    else :
        fout = sys.stdout

    bamfile = pysam.AlignmentFile(args.bam_file, "rb")

    with opener(args.bed_file, 'rt') as f:
        for line in f:
            if line.startswith("#") :
                continue
            chr, start, end, tx, score, strand, source, feature, phase, info = line.rstrip().split("\t")

            conversionDifferences = list()
            for read in bamfile.fetch(chr, int(start), int(end)):
                simulatedConversions = read.qname
                simulatedConversions = int(re.sub(".*_","",re.sub("_\\d$","",simulatedConversions)))
                detectedConversions = read.get_tag("TC")
                conversionDifferences.append(detectedConversions - simulatedConversions)

            mean = np.mean(np.array(conversionDifferences))
            readNumber = len(conversionDifferences)

            print("\t".join([chr, start, end, tx, str(mean), str(readNumber)]))