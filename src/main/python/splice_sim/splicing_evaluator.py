#!/usr/bin/env python3

'''
@author: tobias.neumann@imp.ac.at
'''


from argparse import ArgumentParser, RawDescriptionHelpFormatter, \
    ArgumentTypeError
from collections import *
import csv, datetime, time, logging, sys, os, json
import pysam
import pyranges as pr
import pandas as pd
import numpy as np
from Bio.Seq import Seq
from utils import *
import random
import math
import gzip

#============================================================================
# splice eval readme
#
#============================================================================
VERSION = "0.1"
logo = '''
===================================
Slamstr splice_eval v%s
===================================
''' % VERSION

#============================================================================
usage = '''

  Copyright 2020 Tobias Neumann. All rights reserved.

  Licensed under the Apache License 2.0
  http://www.apache.org/licenses/LICENSE-2.0

  Distributed on an "AS IS" basis without warranties
  or conditions of any kind, either express or implied.

USAGE
'''


class SimulatedRead:

    def __init__(self):
        # Name of the parsed read
        self.name = None
        # Relative starting position of the read alignment
        self.relStart = None
        # Relative end position of the read alignment
        self.relEnd = None
        # Relative position of a sequencing error
        self.relSeqError = []
        # Absolute starting position of the read alignment
        self.absStart = None
        # Absolute end position of the read alignment
        self.absEnd = None
        # Absolute position of a sequencing error
        self.absSeqError = []
        # Splicing status
        self.splicing = None

    def __repr__(self):
        return "\t".join([self.name, str(self.relStart), str(self.relEnd), self.relSeqError, str(self.absStart), str(self.absEnd), self.absSeqError, self.splicing])


class SimulatedReadIterator:

    def __init__(self, readIterator):
        self._readIterator = readIterator

    def __iter__(self):
        return self

    def __next__(self):

        read = self._readIterator.__next__()

        print(read)

        simulatedRead = SimulatedRead()

        return simulatedRead

class SpliceSimFile:

    def __init__(self, bamFile):
        self._bamFile = pysam.AlignmentFile(bamFile, "rb")

    def readsGenomeWide(self):

        return SimulatedReadIterator(self._bamFile.fetch())

    def atoi(self, text):
        return int(text) if text.isdigit() else text

    def natural_keys(self, text):
        '''
        alist.sort(key=natural_keys) sorts in human order
        http://nedbatchelder.com/blog/200712/human_sorting.html
        (See Toothy's implementation in the comments)
        '''
        return [ self.atoi(c) for c in re.split('(\d+)', text) ]

parser = ArgumentParser(description=usage, formatter_class=RawDescriptionHelpFormatter)
parser.add_argument("-b","--bam", type=existing_file, required=True, dest="bamFile", help="Mapped bam file")
parser.add_argument("-o","--out", type=str, required=False, default="outdir", dest="outdir", metavar="outdir", help="Output folder")
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

simFile = SpliceSimFile(args.bamFile)

readIterator = simFile.readsGenomeWide()

for read in readIterator:
    print(read)
    sys.stdin.readline()
