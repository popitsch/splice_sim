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

    def __init__(self, name, chromosome, relStart, relEnd, relSeqError, relConversion, absStart, absEnd, absSeqError, absConversion, splicing):
        # Name of the parsed read
        self.name = name
        # Chromosome
        self.chromosome = chromosome
        # Relative starting position of the read alignment
        self.relStart = relStart
        # Relative end position of the read alignment
        self.relEnd = relEnd
        # Relative position of a sequencing error
        self.relSeqError = relSeqError
        # Relative position of a nucleotide conversion
        self.relConversion = relConversion
        # Absolute starting position of the read alignment
        self.absStart = absStart
        # Absolute end position of the read alignment
        self.absEnd = absEnd
        # Absolute position of a sequencing error
        self.absSeqError = absSeqError
        # Absolute position of a nucleotide conversion
        self.absConversion = absConversion
        # Splicing status
        self.splicing = splicing

    def __repr__(self):
        return "\t".join([self.name, self.chromosome, str(self.absStart), str(self.absEnd), str(self.absSeqError), str(self.absConversion), str(self.splicing)])


class SimulatedReadIterator:

    def __init__(self, readIterator):
        self._readIterator = readIterator

    def __iter__(self):
        return self

    def __next__(self):

        read = self._readIterator.__next__()

        name = read.query_name

        chromosome = read.reference_name

        start = read.reference_start

        end = read.reference_end

        # Fix later to only skip splice site but include Indels
        matches = read.get_aligned_pairs(matches_only=True, with_seq=True)

        conversions = list()
        errors = list()

        for match in matches:
            readPos, refPos, refBase = match

            readBase = read.query_sequence[readPos]

            readBase = readBase.upper()
            refBase = refBase.upper()

            if readBase.upper() != refBase.upper():
                if not self.isTCMismatch(readBase, refBase, read.is_reverse):
                    conversions.append(refPos)
                else :
                    errors.append(refPos)

        spliced = False
        for operation in read.cigartuples:
            if operation[0] == 3:
                spliced = True
                break

        simulatedRead = SimulatedRead(name, chromosome, 0, 0, list(), list(),
        start, end, errors, conversions, spliced)

        return simulatedRead

    def isTCMismatch(self, readBase, refBase, isReverse):
        if(isReverse):
            return readBase == "G" and refBase == "A"
        else:
            return readBase == "C" and refBase == "T"

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
