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

class ReadTruth:

    def __init__(self, name, chromosome, relStart, relEnd, relSeqError, absStart, absEnd, absSeqError, splicing):
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
        # Absolute starting position of the read alignment
        self.absStart = absStart
        # Absolute end position of the read alignment
        self.absEnd = absEnd
        # Absolute position of a sequencing error
        self.absSeqError = absSeqError
        # Splicing status
        self.splicing = splicing

    def __repr__(self):
        return "\t".join([self.name, str(self.chromosome), str(self.absStart), str(self.absEnd), str(self.absSeqError), str(self.splicing)])

class TruthCollection:

    def __init__(self, truthFile):
        self._truthFile = truthFile
        self._truth = None

    def readTruth(self):

        self._truth = pd.read_csv(self._truthFile, delimiter='\t')
        self._truth.index = self._truth.read_name

    def getTruth(self, name):

        row = self._truth.loc[name, :]

        if len(row) > 0:
            return ReadTruth(row.read_name, row.chr_abs, row.start_rel, row.end_rel, row.seq_err_pos_rel, row.start_abs, row.end_abs, row.seq_err_pos_abs, row.read_spliced)
        else:
            return None

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

    def evaluate(self, truth):

        print(self.chromosome,end="\t")
        print(truth.chromosome,end="\t")
        print(self.absStart,end="\t")
        print(truth.absStart,end="\t")
        print(self.absEnd,end="\t")
        print(truth.absEnd,end="\t")
        seqErrorMeasured = self.absSeqError
        seqErrorTruth = truth.absSeqError
        print(len(seqErrorMeasured),end="\t")
        if type(seqErrorTruth) is float:
            seqErrorTruth = list()
        else :
            seqErrorTruth = seqErrorTruth.split(",")
            seqErrorTruth = list(map(int, seqErrorTruth))
        positionalMismatches = 0

        for position in seqErrorMeasured:
            if position in seqErrorTruth:
                positionalMismatches += 1
        print(len(seqErrorTruth),end="\t")
        print(positionalMismatches,end="\t")

        print(self.splicing,end="\t")
        print(truth.splicing == 1)

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
                if self.isTCMismatch(readBase, refBase, read.is_reverse):
                    conversions.append(refPos + 1)
                else :
                    errors.append(refPos + 1)

        spliced = False
        for operation in read.cigartuples:
            if operation[0] == 3:
                spliced = True
                break

        simulatedRead = SimulatedRead(name, chromosome, 0, 0, list(), list(),
        start + 1, end, errors, conversions, spliced)

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
parser.add_argument("-t","--truth", type=existing_file, required=True, dest="truthFile", help="Read truth tsv file")
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

truthCollection = TruthCollection(args.truthFile)

truthCollection.readTruth()

print("Done reading truth collection")

simFile = SpliceSimFile(args.bamFile)

readIterator = simFile.readsGenomeWide()

print("\t".join(["chromsome_observed","chromosome_truth",
                 "start_observed","start_truth",
                 "end_observed","end_truth",
                 "mm_observed","mm_truth",
                 "matching_positions",
                 "splicing_observed","splicing_truth"]))

for read in readIterator:
    truth = truthCollection.getTruth(read.name)
    read.evaluate(truth)
    sys.stdin.readline()
