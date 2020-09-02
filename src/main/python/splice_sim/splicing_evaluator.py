#!/usr/bin/env python3

'''
@author: tobias.neumann@imp.ac.at
'''


from argparse import ArgumentParser, RawDescriptionHelpFormatter, \
    ArgumentTypeError
from collections import *
from intervaltree import Interval, IntervalTree
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
import re

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

# Functions

# Adapted to GFF3 from
## https://github.com/DaehwanKimLab/hisat2/blob/master/hisat2_extract_splice_sites.py


def extract_splice_sites(gff_file, verbose=False):
    genes = defaultdict(list)
    trans = {}

    fh = gzip.open(gff_file, 'rt')

    # Parse valid exon lines from the GFF file into a dict by transcript_id
    for line in fh:

        line = line.strip()
        if not line or line.startswith('#'):
            continue
        if '#' in line:
            line = line.split('#')[0].strip()

        try:
            chrom, source, feature, left, right, score, \
                strand, frame, values = line.split('\t')
        except ValueError:
            continue
        left, right = int(left), int(right)

        if feature != 'exon' or left >= right:
            continue

        values_dict = {}
        for attr in values.split(';'):
            if attr:
                attr, _, val = attr.strip().partition('=')
                values_dict[attr] = val.strip('"')

        if 'gene_id' not in values_dict or \
                'transcript_id' not in values_dict:
            continue

        transcript_id = values_dict['transcript_id']
        if transcript_id not in trans:
            trans[transcript_id] = [chrom, strand, [[left, right]]]
            genes[values_dict['gene_id']].append(transcript_id)
        else:
            trans[transcript_id][2].append([left, right])

    # Sort exons and merge where separating introns are <=5 bps
    for tran, [chrom, strand, exons] in trans.items():
            exons.sort()
            tmp_exons = [exons[0]]
            for i in range(1, len(exons)):
                if exons[i][0] - tmp_exons[-1][1] <= 5:
                    tmp_exons[-1][1] = exons[i][1]
                else:
                    tmp_exons.append(exons[i])
            trans[tran] = [chrom, strand, tmp_exons]

    # Calculate and print the unique junctions
    junctions = dict()
    for transcript, value in trans.items():
        chrom, strand, exons = value

        if not chrom in junctions:
            junctions[chrom] = IntervalTree()
        for i in range(1, len(exons)):
            junctions[chrom][exons[i-1][1]:exons[i][0]] = transcript

    return junctions

# Classes

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

class ReadEvaluation:

    def __init__(self, name, trueChromosome, measuredChromosome, trueStart, measuredStart, trueEnd, measuredEnd, trueSeqErrors, measuredSeqErrors, errorMatches, trueSplicing, measuredSplicing, knownSpliceSites, novelSpliceSites):
        # Name of the read
        self.name = name
        # True chromosome
        self.trueChromosome = trueChromosome
        # True chromosome
        self.measuredChromosome = measuredChromosome
        # True chromosome
        self.trueStart = trueStart
        # True chromosome
        self.measuredStart = measuredStart
        # True chromosome
        self.trueEnd = trueEnd
        # True chromosome
        self.measuredEnd = measuredEnd
        # True chromosome
        self.trueSeqErrors = trueSeqErrors
        # True chromosome
        self.measuredSeqErrors = measuredSeqErrors
        # True chromosome
        self.errorMatches = errorMatches
        # True chromosome
        self.trueSplicing = trueSplicing
        # True chromosome
        self.measuredSplicing = measuredSplicing
        # True chromosome
        self.knownSpliceSites = knownSpliceSites
        # True chromosome
        self.novelSpliceSites = novelSpliceSites

    def __repr__(self):
        return "\t".join([self.name, str(self.measuredChromosome), str(self.trueChromosome), str(self.measuredStart), str(self.trueStart), str(self.measuredEnd),
        str(self.trueEnd), str(self.measuredSeqErrors), str(self.trueSeqErrors),str(self.errorMatches),str(self.measuredSplicing),str(self.trueSplicing),
        str(self.knownSpliceSites), str(self.novelSpliceSites)])

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

    def __init__(self, name, chromosome, relStart, relEnd, relSeqError, relConversion, absStart, absEnd, absSeqError, absConversion, splicing, spliceSites):
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
        # Splicing status
        self.spliceSites = spliceSites

    def hasSpliceSite(self, junction):
        if self.splicing:
            for spliceSite in self.spliceSites:
                res = junction.envelop(spliceSite[0], spliceSite[1])
                if (res) :
                    iv = res.pop()
                    begin, end, data = iv
                    if begin == spliceSite[0] and end == spliceSite[1]:
                        return True
        return False

    def evaluate(self, truth, junctions):

        seqErrorMeasured = self.absSeqError
        seqErrorTruth = truth.absSeqError

        if type(seqErrorTruth) is float:
            seqErrorTruth = list()
        else :
            seqErrorTruth = seqErrorTruth.split(",")
            seqErrorTruth = list(map(int, seqErrorTruth))
        positionalMismatches = 0

        for position in seqErrorMeasured:
            if position in seqErrorTruth:
                positionalMismatches += 1

        knownSpliceSite = 0
        novelSpliceSite = 0

        if self.splicing:
            for spliceSite in self.spliceSites:
                if self.chromosome in junctions:
                    res = junctions[self.chromosome].envelop(spliceSite[0], spliceSite[1])
                    if (res) :
                        iv = res.pop()
                        begin, end, data = iv
                        if begin == spliceSite[0] and end == spliceSite[1] and data == re.sub("_.*","",self.name):
                            knownSpliceSite += 1
                        else :
                            novelSpliceSite += 1
                    else :
                        novelSpliceSite += 1
                else :
                    novelSpliceSite += 1

        return ReadEvaluation(self.name, truth.chromosome, self.chromosome, truth.absStart, self.absStart, truth.absEnd, self.absEnd, len(seqErrorTruth), len(seqErrorMeasured), positionalMismatches, truth.splicing == 1, self.splicing, knownSpliceSite, novelSpliceSite)

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

        spliceSites = list()
        spliced = False

        refPositions = read.get_aligned_pairs(matches_only=False, with_seq=False)
        offset = 0
        for operation in read.cigartuples:
            if operation[0] == 3:
                spliced = True
                spliceSite = (refPositions[offset - 1][1] + 1,refPositions[offset + operation[1]][1] + 1)
                spliceSites.append(spliceSite)
            offset += operation[1]

        simulatedRead = SimulatedRead(name, chromosome, 0, 0, list(), list(),
        start + 1, end, errors, conversions, spliced, spliceSites)

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

    def readsInInterval(self, chromosome, start, end):

        return SimulatedReadIterator(self._bamFile.fetch(contig = chromosome, start = start, stop = end))

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
parser.add_argument("-g","--gff", type=existing_file, required=True, dest="gffFile", help="GFF3 file of gene models")
parser.add_argument("-t","--truth", type=existing_file, required=True, dest="truthFile", help="Read truth tsv file")
parser.add_argument('-i', "--intron-based", action='store_true', dest="introns", help="Quantify per intron, not per read")
parser.add_argument("-o","--out", type=str, required=False, default="outdir", dest="outdir", metavar="outdir", help="Output folder")
args = parser.parse_args()
startTime = time.time()

print(logo, file=sys.stderr)

# output dir (current dir if none provided)
outdir = args.outdir if args.outdir else os.getcwd()
if not outdir.endswith("/"):
    outdir += "/"
if not os.path.exists(outdir):
    print("Creating dir " + outdir, file=sys.stderr)
    os.makedirs(outdir)
tmpdir = outdir + "/tmp/"
if not os.path.exists(tmpdir):
    os.makedirs(tmpdir)

junctions = extract_splice_sites(args.gffFile)

truthCollection = TruthCollection(args.truthFile)

truthCollection.readTruth()

print("Done reading truth collection", file = sys.stderr)

simFile = SpliceSimFile(args.bamFile)

if args.introns:

    print("\t".join(["Transcript", "exon-intron", "intron", "intron-exon", "exon-exon", "exon-exon-false", "garbage"]))

    iv = junctions["6"][122707991:122710920].pop()
    ivTree = IntervalTree()
    ivTree.add(iv)
    readIterator = simFile.readsInInterval("6", iv.begin, iv.end)

    exonintron = 0
    exonexon = 0
    intronexon = 0
    intron = 0
    exonexonfalse = 0
    garbage = 0

    for read in readIterator:
        truth = truthCollection.getTruth(read.name)
        evaluation = read.evaluate(truth, junctions)
        start = read.absStart
        end = read.absEnd

        if start <= iv.begin and end > iv.begin:
            if read.splicing and not read.hasSpliceSite(ivTree):
                exonexonfalse += 1
            elif read.splicing and read.hasSpliceSite(ivTree):
                exonexon += 1
            elif end < iv.end:
                exonintron += 1
            else :
                garbage += 1

        elif start < iv.end and end >= iv.end :
            intronexon += 1
        elif start >= iv.begin and end <= iv.end:
            intron +=1
        else:
            garbage += 1

    print("\t".join([iv.data, str(exonintron), str(intron), str(intronexon), str(exonexon), str(exonexonfalse), str(garbage)]))


else:

    readIterator = simFile.readsGenomeWide()

    print("\t".join(["name","chromsome_observed","chromosome_truth",
                     "start_observed","start_truth",
                     "end_observed","end_truth",
                     "mm_observed","mm_truth",
                     "matching_positions",
                     "splicing_observed","splicing_truth",
                     "known_splice_sites","novel_splice_sites"]))

    for read in readIterator:
        truth = truthCollection.getTruth(read.name)
        evaluation = read.evaluate(truth, junctions)
        print(evaluation)
