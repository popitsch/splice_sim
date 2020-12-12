#!/usr/bin/env python3

'''
@author: tobias.neumann@imp.ac.at
'''

from model import Model, Condition, Isoform, Transcript

from argparse import ArgumentParser, RawDescriptionHelpFormatter, \
    ArgumentTypeError
from collections import *
from intervaltree import Interval, IntervalTree
from joblib import Parallel, delayed
import csv, datetime, time, logging, sys, os, json
import pysam
import pyranges as pr
import pandas as pd
from multiprocessing import Process, Manager
import numpy as np
from Bio.Seq import Seq
from utils import *
import random
import math
import gzip
import re
import shutil
from itertools import zip_longest

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

colors = {"exonexonfalse" : "228,26,28",
          "exonexontrue" : "55,126,184",
          "exonintron" : "77,175,74",
          "garbage" : "255,127,0",
          "intron" : "200,200,200"
}

# TODO: Maybe remove
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

#============================================================================
#    Classes
#============================================================================

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

    def isTruthPositionMatch(self) :
        return self.measuredChromosome == str(self.trueChromosome) and self.measuredStart == self.trueStart and self.measuredEnd == self.trueEnd and self.measuredSplicing == self.trueSplicing

    def __repr__(self):
        return "\t".join([self.name, str(self.measuredChromosome), str(self.trueChromosome), str(self.measuredStart), str(self.trueStart), str(self.measuredEnd),
        str(self.trueEnd), str(self.measuredSeqErrors), str(self.trueSeqErrors),str(self.errorMatches),str(self.measuredSplicing),str(self.trueSplicing),
        str(self.knownSpliceSites), str(self.novelSpliceSites)])

class TruthCollection:

    def __init__(self, truthFile):
        self._truthFile = truthFile
        self._truth = None

    def readTruth(self):

        self._truth = pd.read_csv(self._truthFile, delimiter='\t', dtype={'read_name': "string", 'start_rel' : 'int64', 'end_rel' : 'int64', 'seq_err_pos_rel' : "string", "chr_abs" : "string", "start_abs" : "int64", "end_abs" : "int64", "read_spliced" : 'bool', "seq_err_pos_abs" : "string"})
        self._truth.index = self._truth.read_name

    def getTruth(self, name):

        row = self._truth.loc[name, :]

        if len(row) > 0:
            return ReadTruth(row.read_name, row.chr_abs, row.start_rel, row.end_rel, row.seq_err_pos_rel, row.start_abs, row.end_abs, row.seq_err_pos_abs, row.read_spliced)
        else:
            return None

    def getTruthByInterval(self, chromosome, start, end) :

        subset = self._truth[self._truth.chr_abs.eq(chromosome) &
        (
        ((self._truth.start_abs < start) & (self._truth.end_abs >= start)) |
        ((self._truth.start_abs <= end) & (self._truth.end_abs > end)) |
        ((self._truth.start_abs >= start) & (self._truth.end_abs <= end))
        )
        ]

        #print(start)
        #print(end)
        #print(subset)
        #print(subset["read_name"].tolist())
        #print("ENSMUST00000112580.7_+_pre_4-1258" in subset["read_name"].tolist())
        #print("ENSMUST00000112580.7_+_pre_4-1258" in self._truth["read_name"].tolist())
        #print(self._truth[self._truth.read_name == "ENSMUST00000112580.7_+_pre_4-1258"])
        #print(self._truth[(self._truth.read_name == "ENSMUST00000112580.7_+_pre_4-1258") & ((self._truth.end_abs >= end))])
        #sys.stdin.readline()

        return subset["read_name"].tolist()

class SimulatedRead:

    def __init__(self, name, chromosome, relStart, relEnd, relSeqError, relConversion, absStart, absEnd, absSeqError, absConversion, splicing, spliceSites, bamRead):
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
        # Pysam read object
        self.bamRead = bamRead

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

        if pd.isna(seqErrorTruth):
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
        return "\t".join([self.name, self.chromosome, str(self.absStart), str(self.absEnd), str(self.absSeqError), str(self.absConversion), str(self.splicing), str(self.spliceSites)])


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

        softclippedlist = read.get_reference_positions(full_length=True)

        offsetStart = 0
        offsetEnd = 0

        while (softclippedlist[offsetStart] != start) :
            offsetStart += 1

        softclippedlist.reverse()

        # reference_end points to one past the last aligned residue
        while (softclippedlist[offsetEnd] != (end - 1)) :
            offsetEnd += 1

        start -= offsetStart
        end += offsetEnd

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
        start + 1, end, errors, conversions, spliced, spliceSites, read)

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

def partitionDict(dictionary, chunks=2):

    partition = [dict() for idx in range(chunks)]
    idx = 0
    for k, v in dictionary.items():
        partition[idx][k] = v
        if idx < chunks - 1:  # indexes start at 0
            idx += 1
        else:
            idx = 0
    return partition

def classifyIntron(transcripts, truthCollection, bamFile, outdir, thread):

    simFile = SpliceSimFile(bamFile)

    results = dict()

    headerfile = pysam.AlignmentFile(bamFile, "rb")
    classifiedReads = pysam.AlignmentFile(outdir + "/" + str(thread) + ".bam", "wb", template=headerfile)

    for tid in transcripts:

        results[tid] = list()

        t = transcripts[tid]

        for i in range(0, len(t.introns)):

            txid = t.introns.iloc[i]['transcript_id']
            exonnumber = t.introns.iloc[i]['exon_number']
            chromosome = t.introns.iloc[i]['Chromosome']
            start = t.introns.iloc[i]['Start']
            end = t.introns.iloc[i]['End']

            truthReadSet = truthCollection.getTruthByInterval(chromosome, start, end)

            iv = Interval(start - 1, end + 1)
            ivTree = IntervalTree()
            ivTree.add(iv)

            readIterator = simFile.readsInInterval(chromosome, start - 1, end + 1)

            exonintron = 0
            exonexon = 0
            intronexon = 0
            intron = 0
            exonexonfalse = 0
            garbage = 0

            intronDict = {
                "tpspliced": dict(),
                "fpspliced": dict(),
                "fnspliced": dict(),
                "tpexonintron": dict(),
                "fpexonintron": dict(),
                "fnexonintron": dict(),
                "tpintron": dict(),
                "fpintron": dict(),
                "fnintron": dict()
            }

            for read in readIterator:

                truth = truthCollection.getTruth(read.name)
                evaluation = read.evaluate(truth, {chromosome: ivTree})

                start = read.absStart
                end = read.absEnd

                bamRead = read.bamRead

                color = ""

                if start <= iv.begin and end > iv.begin:

                    if read.splicing and not read.hasSpliceSite(ivTree):
                        intronDict["fpspliced"][read.name] = bamRead
                        # exonexonfalse += 1
                        # color = colors["exonexonfalse"]
                    elif read.splicing and read.hasSpliceSite(ivTree):
                        if read.name in truthReadSet:
                            if evaluation.isTruthPositionMatch():
                                intronDict["tpspliced"][read.name] = bamRead
                            else:
                                intronDict["fpspliced"][read.name] = bamRead
                        else:
                            intronDict["fpspliced"][read.name] = bamRead
                        # exonexon += 1
                        # color = colors["exonexontrue"]
                    elif end < iv.end:
                        if read.name in truthReadSet:

                            if evaluation.isTruthPositionMatch():
                                intronDict["tpexonintron"][read.name] = bamRead
                            else:
                                intronDict["fpexonintron"][read.name] = bamRead
                        else:
                            intronDict["fpexonintron"][read.name] = bamRead
                        # exonintron += 1
                        # color = colors["exonintron"]
                    else:
                        # garbage += 1
                        color = colors["intron"]

                elif start < iv.end and end >= iv.end:
                    if read.splicing and not read.hasSpliceSite(ivTree):
                        intronDict["fpspliced"][read.name] = bamRead
                        # exonexonfalse += 1
                        # color = colors["exonexonfalse"]
                    elif read.splicing and read.hasSpliceSite(ivTree):
                        if read.name in truthReadSet:
                            if evaluation.isTruthPositionMatch():
                                intronDict["tpspliced"][read.name] = bamRead
                            else:
                                intronDict["fpspliced"][read.name] = bamRead
                        else:
                            intronDict["fpspliced"][read.name] = bamRead
                        # exonexon += 1
                        # color = colors["exonexontrue"]
                    else:
                        if read.name in truthReadSet:
                            if evaluation.isTruthPositionMatch():
                                intronDict["tpexonintron"][read.name] = bamRead
                            else:
                                intronDict["fpexonintron"][read.name] = bamRead
                        else:
                            intronDict["fpexonintron"][read.name] = bamRead
                        # intronexon += 1
                        # color = colors["exonintron"]
                elif start >= iv.begin and end <= iv.end:
                    if read.splicing and not read.hasSpliceSite(ivTree):
                        intronDict["fpspliced"][read.name] = bamRead
                        # exonexonfalse += 1
                        # color = colors["exonexonfalse"]
                    elif read.splicing and read.hasSpliceSite(ivTree):
                        if read.name in truthReadSet:
                            if evaluation.isTruthPositionMatch():
                                intronDict["tpspliced"][read.name] = bamRead
                            else:
                                intronDict["fpspliced"][read.name] = bamRead
                        else:
                            intronDict["fpspliced"][read.name] = bamRead
                        # exonexon += 1
                        # color = colors["exonexontrue"]
                    else:
                        if read.name in truthReadSet:
                            if evaluation.isTruthPositionMatch():
                                intronDict["tpintron"][read.name] = bamRead
                            else:
                                intronDict["fpintron"][read.name] = bamRead
                        else:
                            intronDict["fpintron"][read.name] = bamRead
                        # intron +=1
                        # color = colors["intron"]
                else:
                    pass
                    # These are reads that start at base 1 of the exon
                    # garbage += 1
                    # color = colors["intron"]

                # bamRead.setTag('YC', color)
            for truthRead in truthReadSet:
                if not truthRead in intronDict["tpspliced"] and not truthRead in intronDict[
                    "fpspliced"] and not truthRead in intronDict["fnspliced"] and not truthRead in intronDict[
                    "tpexonintron"] and not truthRead in intronDict["fpexonintron"] and not truthRead in intronDict[
                    "fnexonintron"] and not truthRead in intronDict["tpintron"] and not truthRead in intronDict[
                    "fpintron"] and not truthRead in intronDict["fnintron"]:

                    readTruth = truthCollection.getTruth(truthRead)

                    if readTruth.splicing:
                        intronDict["fnspliced"][truthRead] = 0
                    elif readTruth.absStart >= (iv.begin + 1) and readTruth.absEnd <= (iv.end - 2):
                        intronDict["fnintron"][truthRead] = 0
                    else:
                        intronDict["fnexonintron"][truthRead] = 0

            for read in intronDict["tpspliced"]:
                bamRead = intronDict["tpspliced"][read]
                bamRead.setTag('YC', colors["exonexontrue"])
                classifiedReads.write(bamRead)

            for read in intronDict["fpspliced"]:
                bamRead = intronDict["fpspliced"][read]
                bamRead.setTag('YC', colors["exonexonfalse"])
                classifiedReads.write(bamRead)

            for read in intronDict["tpexonintron"]:
                bamRead = intronDict["tpexonintron"][read]
                bamRead.setTag('YC', colors["exonintron"])
                classifiedReads.write(bamRead)

            for read in intronDict["fpexonintron"]:
                bamRead = intronDict["fpexonintron"][read]
                bamRead.setTag('YC', colors["exonexonfalse"])
                classifiedReads.write(bamRead)

            for read in intronDict["tpintron"]:
                bamRead = intronDict["tpintron"][read]
                bamRead.setTag('YC', colors["intron"])
                classifiedReads.write(bamRead)

            for read in intronDict["fpintron"]:
                bamRead = intronDict["fpintron"][read]
                bamRead.setTag('YC', colors["exonexonfalse"])
                classifiedReads.write(bamRead)

            results[tid].append("\t".join(
                [txid + "_" + str(exonnumber), str(len(intronDict["tpspliced"])), str(len(intronDict["fpspliced"])),
                 str(len(intronDict["fnspliced"])), str(len(intronDict["tpexonintron"])),
                 str(len(intronDict["fpexonintron"])), str(len(intronDict["fnexonintron"])),
                 str(len(intronDict["tpintron"])), str(len(intronDict["fpintron"])), str(len(intronDict["fnintron"]))]))

    classifiedReads.close()

    return (results)

if __name__ == '__main__':

    parser = ArgumentParser(description=usage, formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument("-b","--bam", type=existing_file, required=True, dest="bamFile", help="Mapped bam file")
    parser.add_argument("-c","--config", type=existing_file, required=True, dest="confF", help="JSON config file")
    parser.add_argument("-t","--truth", type=existing_file, required=True, dest="truthFile", help="Read truth tsv file")
    parser.add_argument("-p","--parallel", type=int, required=False, default=3, dest="threads", help="Number of threads")
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

    # load + check config
    config = json.load(open(args.confF), object_pairs_hook=OrderedDict)

    truthCollection = TruthCollection(args.truthFile)

    truthCollection.readTruth()

    print("Done reading truth collection", file = sys.stderr)

    if args.introns:

        intronEvalFile = open(os.path.join(outdir, "eval_intron.txt"), "w")

        # TODO: Duplicated code from splicing_simulator.py

        # read transcript data from external file if not in config
        if 'transcripts' not in config:
            assert "transcript_data" in config, "Transcript data needs to be configured either in config file ('transcripts' section) or in an external file referenced via 'transcript_data'"
            tfile = config['transcript_data']
            if not os.path.isabs(tfile):
                tfile = os.path.dirname(os.path.abspath(args.confF))+"/"+tfile
            tdata = json.load(open(tfile), object_pairs_hook=OrderedDict)
            config["transcripts"]=tdata

        model = Model(config)

        print("Evaluating introns", file = sys.stderr)

        chunks = partitionDict(model.transcripts, args.threads)

        print("\t".join(["Transcript", "tp-exon-exon", "fp-exon-exon", "fn-exon-exon", "tp-exon-intron", "fp-exon-intron", "fn-exon-intron", "tp-intron", "fp-intron", "fn-intron"]), file = intronEvalFile)

        results = Parallel(backend = "multiprocessing", n_jobs=args.threads, verbose = 30)(delayed(classifyIntron)(chunks[chunk], truthCollection, args.bamFile, os.path.join(outdir, "tmp"), chunk) for chunk in range(len(chunks)))

        mergedResults = dict()

        for result in results:
            mergedResults.update(result)

        mergeBams = list()
        for chunk in range(len(chunks)):
            mergeBams.append(os.path.join(outdir, "tmp", str(chunk) + ".bam"))

        pysam.merge("-f", os.path.join(outdir, "tmp.bam"), *mergeBams)

        pysam.sort("-o", os.path.join(outdir, "classifiedreads.bam"), os.path.join(outdir, "tmp.bam"))

        pysam.index(os.path.join(outdir, "classifiedreads.bam"))

        os.remove(os.path.join(outdir, "tmp.bam"))

        shutil.rmtree(tmpdir)

        for tid in model.transcripts:
            if tid in mergedResults:
                for quantification in mergedResults[tid]:
                    print(quantification, file = intronEvalFile)

        intronEvalFile.close()

    else :

        junctions = extract_splice_sites(config["gene_gff"])

        simFile = SpliceSimFile(args.bamFile)

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
