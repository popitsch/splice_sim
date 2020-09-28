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

# TODO: Maybe remove
# Functions

# Adapted to GFF3 from
## https://github.com/DaehwanKimLab/hisat2/blob/master/hisat2_extract_splice_sites.py

colors = {"exonexonfalse" : "228,26,28",
          "exonexontrue" : "55,126,184",
          "exonintron" : "77,175,74",
          "garbage" : "255,127,0",
          "intron" : "200,200,200"
}

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

# TODO: Put back
# Methods from splicing_simulator.py

#============================================================================
# utility methods
#============================================================================
def to_region(dat):
    return(dat.Chromosome + ":" + str(dat.Start) + "-" + str(dat.End))

# TODO: Put back
# Classes from splicing_simulator.py

#============================================================================
#    Classes
#============================================================================

class Stat():
    def __init__(self, id, value, cond=None, mapper=None ):
        self.id=id
        self.cond=cond
        self.mapper=mapper
        self.value=value

    def get_header(self):
        ret = ("id\tcond\tmapper\tvalue")
        return (ret)

    def __str__(self):
        ret = ("%s\t%s\t%s\t%s" % (self.id if self.id is not None else "NA",
                                   self.cond if self.cond is not None else "NA",
                                   self.mapper if self.mapper is not None else "NA",
                                   str(self.value)) )
        return (ret)

class Condition():
    def __init__(self, id, timepoint, conversion_rate, coverage ):
        self.id = id
        self.timepoint=timepoint
        self.conversion_rate=conversion_rate
        self.coverage = coverage
        self.bam="?"

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        ret = ("%s: [tp=%i, cr=%f, cv=%f]" % (self.id, self.timepoint, self.conversion_rate, self.coverage ) )
        return (ret)

class Isoform():
    def __init__(self, id, fractions, splicing_status, transcript ):
        self.id = id
        self.fractions = fractions
        self.splicing_status = splicing_status
        self.t = transcript
        # extract alignment blocks
        self.aln_blocks=[]
        bstart=self.t.transcript.Start
        for idx, splicing in list(enumerate(self.splicing_status)):
            if splicing==1:
                bend=self.t.introns.iloc[idx].Start-1 # last exonic 1-based pos
                self.aln_blocks+=[(bstart, bend)]
                bstart=self.t.introns.iloc[idx].End+1 # 1st exonic 1-based pos
        bend=self.t.transcript.End
        self.aln_blocks+=[(bstart, bend)]

    # convert isoform-relative coordinates to genomic coordinates.
    # E.g., a rel_pos==0 will return the 1st position of the transcript in genomic coordinates (1-based)
    def rel2abs_pos(self, rel_pos):
        abs_pos=None
        off = rel_pos
        block_id=0
        for bstart,bend in self.aln_blocks:
            bwidth = bend-bstart+1
            abs_pos = bstart
            if off - bwidth < 0:
                return (abs_pos+off, block_id)
            # next block
            off -= bwidth
            block_id+=1
        # TODO: check WARN: Invalid relative position 3641 / 142903114 in isoform ENSMUST00000100497.10_pre, [0.5], [0,0,0,0,0]
        print("WARN: Invalid relative position %i / %i in isoform %s" % (rel_pos, abs_pos, self))
        return (abs_pos+off, block_id)

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        ret = ("%s_%s, [%s], [%s]" % (self.t.tid, self.id, ",".join(str(x) for x in self.fractions), ",".join(str(x) for x in self.splicing_status)) )
        return (ret)

class Transcript():
    def __init__(self, tid, df, genome, conditions ):
        self.tid = tid
        self.is_valid = True
        self.df = df
        self.transcript = self.df[self.df.Feature=='transcript'].iloc[0]
        self.region = to_region(self.transcript)
        self.transcript_seq = genome.fetch(region=self.region)
        self.exons = self.df[self.df.Feature=='exon']
        self.introns = pd.DataFrame(columns=self.df.columns, index=[0])
        last_end = -1
        idata = []
        for index, ex in self.exons.iterrows():
            if last_end > 0: # skip 1st
                dat={ 'transcript_id': self.transcript.transcript_id,
                      'Feature': 'intron',
                      'Chromosome':self.transcript.Chromosome,
                      'Strand':self.transcript.Strand,
                      'Start':last_end+1,
                      'End':ex.Start-1,
                      'exon_number': int(ex.exon_number)
                      }
                intron = pd.DataFrame(data=dat, columns=self.df.columns, index=[0])
                idata.append(intron)
            last_end = ex.End
        if idata:
            self.introns=pd.concat(idata)
        else:
            self.introns=pd.DataFrame()
        # set abundance and isoforms
        self.cond = conditions
        self.abundance = math.ceil(config['transcripts'][tid]['abundance']) # abundance is float
        self.isoforms={}
        total_frac = [0] * len(self.cond)
        for iso in config['transcripts'][tid]['isoforms'].keys():
            frac=config['transcripts'][tid]['isoforms'][iso]['fractions']
            splicing_status=config['transcripts'][tid]['isoforms'][iso]['splicing_status'] if 'splicing_status' in config['transcripts'][tid]['isoforms'][iso] else []
            if self.transcript.Strand == "-":
                splicing_status = list(reversed(splicing_status)) # so 1st entry in config refers to 1st intron.
            # assert len(splicing_status)==len(self.introns), "misconfigured splicing status for some isoform/conditions: %s (%i vs %i)" % ( self.tid, len(splicing_status), len(self.introns))
            if len(splicing_status)!=len(self.introns):
                print("Misconfigured splicing status for some isoform/conditions of transcript %s (%i vs %i). Skipping transcript..." % ( self.tid, len(splicing_status), len(self.introns)))
                self.is_valid = False
                return
            self.isoforms[iso] = Isoform(iso, frac, splicing_status, self)
            total_frac = [x + y for x, y in zip(total_frac, frac)]
        for x in total_frac:
            # assert with rounding error
            assert x <= 1.00001, "Total fraction of transcript > 1 (%s) for some isoform/conditions: %s" % ( ",".join(str(x) for x in total_frac), self.tid )
        # add mature form if not configured
        if 'mat' not in self.isoforms:
            frac = [(1.0 - x) for x in total_frac]
            iso='mat'
            splicing_status=[1] * len(self.introns)
            self.isoforms[iso] = Isoform(iso, frac, splicing_status, self)
        #print('added mature isoform for %s: %s' % (tid, self.isoforms[iso]) )

    def get_dna_seq(self, splicing_status):
        seq = self.transcript_seq
        for idx, splicing in reversed(list(enumerate(splicing_status))):
            if splicing==1:
                pos1=self.introns.iloc[idx].Start - self.transcript.Start
                pos2=self.introns.iloc[idx].End - self.transcript.Start + 1
                seq = seq[:pos1]+seq[pos2:]
        return (seq)

    def get_rna_seq(self, iso):
        seq = self.get_dna_seq(self.isoforms[iso].splicing_status)
        if self.transcript.Strand == '-':
            seq = str(Seq(seq).reverse_complement())
        return (seq)

    def get_sequences(self):
        ret=OrderedDict()
        for cond_idx, cond in enumerate(self.cond):
            ret[cond]=OrderedDict()
            total_count=0
            for id in self.isoforms.keys():
                frac = self.isoforms[id].fractions[cond_idx]
                count=int(self.abundance * frac)
                total_count+=count
                ret[cond][id]=[
                    self.get_rna_seq(id),
                    count]
            toadd = self.abundance - total_count
            #assert toadd==0 | toadd==1 , "rounding error"
            if toadd > 0: # add 1 to last isoform
                print("corrected counts by %i" % (toadd))
                id = list(self.isoforms.keys())[-1]
                ret[cond][id]=[ret[cond][id][0], ret[cond][id][1] + toadd]
        return (ret)

    def to_bed(self):
        t=self.transcript
        return ("%s\t%i\t%i\t%s\t%i\t%s" % (t.Chromosome, t.Start, t.End, t.ID, 1000, t.Strand) )

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        ret = ("%s" % (self.transcript.transcript_id) )
        return (ret)


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

parser = ArgumentParser(description=usage, formatter_class=RawDescriptionHelpFormatter)
parser.add_argument("-b","--bam", type=existing_file, required=True, dest="bamFile", help="Mapped bam file")
parser.add_argument("-c","--config", type=existing_file, required=True, dest="confF", help="JSON config file")
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

# load + check config
config = json.load(open(args.confF), object_pairs_hook=OrderedDict)

if args.introns:

    # TODO: Duplicated code from splicing_simulator.py

    # read transcript data from external file if not in config
    if 'transcripts' not in config:
        assert "transcript_data" in config, "Transcript data needs to be configured either in config file ('transcripts' section) or in an external file referenced via 'transcript_data'"
        tfile = config['transcript_data']
        if not os.path.isabs(tfile):
            tfile = os.path.dirname(os.path.abspath(args.confF))+"/"+tfile
        tdata = json.load(open(tfile), object_pairs_hook=OrderedDict)
        config["transcripts"]=tdata

    # get conditions and configuration
    conditions=[]
    for id in config["conditions"].keys():
        conditions+=[Condition(id, config["conditions"][id][0], config["conditions"][id][1], config["conditions"][id][2] )]
    print("Configured conditions: %s" % (", ".join(str(x) for x in conditions)) , file = sys.stderr)

    # load + filter gene gff
    print("Loading gene GFF", file = sys.stderr)
    gff = pr.read_gff3(config["gene_gff"])
    gff_df = gff.df
    df = gff_df[gff_df['transcript_id'].isin(list(config['transcripts'].keys()))] # reduce to contain only configured transcripts

    # load genome
    print("Loading genome", file = sys.stderr)
    genome = pysam.FastaFile(config["genome_fa"])
    chrom_sizes = config["chrom_sizes"] if 'chrom_sizes' in config else config["genome_fa"]+".chrom.sizes"

    # instantiate transcripts
    transcripts=OrderedDict()
    for tid in list(config['transcripts'].keys()):
        tid_df = df[df['transcript_id']==tid] # extract data for this transcript only
        if tid_df.empty:
            print("No annotation data found for configured tid %s, skipping..." % (tid), file = sys.stderr)
        else:
            t = Transcript(tid, tid_df, genome, conditions)
            if t.is_valid:
                transcripts[tid] = t

    print("Evaluating introns", file = sys.stderr)

    headerfile = pysam.AlignmentFile(args.bamFile, "rb")
    classifiedReads = pysam.AlignmentFile(outdir + "/classifiedreads.bam", "wb", template=headerfile)

    simFile = SpliceSimFile(args.bamFile)

    # TODO: Wrap in proper function

    print("\t".join(["Transcript", "exon-intron", "intron", "intron-exon", "exon-exon", "exon-exon-false"]))

    for tid in transcripts:
        t = transcripts[tid]
        for i in range(0, len(t.introns)):

            txid = t.introns.iloc[i]['transcript_id']
            exonnumber = t.introns.iloc[i]['exon_number']
            chromosome = t.introns.iloc[i]['Chromosome']
            start = t.introns.iloc[i]['Start']
            end = t.introns.iloc[i]['End']

            #print(txid + "\t" + str(exonnumber) + "\t" + str(chromosome) + "\t" + str(start) + "\t" + str(end))

            iv = Interval(start - 1, end + 2)
            ivTree = IntervalTree()
            ivTree.add(iv)

            readIterator = simFile.readsInInterval(chromosome, start - 1, end + 2)

            exonintron = 0
            exonexon = 0
            intronexon = 0
            intron = 0
            exonexonfalse = 0
            garbage = 0

            for read in readIterator:

                # TODO: Implement check if read mapping actually reflects truth

                #truth = truthCollection.getTruth(read.name)
                #evaluation = read.evaluate(truth, junctions)
                start = read.absStart
                end = read.absEnd

                bamRead = read.bamRead

                color = ""

                if start <= iv.begin and end > iv.begin:

                    if read.splicing and not read.hasSpliceSite(ivTree):
                        exonexonfalse += 1
                        color = colors["exonexonfalse"]
                    elif read.splicing and read.hasSpliceSite(ivTree):
                        exonexon += 1
                        color = colors["exonexontrue"]
                    elif end < iv.end:
                        exonintron += 1
                        color = colors["exonintron"]
                    else :
                        #garbage += 1
                        color = colors["intron"]

                elif start < iv.end and end >= iv.end :
                    if read.splicing and not read.hasSpliceSite(ivTree):
                        exonexonfalse += 1
                        color = colors["exonexonfalse"]
                    elif read.splicing and read.hasSpliceSite(ivTree):
                        exonexon += 1
                        color = colors["exonexontrue"]
                    else :
                        intronexon += 1
                        color = colors["exonintron"]

                elif start >= iv.begin and end <= iv.end:
                    if read.splicing and not read.hasSpliceSite(ivTree):
                        exonexonfalse += 1
                        color = colors["exonexonfalse"]
                    elif read.splicing and read.hasSpliceSite(ivTree):
                        exonexon += 1
                        color = colors["exonexontrue"]
                    else :
                        intron +=1
                        color = colors["intron"]
                else:
                    # These are reads that start at base 1 of the exon
                    #garbage += 1
                    color = colors["intron"]

                bamRead.setTag('YC', color)

                classifiedReads.write(bamRead)

            print("\t".join([txid + "_" + str(exonnumber), str(exonintron), str(intron), str(intronexon), str(exonexon), str(exonexonfalse)]))
    classifiedReads.close()

else :

    junctions = extract_splice_sites(config["gene_gff"])

    truthCollection = TruthCollection(args.truthFile)

    truthCollection.readTruth()

    print("Done reading truth collection", file = sys.stderr)

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
