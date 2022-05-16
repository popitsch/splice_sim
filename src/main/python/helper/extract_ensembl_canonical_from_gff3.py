from argparse import ArgumentParser, RawDescriptionHelpFormatter
import os, sys
import gzip
import pandas as pd

#
# Extract Ensembl canonical tids.tsv from Gencode gff3
#
usage = '''                           

  Copyright (C) 2021 XXX.  All rights reserved.

  Distributed on an "AS IS" basis without warranties
  or conditions of any kind, either express or implied.

USAGE
'''
if __name__ == '__main__':
    parser = ArgumentParser(description=usage, formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument("-g", "--gff3", type=str, required=True, dest="gff3_file",
                        metavar="gff3_file", help="Gencode GFF3 file")
    parser.add_argument("-o", "--outfile", type=str, required=False, dest="outfile", metavar="outfile",
                        help="output file")
    args = parser.parse_args()

    opener = gzip.open if ".gz" in args.gff3_file else open

    if args.outfile:
        fout = open(args.outfile, 'w')
    else :
        fout = sys.stdout

    print("transcript_id", fout)

    with opener(args.gff3_file, 'rt') as f:
        for line in f:
            if line.startswith("#") :
                continue
            chr, source, feature, start, end, score, strand, phase, info = line.rstrip().split("\t")
            if (feature == "transcript"):
                tags = dict(x.split("=") for x in info.split(";"))
                if "Ensembl_canonical" in tags["tag"] :
                    print(tags['transcript_id'], fout)

    if args.outfile:
        fout.close()