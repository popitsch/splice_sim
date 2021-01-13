from argparse import ArgumentParser, RawDescriptionHelpFormatter
import os, sys
import pandas as pd

usage = '''                           

  Copyright (C) 2021 XXX.  All rights reserved.

  Distributed on an "AS IS" basis without warranties
  or conditions of any kind, either express or implied.

USAGE
'''
if __name__ == '__main__':
    parser = ArgumentParser(description=usage, formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument("-i", "--slamstr_transcript_file", type=str, required=True, dest="slamstr_transcript_file", metavar="slamstr_transcript_file", help="slamstr_transcript_file")
    parser.add_argument("-m", "--min_", type=int, required=True, dest="min_tab_180", help="min value in the TAB_ALL_tpm_180min column")
    parser.add_argument("-o", "--outfile", type=str, required=False, dest="outfile", metavar="outfile", help="output file")

    m=args.min_tab_180
    d=pd.read_csv(args.slamstr_transcript_file,delimiter='\t',encoding='utf-8')
    d=d[d['TAB_ALL_tpm_180min']>=m]
    pd.DataFrame(d['transcript_id']).to_csv(args.outfile, '\t', header=True, index=False)
    