'''
@author: niko.popitsch@imba.oeaw.ac.at
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


#============================================================================
# Creates config files from real data
#============================================================================
VERSION = "0.1"
logo = '''
===================================
Slamstr splice_sim config creator v%s
===================================
''' % VERSION

#============================================================================
usage = '''                           

  Copyright 2019 Niko Popitsch. All rights reserved.
  
  Licensed under the Apache License 2.0
  http://www.apache.org/licenses/LICENSE-2.0
  
  Distributed on an "AS IS" basis without warranties
  or conditions of any kind, either express or implied.

USAGE
'''
class Config(object):
    ODICT = "odict"
    def __init__(self):
        self.__dict__[self.ODICT] = OrderedDict()
    def __getattr__(self, item):
        return self.__dict__[self.ODICT][item]
    def __setattr__(self, key, value):
        self.__dict__[self.ODICT][key] = value
        
class Isoform():
    def __init__(self, id, fractions, splicing_status ):
        self.id = id
        self.fractions = fractions
        self.splicing_status = splicing_status
    def __repr__(self):
        return self.__str__()  
    def __str__(self):
        ret = ("%s, [%s], [%s]" % (self.id, ",".join(str(x) for x in self.fractions), ",".join(str(x) for x in self.splicing_status)) )
        return (ret)   
 

#
# This method tries to select the isoform from 'isoforms' that best supports the input read counts ('pre','mat').
# It then calculates the remaining (theoretical) read counts by subtracting the 'used' reads from the input read couts.
# for this, it first calculates which intron provides the min. fraction ('min__frac') and then subtracts relative read counts,
# so the intron with 'min_frac' will get 0 remaining pre or mat reads (depending on the splice status of the selected isoform) and 
# all other intros will keep some positive read counts depending on 'min_frac'. 
#
# Finally, min_frac is used to calculate a fraction using the 'remaining_frac' value which is at the beginning of the recursion 1.0
# but is gradually reduced as introns are selected. In theory, 'remaining_frac' should shrink to zero if all reads can be explained by 
# the given isoforms but often its not. So the final fractions are then calculated relative to 'remaining_frac'.
#
def pick_iso_with_max_support(isoforms, mat, pre, remaining_frac, informative_intron_mask, debug=False):
    if debug:
        print('mat: [%s]' % ', '.join(map(str, mat)))
        print('pre: [%s]' % ', '.join(map(str, pre)))
    new_mat=mat
    new_pre=pre
    f_mat=list(x/(x+y) if (x+y)>0 else 0 for x,y, in zip (mat, pre))
    min_frac = 0
    picked_iso=None
    max_support_sum=0
    for iso in isoforms:
        support = [x if y==1 else 1.0-x for x,y in zip(f_mat, iso.splicing_status) ]
        # remove uninformative introns
        for i in range(len(informative_intron_mask),0,-1):
            if not informative_intron_mask[i-1]:
                del support[i-1]
                i=i-1
        support_sum=sum(support)
        if debug:
            print('%s: [%s] sum=%f' % (iso.id, ', '.join(map(str, support)), support_sum) )
            print(support)
        # pick this isoform for now if support is greater than before and its supported by all introns
        if support_sum>=max_support_sum and min(support) > 0: # this will pick the isoform with max splicing status if equal support sum
            picked_iso = iso
            max_support_sum=support_sum
            min_frac = min(support)
    if min_frac == 0:
        picked_iso=None
    else:
        # calc new read numbers
        new_mat=[m if s==0 else m-min_frac*(m+p) for s,m,p in zip(picked_iso.splicing_status, mat, pre)]
        new_pre=[p if s==1 else p-min_frac*(m+p) for s,m,p in zip(picked_iso.splicing_status, mat, pre)]
    frac = remaining_frac * min_frac
    if debug:
        print('picked: [%s] frac=%f' % (picked_iso.id if picked_iso else "None", frac) )
        
    return (picked_iso, frac, new_mat, new_pre)
    
# calc isoform data 
def calc_iso(tid, transcript, introns, rnk, conditions, times):
    #print(tid)
        
    # calculate isoforms
    isoforms=[]
    for i in range(0,int(rnk)):
        splicing_status = [1] * int(i) + [0] * (int(rnk)-i)
        name = "in"+str(i) if i>0 else "pre"
        isoforms+=[Isoform(name, [], splicing_status)]
    isoforms+=[Isoform('mat', [], [1] * rnk)]
    
    # find uninformative introns (with low EE counts in all conditions)
    sum_informative_mat=None
    for cond in conditions:
        EE = list(introns['EE_TC_'+cond])
        if sum_informative_mat:
            sum_informative_mat = [x+s for x,s in zip(EE,sum_informative_mat)]
        else:
            sum_informative_mat = EE
    informative_intron_mask = [1 if s >= 5 else 0 for s in sum_informative_mat]   
    if sum(informative_intron_mask)==0:
        print("Not enough EE counts for any intron of %s." % (tid))
        return (OrderedDict())
    
    for cond in conditions:
        EE = list(introns['EE_TC_'+cond])
        EI = list(introns['EI_TC_'+cond])
        IE = list(introns['IE_TC_'+cond])
        # calc numbers of reads supporting premature/mature form per intron. 
        # add pseudocount of 1 to each intron
        mat=[x+1 for x in EE]
        pre=list( (x+y)/2+1 for x,y in zip(EI,IE) )
        
        # test (excel)
        # mat=[12 ,   3  ,  2]
        # pre=[4   , 15  ,  18]
        # test Prdx1 30min
        # mat=[0 ,  6 ,  16 , 8, 13]
        # pre=[2.5, 7.5, 3.5, 3.0, 2.5]
        remaining_iso = isoforms.copy()
        iso2frac={}
        sum_frac=0
        remaining_frac=1.0
        debug=False
        #debug=True if tid=='ENSMUST00000112580.7' else False
        if debug:
            print("==========================\n%s\n==========================" % cond )
            print(isoforms)
            #print(sum_informative_mat)
            print(informative_intron_mask)
        while True:
            picked_iso, frac, mat, pre = pick_iso_with_max_support(remaining_iso, mat, pre, remaining_frac, informative_intron_mask, debug=debug)
            if picked_iso:
                remaining_iso.remove(picked_iso)
                iso2frac[picked_iso] = frac
                sum_frac+=frac
                remaining_frac-=frac
            else:
                break

        for iso in isoforms:
            frac = iso2frac[iso]/sum_frac if iso in iso2frac else 0.0
            iso.fractions.append(frac)
            if debug:
             print("%s: %f" % (iso.id, frac))

    ret = OrderedDict()
    for iso in isoforms:
        ret[iso.id] = OrderedDict()
        ret[iso.id]['splicing_status']=iso.splicing_status
        ret[iso.id]['fractions']=iso.fractions
    return (ret)


parser = ArgumentParser(description=usage, formatter_class=RawDescriptionHelpFormatter)
parser.add_argument("-conf","--config", type=existing_file, required=True, dest="confF", help="JSON config file")
parser.add_argument("-m","--mode", type=str, required=True, dest="mode", help="mode: 'from_slamstr', '1:1'")
parser.add_argument("-o","--out", type=str, required=False, dest="outF", metavar="outF", help="Output file (otional). If not set, results will be written to file configured in 'transcript_data' config property.")
args = parser.parse_args()   
startTime = time.time()
print(logo)

    
# load + check config
#config = json.load(open('/Users/niko.popitsch/eclipse-workspace/slamstr/src/test/splicing_simulator/testconfig_data_creator.json'), object_pairs_hook=OrderedDict)
config = json.load(open(args.confF), object_pairs_hook=OrderedDict)
if "random_seed" in config:
    random.seed(config["random_seed"])

# outoput file
outF = args.outF if args.outF else config['transcript_data']
if not os.path.isabs(outF):
    outF = os.path.abspath(outF)      
outdir =  os.path.dirname(outF)                       
if not os.path.exists(outdir):
    print("Creating dir " + outdir)
    os.makedirs(outdir)

# get readlen (needed for abundance calc)
readlen=config["readlen"]

# get conditions and labeling times
conditions = list(config["conditions"].keys())
times = list(config["conditions"].values())
out=OrderedDict()

# ------------------------------------------------------------------------------------
#    from slamstr data
# ------------------------------------------------------------------------------------
if args.mode == 'from_slamstr':
    # get abundance data
    ab_table=pd.read_csv(config["slamstr_abundance_data_table"],delimiter='\t',encoding='utf-8')
    TAB_ALL=list(ab_table["TAB_ALL_"+conditions[-1]])
    transcript_len=list(ab_table["len"])
    abundance_multiplier = config["abundance_multiplier"] if 'abundance_multiplier' in config else 1.0
    abundances = [(x*readlen*abundance_multiplier)/y for x, y in zip(TAB_ALL, transcript_len)]
    min_abundance = config["min_abundance"] if 'min_abundance' in config else 0.0
    print("minimum abundance: %f" % (min_abundance) )
    
    # create output
    transcripts = pd.read_csv(config["slamstr_transcript_table"],delimiter='\t',encoding='utf-8')
    introns = pd.read_csv(config["slamstr_intron_table"],delimiter='\t',encoding='utf-8')
    for i, t in transcripts.iterrows():
        if abundances[i] > min_abundance:
            tid = t['transcript_id']
            rnk = None if math.isnan(t['rnk']) else int(t['rnk']) 
            isoform_data=calc_iso(tid, t, introns[introns['transcript_id']==tid], rnk, conditions, times) if rnk else OrderedDict()
            if rnk and len(isoform_data) == 0:
                # not enough data. skip
                print("Skipped %s due to too-low number of informative reads" % (tid) )
            else:
                out[tid] = OrderedDict()
                out[tid]["gene_name"] = t['gene_name']
                out[tid]["abundance"] = abundances[i]
                out[tid]["isoforms"] = isoform_data
        else:
            print("Skipped %s due to too-low abundance." % (tid) )
# ------------------------------------------------------------------------------------
#    1:1 mat/pre mix with fixed abundance for all transcripts
# ------------------------------------------------------------------------------------
elif args.mode == '1:1':
    # create output
    transcripts = pd.read_csv(config["slamstr_transcript_table"],delimiter='\t',encoding='utf-8')
    introns = pd.read_csv(config["slamstr_intron_table"],delimiter='\t',encoding='utf-8')
    abundances = [10] * len ( transcripts.index )
    min_abundance = config["min_abundance"] if 'min_abundance' in config else 0.0
    for i, t in transcripts.iterrows():
        if abundances[i] > min_abundance:
            tid = t['transcript_id']
            rnk = None if math.isnan(t['rnk']) else int(t['rnk']) 
            isoform_data=OrderedDict()
            isoform_data['pre'] = OrderedDict()
            if rnk:
                isoform_data['pre']['splicing_status']=[0] * rnk
            isoform_data['pre']['fractions']=[0.5] * len(times)
            # output data
            out[tid] = OrderedDict()
            out[tid]["gene_name"] = t['gene_name']
            out[tid]["abundance"] = abundances[i]
            out[tid]["isoforms"] = isoform_data
        else:
            print("Skipped %s due to too-low abundance." % (tid) )
    
else:
    print("Unknown mode %s " % (args.mode))  
    sys.exit(1)
    
with open(outF, 'w') as config_file:
    json.dump(out, config_file, indent=2)
    

print('Done. Results written to %s' % (outF))