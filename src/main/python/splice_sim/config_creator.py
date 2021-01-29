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
from model import Isoform


#============================================================================

def pick_iso_with_max_support(isoforms, mat, pre, remaining_frac, informative_intron_mask, debug=False):
    """ This method tries to select the isoform from 'isoforms' that best supports the input read counts ('pre','mat').
    It then calculates the remaining (theoretical) read counts by subtracting the 'used' reads from the input read couts.
    for this, it first calculates which intron provides the min. fraction ('min__frac') and then subtracts relative read counts,
    so the intron with 'min_frac' will get 0 remaining pre or mat reads (depending on the splice status of the selected isoform) and 
    all other intros will keep some positive read counts depending on 'min_frac'. 
    
    Finally, min_frac is used to calculate a fraction using the 'remaining_frac' value which is at the beginning of the recursion 1.0
    but is gradually reduced as introns are selected. In theory, 'remaining_frac' should shrink to zero if all reads can be explained by 
    the given isoforms but often its not. So the final fractions are then calculated relative to 'remaining_frac'. """

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

#============================================================================
def calc_iso(tid, transcript, introns, rnk, conditions, times):
    """ Calculate isoform data"""
    #print(tid)
        
    # calculate isoforms
    isoforms=[]
    for i in range(0,int(rnk)):
        splicing_status = [1] * int(i) + [0] * (int(rnk)-i)
        name = "in"+str(i) if i>0 else "pre"
        isoforms+=[Isoform(name, [], splicing_status, None)]
    isoforms+=[Isoform('mat', [], [1] * rnk, None)]
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


def calculate_transcript_data(config, config_dir, outdir):
    """ Calculate a data config file """
    # output file
    outF = config['transcript_data']
    if not os.path.isabs(outF):
        outF = config_dir + outF
        
    mode = config['isoform_mode']      
    logging.info("Calculating transcript configuration with mode %s" % mode)
    # get readlen (needed for abundance calc)
    readlen=config["readlen"]
    
    # get conditions and labeling times
    conditions = list(config["conditions"].keys())
    times = list(config["conditions"].values())
    out=OrderedDict()
    
    # ------------------------------------------------------------------------------------
    #    from slamstr data
    # ------------------------------------------------------------------------------------
    if mode == 'from_slamstr':
        assert 'slamstr_abundance_data_table' in config, "For from_slamstr mode, a slamstr_abundance_data_table must be provided in a file referenced by the 'slamstr_abundance_data_table' config property"
        assert 'slamstr_transcript_table' in config, "For from_slamstr mode, a slamstr_transcript_table must be provided in a file referenced by the 'slamstr_transcript_table' config property"
        assert 'slamstr_intron_table' in config, "For from_slamstr mode, a slamstr_intron_table must be provided in a file referenced by the 'slamstr_intron_table' config property"
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
                logging.warn("Skipped %s due to too-low abundance." % (tid) )
    # ------------------------------------------------------------------------------------
    #    1:1 mat/pre mix with fixed abundance for all transcripts.
    #    input: slamstr_transcript_table with columns transcript_id, gene_name, rnk (=#introns)
    # ------------------------------------------------------------------------------------
    elif mode == '1:1_slamstr':
        # create output
        assert 'slamstr_transcript_table' in config, "For 1:1_slamstr mode, a slamstr_transcript_table must be provided in a file referenced by the 'slamstr_transcript_table' config property"
        transcripts = pd.read_csv(config["slamstr_transcript_table"],delimiter='\t',encoding='utf-8')
        base_abundance = config["base_abundance"] if 'base_abundance' in config else 10
        abundances = [base_abundance] * len ( transcripts.index )
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
                logging.warn("Skipped %s due to too-low abundance." % (tid) )
    # ------------------------------------------------------------------------------------
    #    1:1 mat/pre mix with fixed abundance for all transcripts.
    #    input: transcript_ids table with column transcript_id
    #    config: 
    # ------------------------------------------------------------------------------------
    elif mode == '1:1':
        assert 'transcript_ids' in config, "For 1:1 mode, a list of transcript ids must be provided in a file referenced by the 'transcript_ids' config property"
        # load GFF
        gff = pr.read_gff3(config["gene_gff"])
        gff_df = gff.df
        # create output
        tfile = config["transcript_ids"]
        if not os.path.isabs(tfile):
            tfile = config_dir + tfile
        transcripts = pd.read_csv(tfile,delimiter='\t',encoding='utf-8')
        logging.info("read %i tids" % len(transcripts))
        base_abundance = config["base_abundance"] if 'base_abundance' in config else 10
        abundances = [base_abundance] * len ( transcripts.index )
        min_abundance = config["min_abundance"] if 'min_abundance' in config else 0.0
        for i, t in transcripts.iterrows():
            if abundances[i] > min_abundance:
                tid = t['transcript_id']
                df = gff_df[gff_df['transcript_id']==tid] 
                rnk=len(df[df.Feature=='exon'].index)-1
                gene_name = df[df.Feature=='transcript'].iloc[0]['gene_name']
                isoform_data=OrderedDict()
                isoform_data['pre'] = OrderedDict()
                if rnk:
                    isoform_data['pre']['splicing_status']=[0] * rnk
                isoform_data['pre']['fractions']=[0.5] * len(times)
                # output data
                out[tid] = OrderedDict()
                out[tid]["gene_name"] = gene_name
                out[tid]["abundance"] = abundances[i]
                out[tid]["isoforms"] = isoform_data
            else:
                logging.warn("Skipped %s due to too-low abundance." % (tid) )
    else:
        logging.error("Unknown mode %s " % (args.mode))  
        sys.exit(1)     
    with open(outF, 'w') as config_file:
        json.dump(out, config_file, indent=2)
    logging.info('Done. Results written to %s' % (outF))
