"""

Transcriptional model for splice_sim containing transcripts, isoforms and conditions.
Usage:
1) build a model from a config file:
    m = Model.build_model( config_file )
2) save model to disk:
    m.save(f_model)
3) load model from disk:
    m =Model.load_from_file(f_model)
    
"""
import pyranges as pr
import pandas as pd
import numpy as np
import math
import pysam
from collections import OrderedDict, Counter
import os, sys
import pickle
import logging

def localize_config(config):
    """ Recursively iterates json and adds /Volumes prefix to /group paths"""
    for k, v in config.items():
        if isinstance(v, str) and v.startswith('/groups'):
            config[k]='/Volumes'+ config[k]
        elif isinstance(v, dict):
            config[k]=localize_config(v)
        elif isinstance(v, list):
            config[k]=['/Volumes'+x if x.startswith('/groups') else '' for x in v]
    return config

def calculate_transcript_data(config, config_dir, outdir):
    """ Calculate a data config file """ 
    mode = config['isoform_mode']      
    logging.info("Calculating transcript configuration with mode %s" % mode)
    tdata=OrderedDict()
    # ------------------------------------------------------------------------------------
    #    1:1 mat/pre mix with fixed abundance for all transcripts.
    #    input: transcript_ids table with column transcript_id
    #    config: 
    # ------------------------------------------------------------------------------------
    if mode == '1:1':
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
        frac_old_mature=min(1.0,config["frac_old_mature"]) if 'frac_old_mature' in config else 0.0 # fraction of 'old' (unlabeled), mature isoform
        for i, t in transcripts.iterrows():
            if abundances[i] > min_abundance:
                tid = t['transcript_id']
                df = gff_df[gff_df['transcript_id']==tid] 
                rnk=len(df[df.Feature=='exon'].index)-1
                gene_name = df[df.Feature=='transcript'].iloc[0]['gene_name']
                isoform_data=OrderedDict()
                if frac_old_mature>0: # add 'old' rna isoform
                    isoform_data['old'] = OrderedDict()
                    isoform_data['old']['splicing_status']=[1] * rnk
                    isoform_data['old']['fractions']=[frac_old_mature]
                    isoform_data['old']['is_labeled']=[0]
                if frac_old_mature < 1.0:
                    isoform_data['pre'] = OrderedDict()
                    if rnk: # at least one intron
                        isoform_data['pre']['splicing_status']=[0] * rnk
                    isoform_data['pre']['fractions']=[round((1-frac_old_mature) * 0.5,3)]
                    isoform_data['pre']['is_labeled']=[1]
                # output data
                tdata[tid] = OrderedDict()
                tdata[tid]["gene_name"] = gene_name
                tdata[tid]["abundance"] = abundances[i]
                tdata[tid]["isoforms"] = isoform_data
            else:
                logging.warn("Skipped %s due to too-low abundance." % (tid) )
    else:
        logging.error("Unknown mode %s " % (args.mode))  
        sys.exit(1)    
    # Write JSON file
    outF = config['transcript_data']
    if not os.path.isabs(outF):
        outF = config_dir + outF
    with open(outF, 'w') as out:
        json.dump(tdata, out, indent=2)     
    logging.info('Done. Results written to %s' % (outF))

class Condition():
    def __init__(self, name, ref, alt, conversion_rate, base_coverage ):
        self.name=name
        self.ref=ref
        self.alt=alt
        self.conversion_rate=conversion_rate
        self.base_coverage = base_coverage
    def __repr__(self):
        return self.__str__()  
    def __str__(self):
        ret = ("condition %s: %s>%s [cr=%f, cv=%f]" % (self.name, self.ref, self.alt, self.conversion_rate, self.base_coverage ) )
        return (ret)  
 
class Isoform():
    def __init__(self, id, fractions, splicing_status, is_labeled, transcript ):
        self.id = id
        self.fractions = fractions
        self.splicing_status = splicing_status
        self.transcript = transcript
        self.aln_block_len=0
        self.is_labeled=is_labeled
        if self.transcript is not None:
            self.strand = self.transcript.Strand
            # extract alignment blocks. Those contain absolute coords and are always ordered by genomic coords 
            # (e.g., [(16101295, 16101359), (16101727, 16101936), (16102490, 16102599), (16102692, 16102829), (16103168, 16103334), (16103510, 16103684), (16104365, 16104662)])
            self.aln_blocks=[]
            bstart=self.transcript.Start
            for idx, splicing in list(enumerate(self.splicing_status)):
                if splicing==1:
                    bend=self.transcript.introns[idx]['Start']-1 # last exonic 1-based pos
                    self.aln_blocks+=[(bstart, bend)]
                    self.aln_block_len+=(bend-bstart+1)
                    bstart=self.transcript.introns[idx]['End']+1 # 1st exonic 1-based pos
            bend=self.transcript.End
            self.aln_blocks+=[(bstart, bend)]
            self.aln_block_len+=(bend-bstart+1)
        #print("iso %s, alignment blocks: %s, len: %i" % (self.id, self.aln_blocks, self.aln_block_len))       
    def block_contains(self, bstart, bend, abs_pos):
        return abs_pos>=bstart and abs_pos<=bend 
    def __repr__(self):
        return self.__str__()  
    def __str__(self):
        ret = ("%s_%s @ %s, [%s], [%s]" % (self.transcript.tid, self.id, self.transcript.region, ",".join(str(x) for x in self.fractions), ",".join(str(x) for x in self.splicing_status)) )
        return (ret)   
    
class Transcript():
    def __init__(self, config, tid, df, genome, condition, max_ilen=None ):
        self.config = config
        self.tid = tid  
        self.is_valid = True
        self.cond=condition
        self.transcript_data = df[df.Feature=='transcript'].iloc[0].to_dict()
        self.transcript_id=self.transcript_data['transcript_id']
        self.Chromosome=self.transcript_data['Chromosome']
        self.Start=self.transcript_data['Start']
        self.End=self.transcript_data['End']
        self.Strand=self.transcript_data['Strand']
        self.region = "%s:%i-%i" % (self.Chromosome,self.Start,self.End)
        self.len = self.End - self.Start + 1
        self.transcript_seq = genome.fetch(region=self.region)
        self.exons = []
        self.introns = []
        last_end = -1
        for index, ex in df[df.Feature=='exon'].iterrows():
            exon=ex.to_dict()
            self.exons.append(exon)
            if last_end > 0: # skip 1st
                ilen=(exon['Start']-1) - (last_end+1) + 1
                intron={ 'transcript_id': self.transcript_id,
                      'Feature': 'intron',
                      'Chromosome':self.Chromosome, 
                      'Strand':self.Strand,
                      'Start':last_end+1,
                      'End':exon['Start']-1,
                      'exon_number': int(exon['exon_number']),
                      'len': ilen
                      }
                if (max_ilen is not None) and (ilen > max_ilen):
                    logging.warn("Skipping transcript %s due to too long intron (%i vs %i). Skipping transcript..." % ( self.tid, ilen, max_ilen))
                    self.is_valid = False
                    return
                self.introns.append(intron)
            last_end = exon['End']
        # set abundance and isoforms
        self.abundance = math.ceil(self.config['transcripts'][tid]['abundance']) # abundance is float
        self.isoforms=OrderedDict()
        total_frac = [0]
        for iso in self.config['transcripts'][tid]['isoforms'].keys():
            frac=self.config['transcripts'][tid]['isoforms'][iso]['fractions']
            splicing_status=self.config['transcripts'][tid]['isoforms'][iso]['splicing_status'] if 'splicing_status' in self.config['transcripts'][tid]['isoforms'][iso] else []
            is_labeled=self.config['transcripts'][tid]['isoforms'][iso]['is_labeled'] if 'is_labeled' in self.config['transcripts'][tid]['isoforms'][iso] else 1
            if self.Strand == "-":
                splicing_status = list(reversed(splicing_status)) # so 1st entry in config refers to 1st intron.
            if len(splicing_status)!=len(self.introns):
                logging.error("Misconfigured splicing status for some isoform of transcript %s (%i vs %i). Skipping transcript..." % ( self.tid, len(splicing_status), len(self.introns)))
                self.is_valid = False
                return
            self.isoforms[iso] = Isoform(iso, frac, splicing_status, is_labeled, self)    
            total_frac = [x + y for x, y in zip(total_frac, frac)]
        for x in total_frac:
            # assert with rounding error
            assert x <= 1.00001, "Total fraction of transcript > 1 (%s) for some isoform: %s" % ( ",".join(str(x) for x in total_frac), self.tid )
        # add mature, labeled isoform if not configured
        if 'mat' not in self.isoforms:
            frac = [(1.0 - x) for x in total_frac]
            iso='mat'
            splicing_status=[1] * len(self.introns)
            self.isoforms[iso] = Isoform(iso, frac, splicing_status, 1, self) # add labeled, mature isoform
        #print('added mature isoform for %s: %s' % (tid, self.isoforms[iso]) )  
    def get_dna_seq(self, splicing_status):
        seq = self.transcript_seq
        for idx, splicing in reversed(list(enumerate(splicing_status))):
            if splicing==1:
                pos1=self.introns[idx]['Start'] - self.Start
                pos2=self.introns[idx]['End'] - self.Start + 1
                seq = seq[:pos1]+seq[pos2:]
        return (seq)
    def get_rna_seq(self, iso):
        seq = self.get_dna_seq(self.isoforms[iso].splicing_status)
        if self.Strand == '-':
            seq = reverse_complement(seq)
        return (seq)
    def get_sequences(self):
        ret=OrderedDict()
        ret[cond]=OrderedDict()
        total_count=0
        for id in self.isoforms.keys():
            frac = self.isoforms[id].fractions[0]
            count=int(self.abundance * frac)
            total_count+=count
            ret[cond][id]=[
                self.get_rna_seq(id), 
                count]
        toadd = self.abundance - total_count
        #assert toadd==0 | toadd==1 , "rounding error"
        if toadd > 0: # add 1 to last isoform
            logging.debug("corrected counts by %i" % (toadd))
            id = list(self.isoforms.keys())[-1]
            ret[cond][id]=[ret[cond][id][0], ret[cond][id][1] + toadd]
        return (ret)
    def to_bed(self):
        return ("%s\t%i\t%i\t%s\t%i\t%s" % (self.Chromosome, self.Start, self.End, self.transcript_data['ID'], 1000, self.Strand) )
    def __repr__(self):
        return self.__str__()  
    def __str__(self):
        ret = ("%s" % (self.transcript_id) )
        return (ret)   

class Model():
    """ builds a model using the passed configuration. 
        required config keys:
            condition:   the configured experimental condition
            gene_gff:     the GFF3 file containing all gencode-style annotations
            transcripts:  a dict with configured transcript configurations, mapping tids to gene_name, abundance and possible isoforms_+fractions (e.g., as created by splicing_simulator_config_creator)
            genome_fa:    FASTA file of the genome
            chrom_sizes:  optional file containing chrom sizes.
            max_ilen:     maximum intron length. Transcriptrs with a longer intron will be filetered out.
            
    """
    def __init__(self, config):
        self.config = config
        self.threads=config["threads"] if "threads" in config else 1
        # init condition
        self.condition=Condition(
            config["condition"]["name"], 
            config["condition"]["ref"], 
            config["condition"]["alt"], 
            config["condition"]["conversion_rate"], 
            config["condition"]["base_coverage"])
        # load + filter gene gff
        logging.info("Loading gene GFF")
        self.gff = pr.read_gff3(config["gene_gff"])
        self.gff = self.gff[self.gff.transcript_id.isin(list(config['transcripts'].keys()))] # reduce to contain only configured transcripts. 
        self.gff.Start = self.gff.Start+1# correct for pyranges bug?
        df = self.gff.df.sort_values(by=['Chromosome', 'Start','End']) # sort!
        # load genome
        logging.info("Loading genome")
        genome = pysam.FastaFile(config["genome_fa"])
        self.chrom_sizes = config["chrom_sizes"] if 'chrom_sizes' in config else config["genome_fa"]+".chrom.sizes"
        # instantiate transcripts
        self.transcripts=OrderedDict()
        self.max_ilen = config["max_ilen"] if 'max_ilen' in config else None
        self.readlen = config["readlen"] if 'readlen' in config else None
        for tid in list(config['transcripts'].keys()):
            tid_df = df[df['transcript_id']==tid] # extract data for this transcript only
            if tid_df.empty:
                logging.warn("No annotation data found for configured tid %s, skipping..." % (tid))
            else:
                t = Transcript(config, tid, tid_df,genome, self.condition, max_ilen=self.max_ilen) 
                if t.is_valid:
                    self.transcripts[tid] = t
        logging.info("Instantiated %i transcripts" % len(self.transcripts))      
    def write_gff(self, outdir):
        """ Write a filtered GFF file containing all kept/simulated transcripts """
        logging.info("Writing filtered gene GFF")
        out_file=outdir+"gene_anno.gff3"
        d=pd.read_csv(self.config["gene_gff"],delimiter='\t',encoding='utf-8')
        if not files_exist(out_file+".gz"):
            with open(out_file, 'w') as out:
                for index, row in d.iterrows():
                    keep=False
                    for k in row[8].split(';'):
                        if k.startswith('transcript_id='):
                            keep=k[len('transcript_id='):] in self.transcripts.keys()
                    if keep:
                        print('\t'.join(str(x) for x in row), file=out)
            bgzip(out_file, override=True, delinFile=True, index=True, threads=self.threads)
        out_file=out_file+".gz"
        return out_file
    def write_transcript_md(self, outdir):
        """ Write a filtered TSV file containing md about all kept/simulated transcripts """
        logging.info("Writing filtered gene TSV")
        out_file=outdir+"gene_anno.tsv"
        headers=None
        if not files_exist(out_file+".gz"):
            with open(out_file, 'w') as out:
                for t in self.transcripts.values():       
                    if headers is None:
                        headers=[str(h) for h in t.df.columns if h not in ['Source', 'Feature', 'Score', 'Frame', 'ID', 'havana_gene', 'Parent', 'exon_number', 'exon_id']]
                        print('\t'.join(headers), file=out)            
                    print('\t'.join([str(t.transcript[h]) for h in headers]), file=out)
            bgzip(out_file, override=True, delinFile=True, index=False, threads=self.threads)
        out_file=out_file+".gz"
        return out_file
    def write_isoform_truth(self, outdir):
        """ Write a TSV file with isoform truth data"""
        logging.info("Writing isoform truth data")
        out_file=outdir+"isoform_truth.tsv"
        if not files_exist(out_file+".gz"):
            with open(out_file, 'w') as out:
                print('\t'.join([str(x) for x in ['tid', 'iso', 'is_labeled', 'condition_name', 'true_abundance', 'true_fraction', 'splicing_status']]), file=out)
                for t in self.transcripts.values():
                    splice_status=','.join([str(y) for y in iso.splicing_status]) if len(iso.splicing_status)>0 else 'NA'
                    print('\t'.join([str(x) for x in [t.tid, iso.id, iso.is_labeled, self.condition.name, t.abundance, iso.fractions[i], splice_status]]), file=out)
            bgzip(out_file, override=True, delinFile=True, index=False, threads=self.threads)
        out_file=out_file+".gz"
        return out_file    
    def save(self, f_model, recursion_limit=10000):
        """ save model to disk """
        sys.setrecursionlimit(recursion_limit)
        # save as picke
        with open(f_model, 'wb') as handle:
            pickle.dump((self), handle, protocol=pickle.HIGHEST_PROTOCOL)
        print("Model written to %s" % f_model)
    def __repr__(self):
        return self.__str__()  
    def __str__(self):
        ret = ("Model [t: %i]" % (len(self.transcripts) ) )
        return (ret)
    @classmethod
    def load_from_file(cls, model_file):
        return pickle.load( open( model_file, "rb" ) )
    @classmethod
    def build_model(cls, config_file, out_dir=None):
        """ build transcriptional model and save """
        # load config 
        config = json.load(open(config_file), object_pairs_hook=OrderedDict)
        if config_file.startswith('/Volumes'): # localize file paths
            config=localize_config(config)                    
        # ensure out_dir
        if out_dir is None:
            out_dir=os.path.dirname(config_file)+'/'
        if not os.path.exists(out_dir):
            os.mkdirs(out_dir)
        # read transcript data from external file if not in config
        if 'transcripts' in config:
            logging.info("Reading transcript configuration from config file.")
        else:
            assert "transcript_data" in config, "Transcript data needs to be configured either in config file ('transcripts' section) or in an external file referenced via 'transcript_data'"
            config_dir = os.path.dirname(config_file)+'/'
            tfile = config['transcript_data']
            if not os.path.isabs(tfile):
                tfile = config_dir+"/"+tfile
            if not os.path.exists(tfile):
                logging.info("Creating external transcript config file %s" % tfile)
                calculate_transcript_data(config, config_dir, out_dir)               
            logging.info("Reading transcript configuration from external config file %s" % tfile)
            tdata = json.load(open(tfile), object_pairs_hook=OrderedDict)
            config["transcripts"]=tdata
        # build model
        m=cls(config)
        f_model = out_dir+'/'+config["dataset_name"]+".model"
        print("Model saved to " + f_model)
        m.save(f_model)
        return m, f_model

# if __name__ == '__main__':
#     import json, pickle
#     m, f_model = Model.build_model('/Volumes/groups/ameres/Niko/projects/Ameres/splicing/splice_sim/testruns/nf1/config.json')
#     print(m)
#     m2 = Model.load_from_file(f_model)
#     print(m2) 
    
    