import pyranges as pr
import pandas as pd
import numpy as np
import math
import pysam
from collections import *
from utils import *
import logging


class Condition():
    def __init__(self, id, timepoint, ref, alt, conversion_rate, coverage ):
        self.id = id
        self.timepoint=timepoint
        self.ref=ref
        self.alt=alt
        self.conversion_rate=conversion_rate
        self.coverage = coverage
        self.bam="?"
    def __repr__(self):
        return self.__str__()  
    def __str__(self):
        ret = ("%s: [tp=%i, cr=%f, cv=%f]" % (self.id, self.timepoint, self.conversion_rate, self.coverage ) )
        return (ret)  
 
class Isoform():
    def __init__(self, id, fractions, splicing_status, is_labeled, transcript ):
        self.id = id
        self.fractions = fractions
        self.splicing_status = splicing_status
        self.t = transcript
        self.aln_block_len=0
        self.is_labeled=is_labeled
        if self.t is not None:
            self.strand = self.t.transcript.Strand
            # extract alignment blocks. Those contain absolute coords and are always ordered by genomic coords 
            # (e.g., [(16101295, 16101359), (16101727, 16101936), (16102490, 16102599), (16102692, 16102829), (16103168, 16103334), (16103510, 16103684), (16104365, 16104662)])
            self.aln_blocks=[]
            bstart=self.t.transcript.Start
            for idx, splicing in list(enumerate(self.splicing_status)):
                if splicing==1:
                    bend=self.t.introns.iloc[idx].Start-1 # last exonic 1-based pos
                    self.aln_blocks+=[(bstart, bend)]
                    self.aln_block_len+=(bend-bstart+1)
                    bstart=self.t.introns.iloc[idx].End+1 # 1st exonic 1-based pos
            bend=self.t.transcript.End
            self.aln_blocks+=[(bstart, bend)]
            self.aln_block_len+=(bend-bstart+1)
        #print("iso %s, alignment blocks: %s, len: %i" % (self.id, self.aln_blocks, self.aln_block_len))       
    def rel2abs_pos(self, rel_pos):
        """ convert isoform-relative coordinates to genomic coordinates.
            E.g., a rel_pos==0 will return the 1st position of the transcript in genomic coordinates (1-based) """
        ablocks = self.aln_blocks if self.strand == '+' else reversed(self.aln_blocks)
        abs_pos=None
        off = rel_pos
        block_id=0
        for bstart,bend in ablocks:
            bwidth = bend-bstart+1
            abs_pos = bstart if self.strand == '+' else bend
            if off - bwidth < 0:
                ret = abs_pos+off if self.strand == '+' else abs_pos-off
                return (ret, block_id)
            # next block
            off -= bwidth
            block_id+=1
        if off != 0:
            logging.warn("Invalid relative position %i / %i in isoform %s: %i" % (rel_pos, abs_pos, self, off))
        ret = abs_pos+off if self.strand == '+' else abs_pos-off
        return (ret, block_id)
    def block_contains(self, bstart, bend, abs_pos):
        return abs_pos>=bstart and abs_pos<=bend
    @DeprecationWarning
    def calc_cigar(self, abs_start, abs_end):
        """ calculate cigar """
        assert(abs_start<=abs_end), "abs start > abs_end. blocks: %s" % str(self.aln_blocks)
        ablocks = self.aln_blocks
        cigar=""
        gap_start=None
        started=False
        for bstart,bend in ablocks:
            bwidth = bend-bstart+1
            if self.block_contains(bstart, bend, abs_start):
                if self.block_contains(bstart, bend, abs_end): # start and stop in same block
                    return "%iM" % (abs_end-abs_start+1)
                cigar+="%iM" % (bend-abs_start+1)
                gap_start=bend+1
                started=True
            elif self.block_contains(bstart, bend, abs_end):
                if gap_start is not None:
                    cigar+="%iN" % (bstart-gap_start)
                cigar+="%iM" % (abs_end-bstart+1)
                return cigar
            elif started:
                cigar+="%iN%iM" % (bstart-gap_start, bwidth) # gap
                gap_start=bend+1                
        return cigar 
    def calc_cigartuples(self, abs_start, abs_end):      
        """ calculate cigar tuples (comparable to pysam rec.get_blocks() """
        assert(abs_start<=abs_end), "abs start > abs_end. blocks: %s" % str(self.aln_blocks)
        ablocks = self.aln_blocks
        cigartuples=[]
        started=False
        clen=0
        for bstart,bend in ablocks:
            bwidth = bend-bstart+1
            if self.block_contains(bstart, bend, abs_start):
                if self.block_contains(bstart, bend, abs_end): # start and stop in same block
                    cigartuples+=[(abs_start, abs_end)]
                    clen+=abs_end-abs_start+1
                    return cigartuples, clen    
                cigartuples+=[(abs_start, bend)]
                clen+=bend-abs_start+1
                started=True
            elif self.block_contains(bstart, bend, abs_end):
                cigartuples+=[(bstart, abs_end)]
                clen+=abs_end-bstart+1
                return cigartuples, clen    
            elif started:
                cigartuples+=[(bstart, bend)]
                clen+=bend-bstart+1
        return cigartuples, clen    
    def __repr__(self):
        return self.__str__()  
    def __str__(self):
        ret = ("%s_%s @ %s, [%s], [%s]" % (self.t.tid, self.id, self.t.region, ",".join(str(x) for x in self.fractions), ",".join(str(x) for x in self.splicing_status)) )
        return (ret)   

class Transcript():
    def __init__(self, config, tid, df, genome, conditions, max_ilen=None ):
        self.config = config
        self.tid = tid  
        self.is_valid = True
        self.df = df 
        self.transcript = self.df[self.df.Feature=='transcript'].iloc[0]
        self.region = to_region(self.transcript)
        self.len = self.transcript.End - self.transcript.Start + 1
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
                ilen = (ex.Start-1) - (last_end+1) + 1
                if (max_ilen is not None) and (ilen > max_ilen):
                    logging.warn("Skipping transcript %s due to too long intron (%i vs %i). Skipping transcript..." % ( self.tid, ilen, max_ilen))
                    self.is_valid = False
                    return
                idata.append(intron)
            last_end = ex.End
        if idata:
            self.introns=pd.concat(idata) 
        else:
            self.introns=pd.DataFrame()
        # set abundance and isoforms
        self.cond = conditions        
        self.abundance = math.ceil(self.config['transcripts'][tid]['abundance']) # abundance is float
        self.isoforms={}
        total_frac = [0] * len(self.cond)
        for iso in self.config['transcripts'][tid]['isoforms'].keys():
            frac=self.config['transcripts'][tid]['isoforms'][iso]['fractions']
            splicing_status=self.config['transcripts'][tid]['isoforms'][iso]['splicing_status'] if 'splicing_status' in self.config['transcripts'][tid]['isoforms'][iso] else []
            is_labeled=self.config['transcripts'][tid]['isoforms'][iso]['is_labeled'] if 'is_labeled' in self.config['transcripts'][tid]['isoforms'][iso] else 0
            if self.transcript.Strand == "-":
                splicing_status = list(reversed(splicing_status)) # so 1st entry in config refers to 1st intron.
            # assert len(splicing_status)==len(self.introns), "misconfigured splicing status for some isoform/conditions: %s (%i vs %i)" % ( self.tid, len(splicing_status), len(self.introns))
            if len(splicing_status)!=len(self.introns):
                logging.error("Misconfigured splicing status for some isoform/conditions of transcript %s (%i vs %i). Skipping transcript..." % ( self.tid, len(splicing_status), len(self.introns)))
                self.is_valid = False
                return
            self.isoforms[iso] = Isoform(iso, frac, splicing_status, is_labeled, self)    
            total_frac = [x + y for x, y in zip(total_frac, frac)]
        for x in total_frac:
            # assert with rounding error
            assert x <= 1.00001, "Total fraction of transcript > 1 (%s) for some isoform/conditions: %s" % ( ",".join(str(x) for x in total_frac), self.tid )
        # add mature, labeled isoform if not configured
        if 'mat' not in self.isoforms:
            frac = [(1.0 - x) for x in total_frac]
            iso='mat'
            splicing_status=[1] * len(self.introns)
            self.isoforms[iso] = Isoform(iso, frac, splicing_status, 1, self)    
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
            seq = reverse_complement(seq)
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
                logging.debug("corrected counts by %i" % (toadd))
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
    
    
    
class Model():
    """ builds a model using the passed configuration and conditions. 
        required config keys:
            conditions:   the configured experimental conditions
            gene_gff:     the GFF3 file containing all gencode-style annotations
            transcripts:  a dict with configured transcript configurations, mapping tids to gene_name, abundance and possible isoforms_+fractions (e.g., as created by splicing_simulator_config_creator)
            genome_fa:    FASTA file of the genome
            chrom_sizes:  optional file containing chrom sizes.
            max_ilen:     maximum intron length. Transcriptrs with a longer intron will be filetered out.
            
    """
    def __init__(self, config):
        self.config = config
        self.threads=config["threads"] if "threads" in config else 1
        # init conditions
        self.conditions=[]
        for id in config["conditions"].keys():
            self.conditions+=[Condition(id, 
                                        config["conditions"][id][0], 
                                        config["conditions"][id][1], 
                                        config["conditions"][id][2],
                                        config["conditions"][id][3],
                                        config["conditions"][id][4] )]
        # load + filter gene gff
        logging.info("Loading gene GFF")
        self.gff = pr.read_gff3(config["gene_gff"])
        self.gff = self.gff[self.gff.transcript_id.isin(list(config['transcripts'].keys()))] # reduce to contain only configured transcripts. 
        self.gff.Start = self.gff.Start+1# correct for pyranges bug?
        self.df = self.gff.df.sort_values(by=['Chromosome', 'Start','End']) # sort!
        # load genome
        logging.info("Loading genome")
        self.genome = pysam.FastaFile(config["genome_fa"])
        self.chrom_sizes = config["chrom_sizes"] if 'chrom_sizes' in config else config["genome_fa"]+".chrom.sizes"
        # instantiate transcripts
        self.transcripts=OrderedDict()
        self.max_ilen = config["max_ilen"] if 'max_ilen' in config else None
        self.readlen = config["readlen"] if 'readlen' in config else None
        for tid in list(config['transcripts'].keys()):
            tid_df = self.df[self.df['transcript_id']==tid] # extract data for this transcript only
            if tid_df.empty:
                logging.warn("No annotation data found for configured tid %s, skipping..." % (tid))
            else:
                t = Transcript(config, tid, tid_df, self.genome, self.conditions, max_ilen=self.max_ilen) 
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
    def write_isoform_truth(self, outdir):
        """ Write a TSV file with isoform truth data"""
        logging.info("Writing isoform truth data")
        out_file=outdir+"isoform_truth.tsv"
        if not files_exist(out_file+".gz"):
            with open(out_file, 'w') as out:
                print('\t'.join([str(x) for x in ['tid', 'iso', 'is_labeled', 'condition', 'fraction']]), file=out)
                for t in transcripts:
                    for iso in t.isoforms():
                        for i, cond in enumerate(t.cond):
                            print('\t'.join([str(x) for x in [t.tid, iso.id, iso.is_labeled, cond.id, iso.fractions[i]]]), file=out)
            bgzip(out_file, override=True, delinFile=True, index=True, threads=self.threads)
        out_file=out_file+".gz"
        return out_file    
        
    def save(self, f_model, recursion_limit=10000):
        """ save model to disk """
        sys.setrecursionlimit(recursion_limit)
        # save as picke
        with open(f_model, 'wb') as handle:
            pickle.dump((self), handle, protocol=pickle.HIGHEST_PROTOCOL)
        print("Model written to %s" % f_model)


# if __name__ == '__main__':
#     import json, pickle
#      
#     config_dir='/groups/ameres/Niko/projects/Ameres/splicing/splice_sim/testruns/big3_2/'
#     config = json.load(open(config_dir+'config.json'), object_pairs_hook=OrderedDict)
#     # read transcript data from external file if not in config
#     if 'transcripts' in config:
#         logging.info("Reading transcript configuration from config file.")
#     else:
#         assert "transcript_data" in config, "Transcript data needs to be configured either in config file ('transcripts' section) or in an external file referenced via 'transcript_data'"
#         tfile = config['transcript_data']
#         if not os.path.isabs(tfile):
#             tfile = config_dir+"/"+tfile
#         if not files_exist(tfile):
#             logging.info("Creating external transcript config file %s" % tfile)
#             calculate_transcript_data(config, config_dir, outdir)
#               
#         logging.info("Reading transcript configuration from external config file %s" % tfile)
#         tdata = json.load(open(tfile), object_pairs_hook=OrderedDict)
#         config["transcripts"]=tdata
#     m=Model(config)
#     print(m)
#     f_model = config_dir+"big3_2.model"
#     print("Model saved to " + f_model)
#     m.save(f_model)