import pyranges as pr
import pandas as pd
import numpy as np
import maths
from utils import *



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
        self.strand = self.t.transcript.Strand
        # extract alignment blocks. Those contain absolute coords and are always ordered by genomic coords 
        # (e.g., [(16101295, 16101359), (16101727, 16101936), (16102490, 16102599), (16102692, 16102829), (16103168, 16103334), (16103510, 16103684), (16104365, 16104662)])
        self.aln_blocks=[]
        bstart=self.t.transcript.Start
        for idx, splicing in list(enumerate(self.splicing_status)):
            if splicing==1:
                bend=self.t.introns.iloc[idx].Start-1 # last exonic 1-based pos
                self.aln_blocks+=[(bstart, bend)]
                bstart=self.t.introns.iloc[idx].End+1 # 1st exonic 1-based pos
        bend=self.t.transcript.End
        self.aln_blocks+=[(bstart, bend)]
        #print("alignment blocks: %s" % (self.aln_blocks))
        
    # convert isoform-relative coordinates to genomic coordinates.
    # E.g., a rel_pos==0 will return the 1st position of the transcript in genomic coordinates (1-based)
    def rel2abs_pos(self, rel_pos):
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
            print("WARN: Invalid relative position %i / %i in isoform %s: %i" % (rel_pos, abs_pos, self, off))
        ret = abs_pos+off if self.strand == '+' else abs_pos-off
        return (ret, block_id)
    
    def __repr__(self):
        return self.__str__()  
    
    def __str__(self):
        ret = ("%s_%s @ %s, [%s], [%s]" % (self.t.tid, self.id, self.t.region, ",".join(str(x) for x in self.fractions), ",".join(str(x) for x in self.splicing_status)) )
        return (ret)   

class Transcript():
    def __init__(self, tid, df, genome, conditions, max_ilen=None ):
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
                    print("Skipping transcript %s due to too long intron (%i vs %i). Skipping transcript..." % ( self.tid, ilen, max_ilen))
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
