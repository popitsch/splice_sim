"""
    @author: niko.popitsch@univie.ac.at

    Transcriptional model for splice_sim, containing transcripts, isoforms and conditions.
    Usage:
    1) build a model from a config file:
        m = Model.build_model( config_file )
    2) save model to disk:
        m.save(f_model)
    3) load model from disk:
        m =Model.load_from_file(f_model)

"""
import json
import logging
import math
import pickle
import sys
from collections import OrderedDict
from pathlib import Path

import pandas as pd
import numpy as np
import tqdm
from intervaltree import IntervalTree, Interval

from splice_sim.genomic_iterators import TabixIterator
from splice_sim.utils import *

# Necessary for including python modules from a parent directory
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def calculate_transcript_data(config, tid_meta, tid_file):
    """ Calculate a data config file """
    mode = config['isoform_mode'] if 'isoform_mode' in config else '1:1'
    logging.info("Calculating transcript configuration with mode %s" % mode)
    tdata = OrderedDict()
    # ------------------------------------------------------------------------------------
    #    1:1 mat/pre mix with fixed abundance for all transcripts.
    #    input: transcript_ids table with column transcript_id
    #    config: 
    # ------------------------------------------------------------------------------------
    if mode == '1:1':
        assert 'transcript_ids' in config, "For 1:1 mode, a list of transcript ids must be provided in a file referenced by the 'transcript_ids' config property"
        # generate data file
        base_abundance = config["base_abundance"] if 'base_abundance' in config else 10
        frac_old_mature = min(1.0, config[
            "frac_old_mature"]) if 'frac_old_mature' in config else 0.0  # fraction of 'old' (unlabeled), mature isoform
        frac_old_pre = min(1.0, config[
            "frac_old_pre"]) if 'frac_old_pre' in config else 0.0  # fraction of 'old' (unlabeled), premature isoform
        for tid, pinfo in tid_meta.items():
            rnk = pinfo['rnk'] - 1  # number of  introns
            gene_name = pinfo['gene_name']
            isoform_data = OrderedDict()
            if frac_old_mature > 0:  # add 'old' rna isoform
                isoform_data['old'] = OrderedDict()
                isoform_data['old']['splicing_status'] = [1] * rnk
                isoform_data['old']['fraction'] = frac_old_mature
                isoform_data['old']['is_labeled'] = 0
            if frac_old_pre > 0:  # add 'old' rna isoform
                isoform_data['oldpre'] = OrderedDict()
                isoform_data['oldpre']['splicing_status'] = [0] * rnk
                isoform_data['oldpre']['fraction'] = frac_old_pre
                isoform_data['oldpre']['is_labeled'] = 0
            if frac_old_mature < 1.0:  # premature isoform
                isoform_data['pre'] = OrderedDict()
                if rnk:  # at least one intron
                    isoform_data['pre']['splicing_status'] = [0] * rnk
                isoform_data['pre']['fraction'] = round((1 - (frac_old_mature + frac_old_pre)) * 0.5, 3)
                isoform_data['pre']['is_labeled'] = 1
            # output data
            tdata[tid] = OrderedDict()
            tdata[tid]["gene_name"] = gene_name
            tdata[tid]["abundance"] = base_abundance
            tdata[tid]["isoforms"] = isoform_data
    elif mode == 'from_file':
        assert 'transcript_ids' in config, "For 'from_file' mode, a transcript metadata file with headers ['transcript_id','abundance','frac_mature','frac_old_mature'] must be provided in a file referenced by the 'transcript_ids' config property"
        tid_table_data = pd.read_csv(tid_file, delimiter='\t', encoding='utf-8').set_index('transcript_id').to_dict()
        assert 'abundance' in tid_table_data.keys(), "For 'from_file' mode, a transcript metadata file with headers ['transcript_id','abundance','frac_mature','frac_old_mature'] must be provided in a file referenced by the 'transcript_ids' config property"
        assert 'frac_old_mature' in tid_table_data.keys(), "For 'from_file' mode, a transcript metadata file with headers ['transcript_id','abundance','frac_mature','frac_old_mature'] must be provided in a file referenced by the 'transcript_ids' config property"
        base_abundance = config["base_abundance"] if 'base_abundance' in config else 10
        for tid, pinfo in tid_meta.items():
            rnk = pinfo['rnk'] - 1  # number of  introns
            gene_name = pinfo['gene_name']
            isoform_data = OrderedDict()
            frac_mature = tid_table_data['frac_mature'][tid]
            if frac_mature == 'random':
                frac_mature = np.random.uniform(0, 1)
            frac_old_mature = tid_table_data['frac_old_mature'][tid]
            if frac_old_mature == 'random':
                frac_old_mature = np.random.uniform(0, 1)
            frac_old_pre = tid_table_data['frac_old_pre'][tid]
            if frac_old_pre == 'random':
                frac_old_pre = np.random.uniform(0, 1)
            abundance = tid_table_data['abundance'][tid]
            if abundance == 'random':
                abundance = np.random.uniform(0, 1)
            # create per-isoform data
            isoform_data = OrderedDict()
            if frac_old_mature > 0:  # add 'old' rna isoform
                isoform_data['old'] = OrderedDict()
                isoform_data['old']['splicing_status'] = [1] * rnk
                isoform_data['old']['fraction'] = round(frac_old_mature, 3)
                isoform_data['old']['is_labeled'] = 0
            if frac_old_pre > 0:  # add 'old' rna isoform
                isoform_data['oldpre'] = OrderedDict()
                isoform_data['oldpre']['splicing_status'] = [0] * rnk
                isoform_data['oldpre']['fraction'] = round(frac_old_pre, 3)
                isoform_data['oldpre']['is_labeled'] = 0
            if frac_old_mature < 1.0:
                isoform_data['pre'] = OrderedDict()
                if rnk:  # at least one intron
                    isoform_data['pre']['splicing_status'] = [0] * rnk
                isoform_data['pre']['fraction'] = round((1 - (frac_old_mature + frac_old_pre)) * (1 - frac_mature), 3)
                isoform_data['pre']['is_labeled'] = 1
            # output data
            tdata[tid] = OrderedDict()
            tdata[tid]["gene_name"] = gene_name
            tdata[tid]["abundance"] = base_abundance * abundance
            tdata[tid]["isoforms"] = isoform_data
    else:
        logging.error("Unknown mode %s " % (mode))
        sys.exit(1)
        # Write JSON file
    out_file = config['transcript_data']
    with open(out_file, 'w') as out:
        json.dump(tdata, out, indent=2)
    logging.info('Done. Results written to %s' % (out_file))
    return tdata


class Condition():
    def __init__(self, ref, alt, conversion_rates, base_coverage):
        self.ref = ref
        self.alt = alt
        self.conversion_rates = conversion_rates
        self.base_coverage = base_coverage

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        ret = ("condition %s>%s [cr=%f, cv=%f]" % (
        self.name, self.ref, self.alt, ','.join([str(x) for x in self.conversion_rates]), self.base_coverage))
        return (ret)


class Isoform():
    def __init__(self, id, fraction, splicing_status, is_labeled, transcript):
        self.id = id
        self.fraction = fraction
        self.splicing_status = splicing_status
        self.transcript = transcript
        self.aln_block_len = 0
        self.is_labeled = is_labeled
        if self.transcript is not None:
            self.strand = self.transcript.strand
            # extract alignment blocks. Those contain absolute coords and are always ordered by genomic coords 
            # (e.g., [(16101295, 16101359), (16101727, 16101936), (16102490, 16102599), (16102692, 16102829), (16103168, 16103334), (16103510, 16103684), (16104365, 16104662)])
            self.aln_blocks = []
            bstart = self.transcript.start
            for idx, splicing in list(enumerate(self.splicing_status)):
                if splicing == 1:
                    bend = self.transcript.introns[idx]['start'] - 1  # last exonic 1-based pos
                    self.aln_blocks += [(bstart, bend)]
                    self.aln_block_len += (bend - bstart + 1)
                    bstart = self.transcript.introns[idx]['end'] + 1  # 1st exonic 1-based pos
            bend = self.transcript.end
            self.aln_blocks += [(bstart, bend)]
            self.aln_block_len += (bend - bstart + 1)
        # print("iso %s, alignment blocks: %s, len: %i" % (self.id, self.aln_blocks, self.aln_block_len))

    def len(self):
        return self.aln_block_len

    def block_contains(self, bstart, bend, abs_pos):
        return abs_pos >= bstart and abs_pos <= bend

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        ret = ("%s_%s @ %s, [%f], [%s]" % (self.transcript.tid, self.id, self.transcript.region, self.fraction,
                                           ",".join(str(x) for x in self.splicing_status)))
        return (ret)


class Transcript():
    def __init__(self, config, tid, tid_meta, tid_exon, genome, condition, max_ilen=None):
        self.config = config
        self.tid = tid
        self.is_valid = True
        self.cond = condition
        self.transcript_data = tid_meta
        self.transcript_id = self.transcript_data['transcript_id']
        self.chromosome = self.transcript_data['chromosome']
        self.start = self.transcript_data['start']
        self.end = self.transcript_data['end']
        self.strand = self.transcript_data['strand']
        self.region = "%s:%i-%i" % (self.chromosome, self.start, self.end)
        self.len = self.end - self.start + 1
        self.transcript_seq = genome.fetch(region=self.region)
        self.exons = tid_exon
        self.introns = []
        last_end = -1
        for index, exon in enumerate(tid_exon):
            if last_end > 0:  # skip 1st
                ilen = (exon['start'] - 1) - (last_end + 1) + 1
                intron = {'transcript_id': self.transcript_id,
                          'Feature': 'intron',
                          'chromosome': self.chromosome,
                          'strand': self.strand,
                          'start': last_end + 1,
                          'end': exon['start'] - 1,
                          'exon_number': int(exon['exon_number']),
                          'len': ilen
                          }
                if (max_ilen is not None) and (ilen > max_ilen):
                    logging.warning(
                        "Skipping transcript %s due to too long intron (%i vs %i). Skipping transcript..." % (
                        self.tid, ilen, max_ilen))
                    self.is_valid = False
                    return
                self.introns.append(intron)
            last_end = exon['end']
        # set abundance and isoforms
        self.abundance = math.ceil(self.config['transcripts'][tid]['abundance'])  # abundance is float
        self.isoforms = OrderedDict()
        total_frac = 0
        for iso in self.config['transcripts'][tid]['isoforms'].keys():
            frac = self.config['transcripts'][tid]['isoforms'][iso]['fraction']
            splicing_status = self.config['transcripts'][tid]['isoforms'][iso][
                'splicing_status'] if 'splicing_status' in self.config['transcripts'][tid]['isoforms'][iso] else []
            is_labeled = self.config['transcripts'][tid]['isoforms'][iso]['is_labeled'] if 'is_labeled' in \
                                                                                           self.config['transcripts'][
                                                                                               tid]['isoforms'][
                                                                                               iso] else 1
            if self.strand == "-":
                splicing_status = list(reversed(splicing_status))  # so 1st entry in config refers to 1st intron.
            if len(splicing_status) != len(self.introns):
                logging.error(
                    "Misconfigured splicing status for some isoform of transcript %s (%i vs %i). Skipping transcript..." % (
                    self.tid, len(splicing_status), len(self.introns)))
                self.is_valid = False
                return
            self.isoforms[iso] = Isoform(iso, frac, splicing_status, is_labeled, self)
            total_frac += frac
        # assert with rounding error
        assert total_frac <= 1.00001, "Total fraction of isoforms > 1 (%f) for %s" % (total_frac, self.tid)
        # add mature, labeled isoform if not configured
        if 'mat' not in self.isoforms:
            frac = (1.0 - total_frac)
            iso = 'mat'
            splicing_status = [1] * len(self.introns)
            self.isoforms[iso] = Isoform(iso, frac, splicing_status, 1, self)  # add labeled, mature isoform
        # print('added mature isoform for %s: %s' % (tid, self.isoforms[iso]) )

    def get_dna_seq(self, splicing_status):
        seq = self.transcript_seq
        for idx, splicing in reversed(list(enumerate(splicing_status))):
            if splicing == 1:
                pos1 = self.introns[idx]['start'] - self.start
                pos2 = self.introns[idx]['end'] - self.start + 1
                seq = seq[:pos1] + seq[pos2:]
        return (seq)

    def get_rna_seq(self, iso):
        seq = self.get_dna_seq(self.isoforms[iso].splicing_status)
        if self.strand == '-':
            seq = reverse_complement(seq)
        return (seq)

    def get_sequences(self):
        ret = OrderedDict()
        total_count = 0
        last_expressed_iso = 'mat'
        for id in self.isoforms.keys():
            frac = self.isoforms[id].fraction
            if frac > 0:
                last_expressed_iso = id
            count = int(self.abundance * frac)
            total_count += count
            ret[id] = [
                self.get_rna_seq(id),
                count]
        toadd = self.abundance - total_count
        # assert toadd==0 | toadd==1 , "rounding error"
        if toadd > 0:  # add 1 to last isoform
            logging.debug("corrected %s counts by %i" % (last_expressed_iso, toadd))
            ret[last_expressed_iso] = [ret[last_expressed_iso][0], ret[last_expressed_iso][1] + toadd]
        return (ret)

    def to_bed(self):
        return ("%s\t%i\t%i\t%s\t%i\t%s" % (
        self.chromosome, self.start, self.end, self.transcript_data['ID'], 1000, self.strand))

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        ret = ("%s" % (self.transcript_id))
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

    def __init__(self, config, base_dir):
        self.config = config
        self.threads = config["threads"] if "threads" in config else 1
        # init condition
        self.condition = Condition(
            config["condition"]["ref"],
            config["condition"]["alt"],
            config["condition"]["conversion_rates"],
            config["condition"]["base_coverage"])
        # load external list of considered transcripts
        if "transcript_ids" in config:
            tfile = config["transcript_ids"]
            if not os.path.isabs(tfile):
                tfile = base_dir + '/' + tfile
            self.tids = set(pd.read_csv(tfile, delimiter='\t', encoding='utf-8')['transcript_id'])
            logging.info("read %i unique tids" % len(self.tids))
        else:
            self.tids = None
            logging.info("No tid list configured, will use all transcripts from the GFF3")
        # load + filter gene gff
        logging.info("Loading gene GFF")
        tid_meta = {}  # tid 2 metadata
        tid_exon = {}
        f = pysam.TabixFile(self.config["gene_gff"], mode="r")
        for row in tqdm.tqdm(f.fetch(parser=pysam.asTuple()), "Loading transcript data"):
            reference, source, ftype, fstart, fend, score, strand, phase, info = row
            if ftype == 'transcript' or ftype == 'exon':
                pinfo = parse_info(info)
                tid = pinfo['transcript_id'] if 'transcript_id' in pinfo else None
                if self.tids is None or tid in self.tids:
                    pinfo['chromosome'] = reference
                    pinfo['start'] = int(fstart)
                    pinfo['end'] = int(fend)
                    pinfo['strand'] = strand
                    if ftype == 'transcript':
                        tid_meta[pinfo['transcript_id']] = pinfo
                    elif ftype == 'exon':
                        if tid not in tid_exon:
                            tid_exon[tid] = []
                        tid_exon[tid].append(pinfo)
        for tid, exons in tid_exon.items():
            tid_meta[tid]['rnk'] = len(exons)
            # build transcript data file if none existing
        if "transcripts" not in config:
            logging.info("Creating external transcript data file %s" % tfile)
            tdata = calculate_transcript_data(config, tid_meta, tfile)
            config["transcripts"] = tdata
        # load genome
        logging.info("Loading genome")
        genome = pysam.FastaFile(config["genome_fa"])
        self.chrom_sizes = config["chrom_sizes"] if 'chrom_sizes' in config else config["genome_fa"] + ".chrom.sizes"
        # instantiate transcripts
        self.transcripts = OrderedDict()
        self.max_ilen = config["max_ilen"] if 'max_ilen' in config else None
        self.readlen = config["readlen"] if 'readlen' in config else None
        for tid in tid_meta.keys():
            t = Transcript(config, tid, tid_meta[tid], tid_exon[tid], genome, self.condition, max_ilen=self.max_ilen)
            if t.is_valid:
                self.transcripts[tid] = t
        self.snps = None
        self.snp_pos_per_chr = None
        if 'snp_file' in config:
            self.snp_pos_per_chr = {}
            self.snps = {}
            self.snp_file = config['snp_file']
            ti = TabixIterator(self.snp_file)
            for (CHROM, POS, ID, REF, ALT, _, _, INFO) in ti:
                if CHROM.startswith("#"):
                    continue
                if CHROM not in self.snp_pos_per_chr:
                    self.snp_pos_per_chr[CHROM] = set()
                self.snp_pos_per_chr[CHROM].add(int(POS))
                info = self.parse_info(INFO)
                prob = float(info['prob']) if 'prob' in info else 1  # parse from info
                strand_specific = info['strand'] if 'strand' in info else None
                enable_nc = info['enable_nc'] == 'yes' if 'enable_nc' in info else True
                self.snps[CHROM + ':' + POS] = (ALT, prob, strand_specific, enable_nc)
            logging.info("Loaded %i SNPs" % (sum([len(x) for x in self.snp_pos_per_chr.values()])))
        logging.info("Instantiated %i transcripts" % len(self.transcripts))

    def parse_info(self, info):
        """ parse GFF3 info section """
        return {k: v for k, v in [a.split('=') for a in info.split(';') if '=' in a]}

    def write_gff(self, outdir):
        """ Write a filtered GFF file containing all kept/simulated transcripts """
        logging.info("Writing filtered gene GFF")
        out_file = outdir + "gene_anno.gff3"
        f = pysam.TabixFile(self.config["gene_gff"], mode="r")
        with open(out_file, 'w') as out:
            for row in f.fetch(parser=pysam.asTuple()):
                reference, source, ftype, fstart, fend, score, strand, phase, info = row
                pinfo = parse_info(info)
                tid = pinfo['transcript_id'] if 'transcript_id' in pinfo else None
                if tid in self.transcripts.keys():
                    print('\t'.join([str(x) for x in row]), file=out)
            bgzip_and_tabix(out_file, seq_col=0, start_col=3, end_col=4, line_skip=0, zerobased=False)
        out_file = out_file + ".gz"
        return out_file

    def write_transcript_md(self, outdir):
        """ Write a filtered TSV file containing md about all kept/simulated transcripts """
        logging.info("Writing filtered gene TSV")
        out_file = outdir + "gene_anno.tsv"
        headers = None
        if not os.path.exists(out_file + ".gz"):
            with open(out_file, 'w') as out:
                for t in self.transcripts.values():
                    if headers is None:
                        headers = [str(h) for h in t.transcript_data if
                                   h not in ['Source', 'Feature', 'Score', 'Frame', 'ID', 'havana_gene', 'Parent',
                                             'exon_number', 'exon_id']]
                        print('\t'.join(headers), file=out)
                    print('\t'.join([str(t.transcript_data[h]) if h in t.transcript_data else 'NA' for h in headers]),
                          file=out)
            bgzip_and_tabix(out_file, create_index=False)
        out_file = out_file + ".gz"
        return out_file

    def write_isoform_truth(self, outdir):
        """ Write a TSV file with isoform truth data"""
        logging.info("Writing isoform truth data")
        out_file = outdir + self.config["dataset_name"] + ".isoform_truth.tsv"
        if not os.path.exists(out_file + ".gz"):
            with open(out_file, 'w') as out:
                print('\t'.join([str(x) for x in
                                 ['tid', 'iso', 'is_labeled', 'true_abundance', 'true_fraction', 'splicing_status',
                                  'feature_length']]), file=out)
                for t in self.transcripts.values():
                    for iso in t.isoforms.values():
                        splice_status = ','.join([str(y) for y in iso.splicing_status]) if len(
                            iso.splicing_status) > 0 else 'NA'
                        print('\t'.join([str(x) for x in
                                         [t.tid, iso.id, iso.is_labeled, t.abundance, iso.fraction, splice_status,
                                          iso.len()]]), file=out)
            bgzip_and_tabix(out_file, create_index=False)
        out_file = out_file + ".gz"
        return out_file

    def write_isoform_sequences(self, outdir):
        fout = outdir + self.config['dataset_name'] + ".fa"
        buf = []
        with open(fout, 'w') as out:
            for t in self.transcripts.values():
                sequences = t.get_sequences()
                for iso in t.isoforms.keys():
                    for i in range(sequences[iso][1]):
                        buf += [">%s_%s_%s_%i\n" % (t.tid, t.strand, iso, i)]
                        buf += [pad_n(sequences[iso][0], self.config['readlen']) + "\n"]
                        if len(buf) > 10000:
                            out.writelines(buf)
                            buf = []
            if len(buf) > 0:
                out.writelines(buf)
        return fout

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
        ret = ("Model [t: %i]" % (len(self.transcripts)))
        return (ret)

    def build_transcript_iv(self):
        # build interval tree of transcripts
        ivs = {}
        for tid, t in self.transcripts.items():
            if t.chromosome not in ivs:
                ivs[t.chromosome] = IntervalTree()
            ivs[t.chromosome].add(Interval(t.start, t.end, t))
        return ivs

    @classmethod
    def load_from_file(cls, model_file):
        return pickle.load(open(model_file, "rb"))

    @classmethod
    def build_model(cls, config_file, out_dir):
        """ build transcriptional model and save """
        # load config 
        config = json.load(open(config_file), object_pairs_hook=OrderedDict)
        if config_file.startswith('/Volumes'):  # localize file paths (on MacOs)
            config = localize_config(config)
        base_dir = str(Path(config_file).parent.absolute()) + '/'
        # ensure out_dir
        if not os.path.exists(out_dir):
            os.mkdir(out_dir)
        # build model
        m = cls(config, base_dir)
        f_model = out_dir + '/' + config["dataset_name"] + ".model"
        print("Model saved to " + f_model)
        m.save(f_model)
        # write final isoform model truth file
        m.write_isoform_truth(out_dir)
        # write considered transcripts to GFF
        m.write_gff(out_dir)
        # write considered transcripts to tsv
        m.write_transcript_md(out_dir)
        # write isoform sequences
        m.write_isoform_sequences(out_dir)
        return m, f_model
