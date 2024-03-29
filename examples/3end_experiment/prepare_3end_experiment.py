import  os, sys
import pandas as pd
import tqdm
import math
import pysam
import hashlib
import pybedtools
import json
from collections import OrderedDict
from argparse import ArgumentParser, RawDescriptionHelpFormatter

# idea
# 0) subset transcripts to the last 250nt (3' ends)
# 1) simulate reads only from those 3'ends to mimic 3'end sequencing
# 2) do both transcriptome and genome simulation

def sort_bgzip_and_tabix(in_file, out_file=None, create_index=True, del_uncompressed=True, seq_col=0, start_col=1, end_col=1, line_skip=0, zerobased=False):
    """ Will BGZIP the passed file and creates a tabix index with the given params if create_index is True"""
    if out_file is None:
        out_file=in_file+'.gz'
    pre,post=os.path.splitext(in_file)
    sorted_in_file=pre+'.sorted'+post
    #pybedtools.BedTool(in_file).sort().saveas(sorted_in_file)
    os.system("gff3sort.pl --precise " + in_file + " > " + sorted_in_file)
    pysam.tabix_compress(sorted_in_file, out_file, force=True) # @UndefinedVariable
    if create_index:
        pysam.tabix_index(out_file, force=True, seq_col=seq_col, start_col=start_col, end_col=end_col, meta_char='#', line_skip=line_skip, zerobased=zerobased) # @UndefinedVariable
    if del_uncompressed:
        os.remove(in_file)
        
def format_fasta(string, ncol=80):
    return '\n'.join(string[i:i+ncol] for i in range(0, len(string), ncol))

def parse_info(info):
    """ parse GFF3 info section """
    return {k:v for k,v in [a.split('=') for a in info.split(';') if '=' in a]}

def chop3End(tid ,transcriptMeta, transcriptBuffer, nt):
    consumedWindow = nt
    ends = list()
    exons = list(transcriptBuffer[tid])
    txstart = transcriptMeta[tid]['start']
    txend = transcriptMeta[tid]['end']
    if transcriptMeta[tid]['strand'] == "+":
        exons.reverse()
    while consumedWindow > 0 and exons:
        reference, source, ftype, fstart, fend, score, strand, phase, info = exons.pop(0)
        flength = int(fend) - int(fstart)
        if flength - consumedWindow > 0 :
            if transcriptMeta[tid]['strand'] == "+":
                fstart = int(fend) - consumedWindow + 1
                txstart = str(fstart)
            else :
                fend = int(fstart) + consumedWindow - 1
                txend = str(fend)
            consumedWindow = 0
        else :
            consumedWindow = consumedWindow - flength
        ends.append("\t".join([reference, source, ftype, str(fstart), str(fend), score, strand, phase, info]))

    if transcriptMeta[tid]['strand'] == "+":
        ends.reverse()

    return txstart, txend, ends

def extract_transcript_3ends(config, gff_file, tids, out_dir) :
    dataset_name = config.get('dataset_name', 'test')
    print('subsetting transcripts to 3\'ends for ', len(tids), ' tids.')
    out_file_gff3 = out_dir + '/' + dataset_name + '.3ends.gff3'
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    # Collect exons per transcript for 3'end subsetting
    transcriptBuffer = dict()
    transcriptMeta = dict()

    if not 'transcript_end_window' in config:
        config['transcript_end_window'] = 250

    # process data
    with open(out_file_gff3, 'w') as out_gff3:
        f = pysam.TabixFile(gff_file, mode="r")
        for row in tqdm.tqdm(f.fetch(parser=pysam.asTuple()), desc="Processing features"):
            reference, source, ftype, fstart, fend, score, strand, phase, info = row
            pinfo = parse_info(info)
            tid = pinfo['transcript_id'] if 'transcript_id' in pinfo else None
            if tid not in tids:
                continue
            if ftype != 'transcript' and ftype != 'exon':
                continue
            if ftype == 'transcript':
                transcriptMeta[tid] = dict()
                transcriptMeta[tid]['reference'] = reference
                transcriptMeta[tid]['source'] = source
                transcriptMeta[tid]['ftype'] = ftype
                transcriptMeta[tid]['start'] = fstart
                transcriptMeta[tid]['end'] = fend
                transcriptMeta[tid]['score'] = score
                transcriptMeta[tid]['strand'] = strand
                transcriptMeta[tid]['phase'] = phase
                transcriptMeta[tid]['info'] = info
                transcriptBuffer[tid] = list()
            if ftype == "exon":
                transcriptBuffer[tid].append(row)
                if fend == transcriptMeta[tid]['end'] :
                    txstart, txend, ends = chop3End(tid ,transcriptMeta, transcriptBuffer, config['transcript_end_window'])
                    print("\t".join([str(x) for x in [
                        transcriptMeta[tid]['reference'],
                        transcriptMeta[tid]['source'],
                        transcriptMeta[tid]['ftype'],
                        txstart,  # start
                        txend,  # end
                        transcriptMeta[tid]['score'],
                        transcriptMeta[tid]['strand'],
                        transcriptMeta[tid]['phase'],
                        transcriptMeta[tid]['info']
                    ]]), file=out_gff3)
                    for end in ends:
                        print(end, file=out_gff3)
                    del transcriptMeta[tid]
                    del transcriptBuffer[tid]

    # compress + index output files
    sort_bgzip_and_tabix(out_file_gff3, seq_col=0, start_col=3, end_col=4, line_skip=0, zerobased=False)

    print("Created resources:")
    print("GFF file + idx:\t" + out_file_gff3)

    return out_file_gff3 + ".gz"

def build_transcriptome(config, tids, out_dir):
    """ Creates a transcriptome file   
    """
    transcriptome_name=config.get('transcriptome_name', 'test')
    print('extracting sequences for ', len(tids), ' tids.')
    genome = pysam.FastaFile(config['genome_fa'])# @UndefinedVariable
    #dict_chr2idx, dict_idx2chr, dict_chr2len = get_chrom_dicts(fasta_file)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    padding=config['padding'] if 'padding' in config else 100
    amp_extension=config['amp_extension'] if 'amp_extension' in config else 100 # extend amplicon seq updownstream
    out_file_fasta=out_dir+'/'+transcriptome_name+'.fa'
    out_file_chrsize=out_dir+'/'+transcriptome_name+'.fa.chrom.sizes'
    out_file_dict=out_dir+'/'+transcriptome_name+'.fa.dict'
    out_file_gff3=out_dir+'/'+transcriptome_name+'.gff3'
    out_file_gtf=out_dir+'/'+transcriptome_name+'.gtf'
    out_file_bed=out_dir+'/'+transcriptome_name+'.bed'
    tid2start={}
    # process data
    with open(out_file_fasta, 'w') as out_fasta:
        with open(out_file_chrsize, 'w') as out_chrsize:        
            with open(out_file_dict, 'w') as out_dict:        
                with open(out_file_gff3, 'w') as out_gff3:
                    with open(out_file_gtf, 'w') as out_gtf:
                        with open(out_file_bed, 'w') as out_bed:
                            f = pysam.TabixFile(config["gene_gff"], mode="r")
                            for row in tqdm.tqdm(f.fetch(parser=pysam.asTuple()), desc="Processing features"):
                                reference, source,ftype,fstart,fend,score,strand,phase,info=row
                                pinfo=parse_info( info )
                                tid = pinfo['transcript_id'] if 'transcript_id' in pinfo else None
                                if tid not in tids:
                                    continue
                                if ftype != 'transcript':
                                    continue
                                fstart=int(fstart)-1
                                fend=int(fend)
                                ampstart=max(1,fstart-amp_extension)
                                ampend=fend+amp_extension
                                offset=padding+(fstart-ampstart)+1
                                seq='N'*padding+genome.fetch(reference, ampstart, ampend)+'N'*padding
                                new_chrom=tid
                                new_chrom_len=ampend-ampstart+1+2*padding # length of new chrom
                                tid2start[tid]=ampstart
                                # FASTA
                                print('>%s'% new_chrom, file=out_fasta)
                                print(format_fasta(seq), file=out_fasta)
                                # chrom.sizes file
                                print('%s\t%i'% (new_chrom,new_chrom_len), file=out_chrsize)
                                # DICT
                                print('@SQ\tSN:%s\tLN:%i\tM5:%s\tUR:file:%s'% (new_chrom,
                                                                               new_chrom_len,
                                                                               hashlib.md5(seq.encode('utf-8')).hexdigest(), # M5 MD5 checksum of the sequence in the uppercase, excluding spaces but including pads (as ‘*’s)
                                                                               os.path.abspath(out_file_fasta)), file=out_dict)

                                # BED
                                print("\t".join([str(x) for x in [
                                    new_chrom,
                                    offset-1, # start
                                    offset+fend-fstart-1, # end
                                    new_chrom,
                                    '.' if score is None else score,
                                    '.' if strand is None else strand
                                    ]]), file=out_bed)

                            f = pysam.TabixFile(config["gene_gff"], mode="r")
                            for row in tqdm.tqdm(f.fetch(parser=pysam.asTuple()), desc="Processing features"):
                                reference, source,ftype,fstart,fend,score,strand,phase,info=row
                                pinfo=parse_info( info )
                                tid = pinfo['transcript_id'] if 'transcript_id' in pinfo else None
                                if tid not in tids:
                                    continue
                                fstart=int(fstart)-1
                                fend=int(fend)
                                ampstart=tid2start[tid]
                                offset=padding+(fstart-ampstart)+1
                                new_chrom=tid
                                # GFF
                                print("\t".join([str(x) for x in [
                                    new_chrom,
                                    source,
                                    ftype,
                                    offset, # start
                                    offset+fend-fstart-1, # end
                                    score,
                                    strand,
                                    phase,
                                    info
                                    ]]), file=out_gff3)
                                # GTF
                                if ftype in ['transcript', 'exon']:
                                    gtf_info = 'transcript_id "%s"; gene_id "%s"; gene_name "%s";' % (pinfo['transcript_id'], pinfo['gene_id'], pinfo['gene_name'])
                                    print("\t".join([str(x) for x in [
                                        new_chrom,
                                        source,
                                        ftype,
                                        offset,  # start
                                        offset + fend - fstart - 1,  # end
                                        score,
                                        strand,
                                        phase,
                                        gtf_info
                                    ]]), file=out_gtf)

    # compress + index output files            
    sort_bgzip_and_tabix(out_file_gff3, seq_col=0, start_col=3, end_col=4, line_skip=0, zerobased=False)
    sort_bgzip_and_tabix(out_file_bed, seq_col=0, start_col=1, end_col=2, line_skip=0, zerobased=True)
    pysam.faidx(out_file_fasta)# @UndefinedVariable
    
    print("Created resources:")
    print("FASTA file + idx:\t"+out_file_fasta)
    print("CHROMSIZE file:\t"+out_file_chrsize)
    print("DICT file:\t"+out_file_dict)
    print("GFF file + idx:\t"+out_file_gff3)
    print("GTF file:\t"+out_file_gtf)
    print("BED file + idx:\t"+out_file_bed)
                    
usage = '''                           

  Copyright 2021 Niko Popitsch. All rights reserved.
  
  Licensed under the Apache License 2.0
  http://www.apache.org/licenses/LICENSE-2.0
  
  Distributed on an "AS IS" basis without warranties
  or conditions of any kind, either express or implied.

'''
if __name__ == '__main__':
    
    parser = {}
    
    parser["prepare_3end_experiment"] = ArgumentParser(description=usage, formatter_class=RawDescriptionHelpFormatter)
    parser["prepare_3end_experiment"].add_argument("-c", "--config", type=str, required=True, dest="config_file", metavar="config_file", help="JSON config file")
    parser["prepare_3end_experiment"].add_argument("-t", "--transcript", action='store_true', dest="transcript", help="Simulate transcript-wide read noise")
    parser["prepare_3end_experiment"].add_argument("-o", "--outdir", type=str, required=False, dest="outdir", metavar="outdir", help="output directory (default is current dir)")

    parser["prepare_3end_experiment_genome"] = ArgumentParser(description=usage, formatter_class=RawDescriptionHelpFormatter)
    parser["prepare_3end_experiment_genome"].add_argument("-c", "--config", type=str, required=True, dest="config_file", metavar="config_file", help="JSON config file")
    parser["prepare_3end_experiment_genome"].add_argument("-t", "--transcript", action='store_true', dest="transcript", help="Simulate transcript-wide read noise")

    parser["prepare_3end_experiment_genome"].add_argument("-o", "--outdir", type=str, required=False, dest="outdir", metavar="outdir", help="output directory (default is current dir)")
#============================================================================    
    if len(sys.argv) <= 1 or sys.argv[1] in ['-h', '--help']:
        print("usage: splice_sim.py [-h] " + ",".join(parser.keys()))
        sys.exit(1)
    mod = sys.argv[1]  
    if mod not in parser.keys():
        print("Invalid module '%s' selected. Please use one of %s" % (mod, ",".join(parser.keys())))
        sys.exit(1)

    args = parser[mod].parse_args(sys.argv[2:])
    #============================================================================

    # output dir (current dir if none provided)
    out_dir = os.path.abspath(args.outdir if args.outdir else os.getcwd())+'/'
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
        
    if mod == "prepare_3end_experiment":
        # load and check onfig
        config=json.load(open(args.config_file), object_pairs_hook=OrderedDict)
    
        # read a list of tids and metadata:
        # headers: transcript_id, k
        tab=pd.read_csv(config['isoform_config'],delimiter='\t',encoding='utf-8')
        print(tab)

        if args.transcript:

            # build a reference transcriptome
            build_transcriptome(config, tab['transcript_id'].tolist(), out_dir + 'ref/')

            transcriptome_gff = out_dir + 'ref/%s.gff3.gz' % (config['transcriptome_name'])

            # build a 3'end reference from transcriptome
            end_gff3 = extract_transcript_3ends(config, transcriptome_gff, tab['transcript_id'].tolist(), out_dir + 'ref/')

        else:
            # build a 3'end reference
            end_gff3 = extract_transcript_3ends(config, config["gene_gff"], tab['transcript_id'].tolist(), out_dir + 'ref/')

            config["gene_gff"] = end_gff3

            # build a reference transcriptome
            build_transcriptome(config, tab['transcript_id'].tolist(), out_dir + 'ref/')
        
        header="""
#!/usr/bin/bash
ml purge &> /dev/null
ml load java/1.8.0_212
ml load nextflow/21.04.1
set -e

export NXF_WORK="/scratch/nextflow/"$USER"/nxf_work"
export NXF_TEMP="/scratch/nextflow/"$USER"/nxf_temp"
export NXF_ANSI_LOG=false
export NXF_OPTS="-Xms50m -Xmx500m"
mkdir -p logs
profile="standard"
if [ -n "$1" ]; then
  profile=$1
fi
echo "Starting splice_sim pipeline with profile $profile"
            """
           
        with open(out_dir+'run_splice_sim.sh', 'w') as script_out:
            with open(out_dir+'run_splice_sim_eva.sh', 'w') as script_out2:
                print(header, file=script_out)
                print(header, file=script_out2)
                splice_sim_config={
                    "dataset_name": config['transcriptome_name'],
                    "splice_sim_cmd": config['splice_sim_cmd'],
                    "genome_fa": out_dir+'ref/%s.fa' % (config['transcriptome_name']),
                    "genome_chromosome_sizes": out_dir+'ref/%s.fa.chrom.sizes' % (config['transcriptome_name']),
                    "gene_gff": out_dir+'ref/%s.gff3.gz' % (config['transcriptome_name']),
                    "create_bams": True,
                    "transcript_ids": out_dir+'ref/%s.tids.tsv' % (config['transcriptome_name']),
                    "transcript_data": out_dir + 'ref/%s.data.config.json' % (config['transcriptome_name']),
                    "isoform_mode": "1:1",
                    "frac_old_mature": 0,
                    "genome_mappability": out_dir+'ref/%s.genmap.bedgraph.gz' % (config['transcriptome_name']),
                    "condition": {
                        "ref": "T",
                        "alt": "C",
                        "conversion_rates": config['conversion_rates'],
                        "base_coverage": 10
                    },
                    "mappers": {
                        "STAR": {
                            "star_cmd": "STAR",
                            "star_genome_idx": out_dir+'ref/star_2.7.1_index/', # NB must be built externally
                            "star_splice_gtf": out_dir+'ref/%s.gtf' % (config['transcriptome_name']),
                            },
                        "HISAT3N": {
                            "hisat3n_cmd": "singularity exec /groups/ameres/Niko/software/SIF/hisat-3n.sif /hisat-3n/hisat-3n",
                            "hisat3n_idx": out_dir+'ref/hisat2-3n_index/%s' % (config['transcriptome_name']), # NB must be built externally
                            "hisat3n_kss": out_dir+'ref/%s.gtf.hisat2_splice_sites.txt' % (config['transcriptome_name'])
                            }
                        },
                    "max_ilen": 100000,
                    "min_abundance": 1,
                    "random_seed": 1234,
                    "readlen": 100,
                    "create_tdf": True,
                    "write_reads": False,
                    "write_intron_bam": False,
                    "rseqc": False
                    }
                with open(out_dir+'splice_sim.config.json', 'w') as out:
                    json.dump(splice_sim_config, out, indent=4)

                if (args.transcript) :
                    splice_sim_eva_config = {
                        "dataset_name": config['transcriptome_name'],
                        "splice_sim_cmd": config['splice_sim_cmd'],
                        "genome_fa": out_dir + 'ref/%s.fa' % (config['transcriptome_name']),
                        "genome_chromosome_sizes": out_dir + 'ref/%s.fa.chrom.sizes' % (config['transcriptome_name']),
                        "gene_gff": end_gff3,
                        "create_bams": True,
                        "transcript_ids": out_dir + 'ref/%s.tids.tsv' % (config['transcriptome_name']),
                        "transcript_data": out_dir + 'ref/%s.data.config.json' % (config['transcriptome_name']),
                        "isoform_mode": "1:1",
                        "frac_old_mature": 0,
                        "genome_mappability": out_dir + 'ref/%s.genmap.bedgraph.gz' % (config['transcriptome_name']),
                        "condition": {
                            "ref": "T",
                            "alt": "C",
                            "conversion_rates": config['conversion_rates'],
                            "base_coverage": 10
                        },
                        "mappers": {
                            "STAR": {
                                "star_cmd": "STAR",
                                "star_genome_idx": out_dir + 'ref/star_2.7.1_index/',  # NB must be built externally
                                "star_splice_gtf": out_dir + 'ref/%s.gtf' % (config['transcriptome_name']),
                            },
                            "HISAT3N": {
                                "hisat3n_cmd": "singularity exec /groups/ameres/Niko/software/SIF/hisat-3n.sif /hisat-3n/hisat-3n",
                                "hisat3n_idx": out_dir + 'ref/hisat2-3n_index/%s' % (config['transcriptome_name']),
                                # NB must be built externally
                                "hisat3n_kss": out_dir + 'ref/%s.gtf.hisat2_splice_sites.txt' % (
                                config['transcriptome_name'])
                            }
                        },
                        "max_ilen": 100000,
                        "min_abundance": 1,
                        "random_seed": 1234,
                        "readlen": 100,
                        "create_tdf": True,
                        "write_reads": False,
                        "write_intron_bam": False,
                        "rseqc": False
                    }
                    with open(out_dir + 'splice_sim.config.eva.json', 'w') as out:
                        json.dump(splice_sim_eva_config, out, indent=4)

                # write config per timepoint
                with open(out_dir+'ref/%s.tids.tsv' % (config['transcriptome_name']), 'w') as out:
                    print('transcript_id', file=out)
                    for tid in tab['transcript_id'].tolist():
                        print(tid, file = out)

                # write run commands
                print("nextflow -log logs/nextflow.log run "+config["splice_sim_nf"]+" -params-file %s -resume -profile $profile" % (out_dir+'splice_sim.config.json'), file=script_out)
                if (args.transcript):
                    print("nextflow -log logs/nextflow.log run "+config["splice_sim_eva_nf"]+" -params-file %s -resume -profile $profile" % (out_dir+'splice_sim.config.eva.json'), file=script_out2)
                else :
                    print("nextflow -log logs/nextflow.log run " + config["splice_sim_eva_nf"] + " -params-file %s -resume -profile $profile" % (out_dir + 'splice_sim.config.json'), file=script_out2)
#============================================================================
    if mod == "prepare_3end_experiment_genome":
        # load and check config
        config=json.load(open(args.config_file), object_pairs_hook=OrderedDict)
    
        # read a list of tids and metadata:
        # headers: transcript_id, k
        tab=pd.read_csv(config['isoform_config'],delimiter='\t',encoding='utf-8')
        print(tab)
        # build a 3'end reference
        endsGff = extract_transcript_3ends(config, config["gene_gff"], tab['transcript_id'].tolist(), out_dir + 'ref/')
    
        header="""
#!/usr/bin/bash
ml purge &> /dev/null
ml load java/1.8.0_212
ml load nextflow/21.04.1
set -e

export NXF_WORK="/scratch/nextflow/"$USER"/nxf_work"
export NXF_TEMP="/scratch/nextflow/"$USER"/nxf_temp"
export NXF_ANSI_LOG=false
export NXF_OPTS="-Xms50m -Xmx500m"
mkdir -p logs
profile="standard"
if [ -n "$1" ]; then
  profile=$1
fi
echo "Starting splice_sim pipeline with profile $profile"
            """

        # create config
        splice_sim_config={
            "dataset_name": config['dataset_name'],
            "splice_sim_cmd": config['splice_sim_cmd'],
            "genome_fa": config['genome_fa'],
            "genome_chromosome_sizes": config['genome_chromosome_sizes'],
            "gene_gff": endsGff,
            "create_bams": True,
            "frac_old_mature": 0,
            "transcript_ids": config['isoform_config'],
            "transcript_data": out_dir + '/' + config['dataset_name'] + '.data.config.json',
            "isoform_mode": "1:1",
            "genome_mappability": config['genome_mappability'],
            "condition": {
                "ref": "T",
                "alt": "C",
                "conversion_rates": config['conversion_rates'],
                "base_coverage": 10
            },
            "mappers": {
                "STAR": {
                    "star_cmd": "STAR",
                    "star_genome_idx": "/groups/ameres/Niko/ref/genomes/mm10/indices/star_2.7.1",
                    "star_splice_gtf": "/groups/ameres/Niko/ref/genomes/mm10/annotation/GRCm38.p6.annotation.gtf"
                },
                "HISAT3N": {
                    "hisat3n_cmd": "singularity exec /groups/ameres/Niko/software/SIF/hisat-3n.sif /hisat-3n/hisat-3n",
                    "hisat3n_idx": "/groups/ameres/Niko/ref/genomes/mm10/indices/hisat2-3n/Mus_musculus.GRCm38.dna.primary_assembly",
                    "hisat3n_kss": "/groups/ameres/Niko/ref/genomes/mm10/annotation/GRCm38.p6.annotation.gtf.hisat2_splice_sites.txt"
                    }
                },
            "max_ilen": 100000,
            "min_abundance": 1,
            "random_seed": 1234,
            "readlen": 100,
            "create_tdf": True,
            "write_reads": False,
            "write_intron_bam": False,
            "rseqc": False
        }

        with open(out_dir + 'splice_sim.config.json', 'w') as out:
            json.dump(splice_sim_config, out, indent=4)

        if (args.transcript):
            # create config
            splice_sim_simulator_config = {
                "dataset_name": config['dataset_name'],
                "splice_sim_cmd": config['splice_sim_cmd'],
                "genome_fa": config['genome_fa'],
                "genome_chromosome_sizes": config['genome_chromosome_sizes'],
                "gene_gff": config["gene_gff"],
                "create_bams": True,
                "frac_old_mature": 0,
                "transcript_ids": config['isoform_config'],
                "transcript_data": out_dir + '/' + config['dataset_name'] + '.data.config.json',
                "isoform_mode": "1:1",
                "genome_mappability": config['genome_mappability'],
                "condition": {
                    "ref": "T",
                    "alt": "C",
                    "conversion_rates": config['conversion_rates'],
                    "base_coverage": 10
                },
                "mappers": {
                    "STAR": {
                        "star_cmd": "STAR",
                        "star_genome_idx": "/groups/ameres/Niko/ref/genomes/mm10/indices/star_2.7.1",
                        "star_splice_gtf": "/groups/ameres/Niko/ref/genomes/mm10/annotation/GRCm38.p6.annotation.gtf"
                    },
                    "HISAT3N": {
                        "hisat3n_cmd": "singularity exec /groups/ameres/Niko/software/SIF/hisat-3n.sif /hisat-3n/hisat-3n",
                        "hisat3n_idx": "/groups/ameres/Niko/ref/genomes/mm10/indices/hisat2-3n/Mus_musculus.GRCm38.dna.primary_assembly",
                        "hisat3n_kss": "/groups/ameres/Niko/ref/genomes/mm10/annotation/GRCm38.p6.annotation.gtf.hisat2_splice_sites.txt"
                    }
                },
                "max_ilen": 100000,
                "min_abundance": 1,
                "random_seed": 1234,
                "readlen": 100,
                "create_tdf": True,
                "write_reads": False,
                "write_intron_bam": False,
                "rseqc": False
            }

            with open(out_dir + 'splice_sim.config.simulator.json', 'w') as out:
                json.dump(splice_sim_simulator_config, out, indent=4)

        # write run commands
        with open(out_dir+'run_splice_sim.sh', 'w') as script_out:
            print(header, file=script_out)
            if args.transcript:
                print("nextflow -log logs/nextflow.log run " + config["splice_sim_nf"] + " -params-file %s -resume -profile $profile" % (out_dir + 'splice_sim.config.simulator.json'), file=script_out)
            else :
                print("nextflow -log logs/nextflow.log run "+config["splice_sim_nf"]+" -params-file %s -resume -profile $profile" % (out_dir + 'splice_sim.config.json'), file=script_out)
        with open(out_dir+'run_splice_sim_eva.sh', 'w') as script_out:
            print(header, file=script_out)
            print("nextflow -log logs/nextflow.log run "+config["splice_sim_eva_nf"]+" -params-file %s -resume -profile $profile" % (out_dir + 'splice_sim.config.json'), file=script_out)
        with open(out_dir+'run_splice_sim_count.sh', 'w') as script_out:
            print(header, file=script_out)
            print("nextflow -log logs/nextflow.log run "+config["splice_sim_count_nf"]+" -params-file %s -resume -profile $profile" % (out_dir + 'splice_sim.config.json'), file=script_out)