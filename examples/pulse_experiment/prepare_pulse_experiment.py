import  os, sys
import pandas as pd
import tqdm
import pysam
import hashlib
import pybedtools
import json
from collections import OrderedDict
from argparse import ArgumentParser, RawDescriptionHelpFormatter

def sort_bgzip_and_tabix(in_file, out_file=None, create_index=True, del_uncompressed=True, seq_col=0, start_col=1, end_col=1, line_skip=0, zerobased=False):
    """ Will BGZIP the passed file and creates a tabix index with the given params if create_index is True"""
    if out_file is None:
        out_file=in_file+'.gz'
    pre,post=os.path.splitext(in_file)
    sorted_in_file=pre+'.sorted'+post
    pybedtools.BedTool(in_file).sort().saveas(sorted_in_file)    
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
            
def build_transcriptome(config, tab, out_dir):
    """ Creates a transcriptome file   
    """
    transcriptome_name=config.get('transcriptome_name', 'test')
    tids=tab['mode'].keys()
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
    out_file_bed=out_dir+'/'+transcriptome_name+'.bed'
    tid2start={}
    # process data
    with open(out_file_fasta, 'w') as out_fasta:
        with open(out_file_chrsize, 'w') as out_chrsize:        
            with open(out_file_dict, 'w') as out_dict:        
                with open(out_file_gff3, 'w') as out_gff3:
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

    # compress + index output files            
    sort_bgzip_and_tabix(out_file_gff3, seq_col=0, start_col=3, end_col=4, line_skip=0, zerobased=False)
    sort_bgzip_and_tabix(out_file_bed, seq_col=0, start_col=1, end_col=2, line_skip=0, zerobased=True)
    pysam.faidx(out_file_fasta)# @UndefinedVariable
    
    print("Created resources:")
    print("FASTA file + idx:\t"+out_file_fasta)
    print("CHROMSIZE file:\t"+out_file_chrsize)
    print("DICT file:\t"+out_file_dict)
    print("GFF file + idx:\t"+out_file_gff3)
    print("BED file + idx:\t"+out_file_bed)
                    

if __name__ == '__main__':
    
    parser = {}
    
    parser["prepare_pulse_experiment"] = ArgumentParser(description=usage, formatter_class=RawDescriptionHelpFormatter)
    parser["prepare_pulse_experiment"].add_argument("-c", "--config", type=str, required=True, dest="config_file", metavar="config_file", help="JSON config file")
    parser["prepare_pulse_experiment"].add_argument("-o", "--outdir", type=str, required=False, dest="outdir", metavar="outdir", help="output directory (default is current dir)")
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
        
    if mod == "prepare_pulse_experiment":
        # load and check onfig
        config=json.load(open(args.config_file), object_pairs_hook=OrderedDict)
    
        # read a list of tids and the respective mode
        tab=pd.read_csv(config['isoform_config'],delimiter='\t',encoding='utf-8').set_index('transcript_id').to_dict()
        print(tab)
        # build a reference transcriptome
        build_transcriptome(config, tab, out_dir+'ref/')
        
        # create timepoint configs
        tp_dir=out_dir+'/conf/'
        if not os.path.exists(tp_dir):
            os.makedirs(tp_dir)
           
        with open(out_dir+'run_splice_sim.sh', 'w') as script_out:
            print("""
#!/usr/bin/bash
ml purge &> /dev/null
ml load java/1.8.0_212
ml load nextflow/21.04.1

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
            """, file=script_out)
            for cr_mode in ['cr_const', 'cr_inc', 'cr_dec']:
                if cr_mode=='cr_const':
                    conversion_rates = [0.035]*config['timepoints']
                elif cr_mode=='cr_inc':
                    conversion_rates=[0.02, 0.025, 0.03, 0.035, 0.04, 0.045, 0.05]
                elif cr_mode=='cr_dec':
                    conversion_rates=list(reversed([0.02, 0.025, 0.03, 0.035, 0.04, 0.045, 0.05]))
                for tp in range(1, config['timepoints']+1):
                    tp_config={
                        "dataset_name": config['transcriptome_name']+'_tp%i'%(tp),
                        "splice_sim_cmd": config['splice_sim_cmd'],
                        "genome_fa": out_dir+'ref/%s.fa' % (config['transcriptome_name']),
                        "genome_chromosome_sizes": out_dir+'ref/%s.fa.chrom.sizes' % (config['transcriptome_name']),
                        "gene_gff": out_dir+'ref/%s.gff.gz' % (config['transcriptome_name']),
                        "create_bams": False,
                        "transcript_ids": "conf/isoform_data.tsv",
                        "isoform_mode": "from_file",
                        "condition": {
                            "ref": "T",
                            "alt": "C",
                            "conversion_rates": [ conversion_rates[tp-1] ],
                            "base_coverage": 10
                        },
                        "max_ilen": 100000,
                        "min_abundance": 1,
                        "random_seed": 1234,
                        "readlen": 100,
                        "create_tdf": False,
                        "write_reads": False
                        }
                    with open(tp_dir+'pulse_experiment.'+cr_mode+'.'+str(tp)+'.config.json', 'w') as out:
                        json.dump(tp_config, out, indent=4)
            
                    # write config per timepoint
                    with open(tp_dir+'pulse_experiment.'+cr_mode+'.'+str(tp)+'.isoform_data.tsv', 'w') as out:
                        print('transcript_id\tabundance\tfrac_mature\tfrac_old_mature', file=out)
                        for tid in tab['mode'].keys():
                            if tab['mode'][tid]==1: # constant (steady-state) abundance; decreasing fraction of unlabeled rnas
                                abundance=1
                                frac_mat=1
                                frac_old_mat=(config['timepoints']-tp)/(config['timepoints']-1)
                            elif tab['mode'][tid]==2: # transcription starts at half tp; no unlabeled rnas
                                abundance=0 if tp <= config['timepoints']/2 else 1
                                frac_mat=1
                                frac_old_mat=0
                            print('\t'.join([str(x) for x in [tid, abundance, frac_mat, frac_old_mat]]), file=out)
        
                    # write run command
                    print("nextflow -log logs/nextflow.log run ../nf1/splice_sim/splice_sim.nf -params-file %s -resume -profile $profile" % (tp_dir+'pulse_experiment.'+cr_mode+'.'+str(tp)+'.config.json'), file=script_out)    
    