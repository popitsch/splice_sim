#!/usr/bin/env nextflow

Channel.fromFilePairs("${params.bam_dir}/*.{bam,bai}", flat:true) { file -> file.name.replaceAll(/.bam|.bai$/,'') }.set{ input_bams }

/*
 * count tc reads
 */
process highlight_tc_reads {
	tag "$name" 
    cpus 1
    label "long"
    time 8.h
    module 'python/3.7.2-gcccore-8.2.0:sambamba/0.6.6:htslib/1.10.2-gcccore-7.3.0'
    publishDir "counts", mode: 'copy'
    input:   
    	set name, file(bam), file(bai) from input_bams
    output: 
		file '*+tc.bam' into highlight_tc_bam
		file '*+tc.bam.bai' into highlight_tc_bai
		file '*_tc_only.bam' into highlight_tc_only_bam
		file '*_tc_only.bam.bai' into highlight_tc_only_bai
		file '*tc_snps.vcf' into highlight_tc_bam_vcf
		file '*tc_hist.tsv' into highlight_tc_bam_hist
		file '*tcreadpos_hist.tsv' into highlight_tc_bam_readposhist
    	file("stats/*") into highlight_tc_bam_stats
    script:
	    """
	    # highlight  TC conversions
	    python ${params.splice_sim_count_cmd} annotate \
    		 --bam ${bam} \
			 --threads ${task.cpus} \
			 --config_string '{ "region_width": 50000, "min_base_quality": 0,  "bam_tags": { "num_tc": "yc","num_tt": "yt","density_tc": "yd","pos_tc": "yp","is_tc_duplicate": "yx","is_high_read_stack": "yh"} }' \
			 --out_dir .
		gzip ${name}.tc_hist.tsv
		gzip ${name}.tcreadpos_hist.tsv
		bgzip ${name}.tc_snps.vcf	
    	mkdir -p stats
		mv *tc_stats.tsv.gz stats
		samtools flagstat ${name}_tc_only.bam > stats/${name}_tc_only.bam.flagstat
	    samtools flagstat ${name}+tc.bam > stats/${name}+tc.bam.flagstat    
	    """
	} 	

/* TODO: usef featurecounts */
/*process count_tc_reads {
# count
python ${params.splice_sim_count_cmd}  count \
	 --bam ${name}+tc.bam \
	 --gff ${params.gene_gff} \
	 --threads ${task.cpus} \
	 --config_string '{ "output_format": "wide", "dataset_names": [ "${name}" ], "min_mapping_quality": 0, "bam_tags": { "num_tc": "yc","num_tt": "yt","density_tc": "yd","pos_tc": "yp","is_tc_duplicate": "yx","is_high_read_stack": "yh"} }' \
	 --out ${name}.counts.tsv
gzip ${name}.counts.tsv
}*/
