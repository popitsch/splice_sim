#!/usr/bin/env nextflow

Channel.fromFilePairs("${params.bam_dir}/*.{bam,bai}", flat:true) { file -> file.name.replaceAll(/.bam|.bai$/,'') }.set{ input_bams }

/*
 * count tc reads
 */
process count_tc_reads {
	tag "$name" 
    cpus 8
    module 'python/3.7.2-gcccore-8.2.0:sambamba/0.6.6:htslib/1.10.2-gcccore-7.3.0'
    publishDir "counts", mode: 'copy'
    input:   
    	set name, file(bam), file(bai) from input_bams
    output: 
    	file("*") into counts
    script:
	    """
	    # highlight  TC conversions
	    python ${params.splice_sim_count_cmd} annotate \
    		 --bam ${bam} \
			 --threads ${task.cpus} \
			 --config_string '{ "region_width": 50000, "bam_tags": { "num_tc": "yc","num_tt": "yt","density_tc": "yd","pos_tc": "yp","is_tc_duplicate": "yx","is_high_read_stack": "yh"} }' \
			 --out_dir .
		# count
    	python ${params.splice_sim_count_cmd}  count \
			 --bam ${name}+tc.bam \
			 --gff ${params.gene_gff} \
			 --threads ${task.cpus} \
			 --config_string '{ "output_format": "wide", "dataset_names": [ "${name}" ], "min_mapping_quality": 0, "bam_tags": { "num_tc": "yc","num_tt": "yt","density_tc": "yd","pos_tc": "yp","is_tc_duplicate": "yx","is_high_read_stack": "yh"} }' \
			 --out ${name}.counts.tsv
		
    	gzip ${name}.counts.tsv
		gzip ${name}.tc_hist.tsv
		gzip ${name}.tcreadpos_hist.tsv
		bgzip ${name}.tc_snps.vcf	    
	    """
	} 	
