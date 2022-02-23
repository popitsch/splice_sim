#!/usr/bin/env nextflow

params.final_bam_dir="${workflow.launchDir}/sim/final_bams/"
params.truth_bam_dir="${workflow.launchDir}/sim/bams_truth/"

Channel.fromFilePairs("${params.final_bam_dir}*.{bam,bai}", flat:true) { file -> file.name.replaceAll(/.bam|.bai$/,'') }.
into { final_bams; final_bams2 }

Channel.fromFilePairs("${params.truth_bam_dir}*.{bam,bai}", flat:true) { file -> file.name.replaceAll(/.bam|.bai$/,'') }.
into { truth_bams }

final_bams.mix(truth_bams).
into { all_bams; all_bams2 }

/*
 * count tc reads
 */
process count_tc_reads {
	tag "$name" 
    cpus 1
    module 'python/3.7.2-gcccore-8.2.0:sambamba/0.6.6:htslib/1.10.2-gcccore-7.3.0'
    publishDir "counts", mode: 'copy'
    input:   
    	set name, file(bam), file(bai) from all_bams
    output: 
    	file("*") into counts
    script:
	    """
	    # highlight  TC conversions
	    python /groups/ameres/Niko/workspace/genomic_iterators/tools/annotate_tc_conversions.py annotate \
    		 --bam ${bam} \
			 --threads ${task.cpus} \
			 --config_string '{ "region_width": 50000 }' \
			 --out_dir .
		# count
    	python /groups/ameres/Niko/workspace/genomic_iterators/tools/annotate_tc_conversions.py count \
			 --bam ${name}+tc.bam \
			 --gff ${params.gene_gff} \
			 --threads ${task.cpus} \
			 --config_string '{ "output_format": "wide", "dataset_names": [ "${name}" ], "min_mapping_quality": 0 }' \
			 --out ${name}.counts.tsv
		
    	gzip ${name}.counts.tsv
		gzip ${name}.tc_hist.tsv
		gzip ${name}.tcreadpos_hist.tsv
		bgzip ${name}.tc_snps.vcf	    
	    """
	} 	
