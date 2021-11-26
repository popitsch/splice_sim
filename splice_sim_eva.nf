#!/usr/bin/env nextflow

params.config_file = "${workflow.launchDir}/splice_sim.config.json"
params.model=        "${workflow.launchDir}/sim/reference_model/${params.dataset_name}.model"
params.final_bam_dir="${workflow.launchDir}/sim/final_bams/"
params.truth_bam_dir="${workflow.launchDir}/sim/bams_truth/"

final_bams = Channel.fromFilePairs("${params.final_bam_dir}*.{bam,bai}", flat:true) { file -> file.name.replaceAll(/.bam|.bai$/,'') }
truth_bams = Channel.fromFilePairs("${params.truth_bam_dir}*.{bam,bai}", flat:true) { file -> file.name.replaceAll(/.bam|.bai$/,'') }
all_bams=final_bams.mix(truth_bams)

log.info "====================================="
log.info "Config file : ${params.config_file}"
log.info "Dataset : ${params.dataset_name}"
log.info "final_bam_dir : ${params.final_bam_dir}"
log.info "truth_bam_dir : ${params.truth_bam_dir}"
log.info "====================================="
log.info "\n"

/*
 * evaluate_overall_performance
 */
process evaluate_bam_performance {
	tag "$name" 
    cpus 1
    time 2.h
    module 'python/3.7.2-gcccore-8.2.0:sambamba/0.6.6'
    publishDir "eva/overall_performance", mode: 'copy'
    input:   
    	set name, file(bam), file(bai) from all_bams
    output: 
    	file("*") into overall_performance
    script:
	    """
    		${params.splice_sim_cmd} evaluate_bam_performance \
    			--config ${params.config_file} \
    			--model ${params.model} \
    			--bam_file ${bam} \
    			--outdir .
    			
    		# sort output bams
    		shopt -s nullglob
    		for bam in *.mismapped.bam.tmp.bam
	    	do
	    		final_bam="\${bam%.tmp.bam}"
	    		sambamba sort -t ${task.cpus}  -o \${final_bam} \${bam} && rm \${bam}
	    	done
    	    
	    """
	} 	

