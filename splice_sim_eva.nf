#!/usr/bin/env nextflow

params.config_file = "${workflow.launchDir}/splice_sim.config.json"
params.model="${workflow.launchDir}/sim/reference_model/${params.dataset_name}.model"
params.final_bam_dir="${workflow.launchDir}/sim/final_bams/"
params.truth_bam_dir="${workflow.launchDir}/sim/bams_truth/"

final_bams = Channel.fromPath( "${params.final_bam_dir}/*.ba{m,i}" )
truth_bams = Channel.fromPath( "${params.truth_bam_dir}/*.ba{m,i}" )

log.info "====================================="
log.info "Config file : ${params.config_file}"
log.info "Dataset : ${params.dataset_name}"
log.info "====================================="
log.info "\n"

/*
 * evaluate_overall_performance
 */
process evaluate_overall_performance {
    cpus 1
    time 2.h
    module 'python/3.7.2-gcccore-8.2.0:sambamba/0.6.6'
    publishDir "eva/overall_performance", mode: 'copy'
    input:   
    	file("*") from final_bams.collect()
    	file("*") from truth_bams.collect()
    output: 
    	file("*") into overall_performance
    script:
	    """
    		${params.splice_sim_cmd} evaluate_overall_performance \
    			--config ${params.config_file} \
    			--model ${params.model} \
    			--final_bam_dir ${params.final_bam_dir} \
    			--truth_bam_dir ${params.truth_bam_dir} \
    			--outdir .
    			
    		# sort bams
    		for bam in *.mismapped.bam.tmp.bam
	    	do
	    		final_bam="\${bam%.tmp.bam}"
	    		sambamba sort -t ${task.cpus}  -o \${final_bam} \${bam} && rm \${bam}
	    	done
    	    
	    """
	} 	

