#!/usr/bin/env nextflow

params.config_file = "${workflow.workDir}/splice_sim.config.json"

log.info "====================================="
log.info "Config file : ${params.config_file}"
log.info "Dataset : ${params.dataset_name}"
log.info "====================================="
log.info "\n"

/*
 * build model
 */
process build_model {
	tag "$name" // custom label for trace
    cpus 1
    memory '64 GB'
    time 2.h
    //cache false 
    module 'python/3.7.2-gcccore-8.2.0'
    publishDir "01_model", mode: 'copy'
    input: 
    output: 
    	set file("${params.dataset_name}.model") into model

    	script:
	    """
    		${params.executables.splice_sim_cmd} build_model --config ${params.config_file} --outdir .
	    """
	} 	
