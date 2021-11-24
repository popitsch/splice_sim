#!/usr/bin/env nextflow

params.config_file = "${workflow.launchDir}/splice_sim.config.json"

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
    	file("${params.dataset_name}.model") into model
    	file("${params.dataset_name}.fa") into isoform_fasta
    	file("gene_anno.gff3.gz") into gff3
    	file("gene_anno.gff3.gz.tbi") into gff3_tbi
    	file("*") into model_files
    script:
	    """
    		${params.executables.splice_sim_cmd} build_model --config ${params.config_file} --outdir .
	    """
	} 	


/*
 * simulate reads
 */
process simulate_reads {
	tag "$name" // custom label for trace
    cpus 1
    memory '64 GB'
    time 2.h
    //cache false 
    module 'art/2016.06.05-gcc-7.3.0-2.30:python/3.7.2-gcccore-8.2.0'
    publishDir "02_simulated_reads", mode: 'copy'
    input: 
    	file(isoform_fasta) from isoform_fasta
    output: 
    	file("*") into simulate_reads
    params.double_cov=params.condition.base_coverage*2
    if ( ! params.additional_art_params ) {
    	params.additional_art_params=""
    }
    if ( params.random_seed ) {
    	params.seed="-rs " + params.random_seed
    } else {
    	params.seed=""
    }
    script:
	    """
	    	${params.executables.art_cmd} -ss HS25 -i ${isoform_fasta} -l ${params.readlen} -f ${params.double_cov} \
	    	-na --samout -o ${params.dataset_name} ${params.seed} ${params.additional_art_params} > art.log
	    """
	} 

