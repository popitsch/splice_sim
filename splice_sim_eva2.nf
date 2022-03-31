#!/usr/bin/env nextflow

params.config_file = "${workflow.launchDir}/splice_sim.config.json"
params.model=        "${workflow.launchDir}/sim/reference_model/${params.dataset_name}.model"
params.final_bam_dir="${workflow.launchDir}/sim/final_bams/"
params.truth_bam_dir="${workflow.launchDir}/sim/bams_truth/"
log.info "====================================="
log.info "Config file : ${params.config_file}"
log.info "Dataset : ${params.dataset_name}"
log.info "final_bam_dir : ${params.final_bam_dir}"
log.info "truth_bam_dir : ${params.truth_bam_dir}"
log.info "====================================="
log.info "\n"


Channel.fromFilePairs("${params.final_bam_dir}*.{bam,bai}", flat:true) { file -> file.name.replaceAll(/.bam|.bai$/,'') }.
	into { final_bams }


/*
 * evaluate
 */
process evaluate_bams {
    tag "$name"
    label "long"
    module 'python/3.7.2-gcccore-8.2.0:htslib/1.10.2-gcccore-7.3.0'
    publishDir "eva/counts", mode: 'copy'
    input:
    	set name, file(bam), file(bai) from final_bams
    output:
    	file("*") into evaluate_results
		val true into done_evaluate
    script:
	    """
    		${params.splice_sim_cmd} evaluate \
    			--config ${params.config_file} \
    			--model ${params.model} \
    			--bam_file ${bam} \
    			--threads ${task.cpus} \
    			--outdir .

    		# bgzip TSV tables
    		shopt -s nullglob
    		for tsv in *.tsv
	    	do
	    		bgzip \${tsv}
	    	done
	    """
	}


/*
 * calc metadata
 */
process extract_feature_metadata {
	tag "$name"
    cpus 1
    module 'python/3.7.2-gcccore-8.2.0:htslib/1.10.2-gcccore-7.3.0'
    publishDir "eva/meta", mode: 'copy'
    cache true
    output:
    	file("*") into extract_feature_metadata_results
    	val true into done_extract_feature_metadata
    script:
	    """
	    	${params.splice_sim_cmd} extract_feature_metadata \
				--config ${params.config_file} \
				--model ${params.model} \
				--outdir .
    		# bgzip TSV tables
    		shopt -s nullglob
    		for tsv in *.tsv
	    	do
	    		bgzip \${tsv}
	    	done
		"""
	}

/*
 * preprocess results
 */
process preprocess_results {
    cpus 1
    module 'r/4.0.2-foss-2018b'
    publishDir "eva/results", mode: 'copy'
    cache true
    input:
    	val flag1 from done_evaluate.collect()
    	val flag2 from done_extract_feature_metadata.collect()
    output:
    	file("*") into preprocess_results_results
    script:
	    """
	    	${params.splice_eva_preprocess_cmd} ${params.config_file} . 
		"""
}
