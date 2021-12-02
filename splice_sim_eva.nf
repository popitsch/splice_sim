#!/usr/bin/env nextflow

params.config_file = "${workflow.launchDir}/splice_sim.config.json"
params.model=        "${workflow.launchDir}/sim/reference_model/${params.dataset_name}.model"
params.final_bam_dir="${workflow.launchDir}/sim/final_bams/"
params.truth_bam_dir="${workflow.launchDir}/sim/bams_truth/"

final_bams = Channel.fromFilePairs("${params.final_bam_dir}*.{bam,bai}", flat:true) { file -> file.name.replaceAll(/.bam|.bai$/,'') }
truth_bams = Channel.fromFilePairs("${params.truth_bam_dir}*.{bam,bai}", flat:true) { file -> file.name.replaceAll(/.bam|.bai$/,'') }
all_bams=final_bams.mix(truth_bams)
final_bams2 = Channel.fromFilePairs("${params.final_bam_dir}*.{bam,bai}", flat:true) { file -> file.name.replaceAll(/.bam|.bai$/,'') }

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
    	file("*") into bam_performance
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

/*
 * evaluate_splice_sites
 */
process evaluate_splice_sites_performance {
	tag "$name" 
    cpus 1
    time 2.h
    module 'python/3.7.2-gcccore-8.2.0:sambamba/0.6.6'
    publishDir "eva/splice_sites_performance", mode: 'copy'
    input:   
    	set name, file(bam), file(bai) from all_bams
    output: 
    	file("*") into splice_sites_performance
    script:
	    """
    		${params.splice_sim_cmd} evaluate_splice_sites_performance \
    			--config ${params.config_file} \
    			--model ${params.model} \
    			--bam_file ${bam} \
    			--outdir .
    			
    		# sort output bams
    		shopt -s nullglob
    		for bam in *.intron.bam.tmp.bam
	    	do
	    		final_bam="\${bam%.tmp.bam}"
	    		sambamba sort -t ${task.cpus}  -o \${final_bam} \${bam} && rm \${bam}
	    	done
    	    
	    """
	} 	

/*
 * evaluate_splice_sites
 */
process calculate_splice_site_mappability {
	tag "$name" 
    cpus 1
    time 2.h
    module 'python/3.7.2-gcccore-8.2.0:bedtools/2.27.1-foss-2018b'
    publishDir "eva/splice_site_mappability", mode: 'copy'
    input:   
    	set name, file(bam), file(bai) from final_bams2
    output: 
    	file("*") into splice_site_mappability
    script:
	    """
    		${params.splice_sim_cmd} evaluate_splice_sites_performance \
    			--config ${params.config_file} \
    			--model ${params.model} \
    			--truth_bam_dir ${params.truth_bam_dir} \
    			--bam_file ${bam} \
    			--outdir .
    			
    		# sort output bed
    		shopt -s nullglob
    		for bed in *.bedGraph.tmp
	    	do
	    		final_bed="\${bed%.tmp}"
	    		bedtools sort -i \${bed} | bedtools genomecov -bg -i stdin -g ${params.genome_chromosome_sizes}
	    	done
    	    
	    """
	} 	


