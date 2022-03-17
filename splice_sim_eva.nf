#!/usr/bin/env nextflow

params.config_file = "${workflow.launchDir}/splice_sim.config.json"
params.model=        "${workflow.launchDir}/sim/reference_model/${params.dataset_name}.model"
params.final_bam_dir="${workflow.launchDir}/sim/final_bams/"
params.truth_bam_dir="${workflow.launchDir}/sim/bams_truth/"


Channel.fromFilePairs("${params.final_bam_dir}*.{bam,bai}", flat:true) { file -> file.name.replaceAll(/.bam|.bai$/,'') }.
	into { final_bams; final_bams2 }

Channel.fromFilePairs("${params.truth_bam_dir}*.{bam,bai}", flat:true) { file -> file.name.replaceAll(/.bam|.bai$/,'') }.
	into { truth_bams }

final_bams.mix(truth_bams).
	into { all_bams; all_bams2; all_bams3 }

log.info "====================================="
log.info "Config file : ${params.config_file}"
log.info "Dataset : ${params.dataset_name}"
log.info "final_bam_dir : ${params.final_bam_dir}"
log.info "truth_bam_dir : ${params.truth_bam_dir}"
log.info "====================================="
log.info "\n"

/*
 * extract transcript features
 */
/*process extract_transcript_features {
	tag "${params.dataset_name}"
    module 'python/3.7.2-gcccore-8.2.0'
    publishDir "eva/meta", mode: 'copy'
    input:
		file(config) from Channel.fromPath("${params.config_file}")
    	file(model) from Channel.fromPath("${params.model}")
    output:
    	file("*transcript.metadata.tsv") into transcript_metadata
		val true into done_trancript_feature_extraction
    script:
	    """
    		${params.splice_sim_cmd} extract_transcript_features \
    			--config $config \
    			--model $model \
    			--outdir .
	    """
}*/

/*
 * extract splice-junction features
 */
/*process extract_splice_junction_features {
	tag "${params.dataset_name}"
    module 'python/3.7.2-gcccore-8.2.0'
    publishDir "eva/meta", mode: 'copy'
    input:
		file(config) from Channel.fromPath("${params.config_file}")
    	file(model) from Channel.fromPath("${params.model}")
    output:
    	file("*SJ.metadata.tsv") into sj_metadata
		val true into done_splice_junction_feature_extraction
    script:
	    """
    		${params.splice_sim_cmd} extract_splice_site_features \
    			--config $config \
    			--model $model \
    			--outdir .
	    """
}*/

/*
 * evaluate_overall_performance
 */
/*process evaluate_bam_performance {
	tag "$name"
	label "medium"
    module 'python/3.7.2-gcccore-8.2.0:sambamba/0.6.6'
    publishDir "eva/overall_performance", mode: 'copy'
    input:
    	set name, file(bam), file(bai) from all_bams
    output:
    	file("*") into bam_performance
    	set name, file("*mismapped.bam"),file("*mismapped.bam.bai") optional true into mismapped_bams
		val true into done_bam_performance
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
	}*/

/*
 * evaluate_splice_site_performance
 */
/*process evaluate_splice_site_performance {
	tag "$name"
    module 'python/3.7.2-gcccore-8.2.0:sambamba/0.6.6'
    publishDir "eva/splice_site_performance", mode: 'copy'
    //cache false
    input:
    	set name, file(bam), file(bai) from all_bams2
    output:
    	file("*") into splice_sites_performance
		val true into done_splice_site_performance
    script:
	    """
    		${params.splice_sim_cmd} evaluate_splice_site_performance \
    			--config ${params.config_file} \
    			--model ${params.model} \
    			--truth_bam_dir ${params.truth_bam_dir} \
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
	}*/

/*
 * calculate_splice_site_mappability
 */
/*process calculate_splice_site_mappability {
	tag "$name"
    module 'python/3.7.2-gcccore-8.2.0:bedtools/2.27.1-foss-2018b:htslib/1.10.2-gcccore-7.3.0'
    publishDir "eva/splice_site_mappability", mode: 'copy'
    input:
    	set name, file(bam), file(bai) from final_bams2
    output:
    	file("*") into splice_site_mappability
		val true into done_splice_site_mappability
    script:
	    """
    		${params.splice_sim_cmd} calculate_splice_site_mappability \
    			--config ${params.config_file} \
    			--model ${params.model} \
    			--truth_bam_dir ${params.truth_bam_dir} \
    			--bam_file ${bam} \
    			--outdir .

    		# sort output bedgraph files
    		shopt -s nullglob
    		for bed in *.bedGraph.tmp
	    	do
	    		final_bed="\${bed%.tmp}"
    			bedtools sort -i \${bed} > \${bed}.sorted
    			bedtools genomecov -bg -i \${bed}.sorted -g ${params.genome_chromosome_sizes} | bedtools map -a stdin -b \${bed}.sorted -o max | cut -f1,2,3,5 > \${final_bed} && rm \${bed}.sorted
    			bgzip \${final_bed} && tabix \${final_bed}.gz
	    	done

	    """
	}*/

Channel.fromFilePairs("${workflow.launchDir}/eva/overall_performance/*.{bam,bai}", flat:true) { file -> file.name.replaceAll(/.bam|.bai$/,'') }.
set { mismapped_bams }


/*
 * calc_feature_overlap
 */
process calc_feature_overlap {
	tag "$name"
    cpus 1
    module 'python/3.7.2-gcccore-8.2.0:sambamba/0.6.6'
    publishDir "eva/feature_overlap", mode: 'copy'
    cache false
    input:
    	set name, file(bam), file(bai) from mismapped_bams
    output:
    	file("*") into feature_overlaps
		val true into done_feature_overlap_mappability
    script:
	    """
    		${params.splice_sim_cmd} calc_feature_overlap \
    			--config ${params.config_file} \
    			--model ${params.model} \
    			--bam_file ${bam} \
    			--outdir .
	    """
	}


/*process run_featureCounts {
	tag "${params.dataset_name}"
	cpus 4
	module 'subread/2.0.1-gcc-7.3.0-2.30'
	publishDir "eva/feature_counts", mode: 'copy'
	cache false
	when:
		params.feature_counts
	input:
    	set name, file(bam), file(bai) from all_bams2
	output:
		file "*" into ch_feature_counts
		val true into done_run_featureCounts
	script:
	"""
	featureCounts -T ${task.cpus} \
		-a ${params.gene_gff} \
		-O \
		-M \
		--primary \
		-t exon \
		-f \
		-g ID \
		-o ${name}.exon.featureCounts.txt \
		-s 0 \
		${bam}
	featureCounts -T ${task.cpus} \
		-a ${params.intron_gff} \
		-O \
		-M \
		--primary \
		-t intron \
		-f \
		-g ID \
		-o ${name}.intron.featureCounts.txt \
		-s 0 \
		${bam}
	cp .command.log ${name}.featureCounts.log
	"""
}
*/
/*
 * collect_splice_sim_eva_results and build parquet db
 */
/*process collect_splice_sim_eva_results {
	tag "$name"
    cpus 1
    module 'python/3.7.2-gcccore-8.2.0:sambamba/0.6.6'
    publishDir ".", mode: 'copy'
    cache false
    input:
    	val flag1 from done_trancript_feature_extraction.collect()
    	val flag2 from done_splice_junction_feature_extraction.collect()
    	val flag3 from done_bam_performance.collect()
    	val flag4 from done_splice_site_performance.collect()
    	val flag5 from done_splice_site_mappability.collect()
    	val flag6 from done_feature_overlap_mappability.collect()
    	val flag7 from done_run_featureCounts.collect()
    output:
    	file("results") into splice_sim_eva_results
    script:
	    """
	    		${params.splice_sim_cmd} write_parquet_db \
	    			--config ${params.config_file} \
	    			--indir $PWD \
	    			--outdir results
		"""
	}*/