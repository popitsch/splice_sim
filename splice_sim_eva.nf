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
	into { all_bams; all_bams2 }

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
process extract_transcript_features {
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
}

/*
 * extract splice-junction features
 */
process extract_splice_junction_features {
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
}

/*
 * evaluate_overall_performance
 */
process evaluate_bam_performance {
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
	}

/*
 * evaluate_splice_site_performance
 */
process evaluate_splice_site_performance {
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
	}

/*
 * calculate_splice_site_mappability
 */
process calculate_splice_site_mappability {
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
	}


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

	/*
	 * calc_feature_overlap
	 */
	process collect_splice_sim_eva_results {
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
	    output:
	    	file("results") into splice_sim_eva_results
	    script:
		    """
	    		${params.splice_sim_cmd} write_parquet_db \
	    			--config ${params.config_file} \
	    			--indir $PWD \
	    			--outdir results
		    """
		}

	/*
	 * Prepare rseqc input bams
	 */
	Channel.
		fromPath( "${workflow.launchDir}/sim/bams_truth/*bam" ).
		set{ truth_bams_channel }

	Channel.
		fromPath( "${workflow.launchDir}/sim/final_bams/*bam" ).
		set{ mapped_bams_channel }

	truth_bams_channel.
	  mix(mapped_bams_channel).
		map { file -> tuple(file.baseName, file) }.
		into { rseqc_clipping_input;
			rseqc_mismatch_input;
			rseqc_deletion_input;
			rseqc_insertion_input;
			rseqc_junction_annotation_input;
			rseqc_junction_saturation_input;
			rseqc_read_distribution_input ;
			rseqc_read_duplication_input ;
			rseqc_read_quality_input
		}

	Channel.
		fromPath( "${workflow.launchDir}/sim/bams_truth/*bam*" ).
		set{ truth_bams_bais_channel }

	Channel.
		fromPath( "${workflow.launchDir}/sim/final_bams/*bam*" ).
		set{ mapped_bams_bais_channel }

	truth_bams_bais_channel.
		mix(mapped_bams_bais_channel).
		set { rseqc_genebody_input }

	/*
	 * Prepare rseqc input bed
	 */
	process gff_to_bed {
		tag "${params.dataset_name}"
	    cpus 1
	    module 'transdecoder/5.5.0-foss-2018b-perl-5.28.0'
	    cache false

			when:
    	params.rseqc

	    input:

	    output:
	    	file("gene_anno.bed") into rseqcbed_genebody_profile, rseqcbed_junction_saturation, rseqcbed_junction_annotation, rseqcbed_read_distribution
	    script:
		    """
				  gunzip -c ${params.gene_gff} > gene_anno.gff3
					sed -i 's/\ttranscript\t/\tmRNA\t/g' gene_anno.gff3
					gff3_file_to_bed.pl gene_anno.gff3 > gene_anno.bed
					sed -i '1d' gene_anno.bed
				  #gff2bed < gene_anno.gff3 > gene_anno.bed
		    """
		}

		/*
		 * Calculate rseqc clipping profile
		 */
		process rseqc_clipping_profile {
			tag "$name"
			cpus 1
			module 'rseqc/2.6.5-foss-2018b-python-2.7.15'
			publishDir "eva/rseqc/clipping_profile", mode: 'copy'
			cache false

			when:
			params.rseqc

			input:
				set val(name), file(bam) from rseqc_clipping_input
			output:
				file("${name}*") into rseqc_clipping_output
			script:
				"""
				  clipping_profile.py -i $bam -o $name -s SE
				"""
		}

		/*
		 * Calculate rseqc deletion profile
		 */
		process rseqc_mismatch_profile {
			tag "$name"
			cpus 1
			module 'rseqc/2.6.5-foss-2018b-python-2.7.15'
			publishDir "eva/rseqc/mismatch_profile", mode: 'copy'
			cache false

			when:
			params.rseqc

			input:
				set val(name), file(bam) from rseqc_mismatch_input
			output:
				file("${name}*") into rseqc_mismatch_output
			script:
				"""
				  mismatch_profile.py -i $bam -o $name -l ${params.readlen}
				"""
		}

		/*
		 * Calculate rseqc deletion profile
		 */
		process rseqc_insertion_profile {
			tag "$name"
	    cpus 1
	    module 'rseqc/2.6.5-foss-2018b-python-2.7.15'
			publishDir "eva/rseqc/insertion_profile", mode: 'copy'
	    cache false

			when:
    	params.rseqc

	    input:
	    	set val(name), file(bam) from rseqc_insertion_input
	    output:
	    	file("${name}*") into rseqc_insertion_output
	    script:
		    """
				  insertion_profile.py -i $bam -o $name -s SE
		    """
		}

		/*
		 * Calculate rseqc deletion profile
		 */
		process rseqc_deletion_profile {
			tag "$name"
	    cpus 1
	    module 'rseqc/2.6.5-foss-2018b-python-2.7.15'
			publishDir "eva/rseqc/deletion_profile", mode: 'copy'
	    cache false

			when:
    	params.rseqc

	    input:
	    	set val(name), file(bam) from rseqc_deletion_input
	    output:
	    	file("${name}*") into rseqc_deletion_output
	    script:
		    """
				  deletion_profile.py -i $bam -o $name -l ${params.readlen}
		    """
		}

		/*
		 * Calculate rseqc read duplication
		 */
		process rseqc_read_duplication {
			tag "$name"
	    cpus 1
	    module 'rseqc/2.6.5-foss-2018b-python-2.7.15'
			publishDir "eva/rseqc/read_duplication", mode: 'copy'
	    cache false

			label 'highmem'

			when:
    	params.rseqc

	    input:
	    	set val(name), file(bam) from rseqc_read_duplication_input
	    output:
	    	file("${name}*") into rseqc_read_duplication_output
	    script:
		    """
				  read_duplication.py -i $bam -o $name
		    """
		}

		/*
		 * Calculate rseqc read distribution
		 */
		process rseqc_read_distribution {
			tag "$name"
			cpus 1
			module 'rseqc/2.6.5-foss-2018b-python-2.7.15'
			publishDir "eva/rseqc/read_distribution", mode: 'copy'
			cache false

			when:
			params.rseqc

			input:
				set val(name), file(bam) from rseqc_read_distribution_input
				file(reference) from rseqcbed_read_distribution
			output:
				file("${name}*") into rseqc_read_distribution_output
			script:
				"""
				  read_distribution.py -i $bam -r $reference > ${name}_readDistribution.txt
				"""
			}

			/*
			 * Calculate rseqc read quality
			 *
			process rseqc_read_quality {
				tag "$name"
				cpus 1
				module 'rseqc/2.6.5-foss-2018b-python-2.7.15'
				publishDir "eva/rseqc/read_quality", mode: 'copy'
				cache false

				label 'highmem'

				when:
				params.rseqc

				input:
					set val(name), file(bam) from rseqc_read_quality_input
				output:
					file("${name}*") into rseqc_read_quality_output
				script:
					"""
					  read_quality.py -i $bam -o $name
					"""
				}
				*/

			/*
			 * Calculate rseqc junction annotation
			 */
			process rseqc_junction_annotation {
				tag "$name"
		    cpus 1
		    module 'rseqc/2.6.5-foss-2018b-python-2.7.15'
				publishDir "eva/rseqc/junction_annotation", mode: 'copy'
		    cache false

				when:
	    	params.rseqc

		    input:
		    	set val(name), file(bam) from rseqc_junction_annotation_input
					file(reference) from rseqcbed_junction_annotation
		    output:
		    	file("${name}*") into rseqc_junction_annotation_output
		    script:
			    """
					  junction_annotation.py -i $bam -r $reference -o $name 2> ${name}_junction_annotation.txt
			    """
			}

			/*
			 * Calculate rseqc junction saturation
			 */
			process rseqc_junction_saturation {
				tag "$name"
		    cpus 1
		    module 'rseqc/2.6.5-foss-2018b-python-2.7.15'
				publishDir "eva/rseqc/junction_saturation", mode: 'copy'
		    cache false

				when:
	    	params.rseqc

		    input:
		    	set val(name), file(bam) from rseqc_junction_saturation_input
					file(reference) from rseqcbed_junction_saturation
		    output:
		    	file("${name}*") into rseqc_junction_saturation_output
		    script:
			    """
					  junction_saturation.py -i $bam -r $reference -o $name
			    """
			}

			/*
			 * Calculate rseqc genebody profile
			 */
			process rseqc_genebody_profile {
				tag "${params.dataset_name}"
				time='8.h'
		    module 'rseqc/2.6.5-foss-2018b-python-2.7.15'
				publishDir "eva/rseqc/geneBody_coverage", mode: 'copy'
		    cache false

				label 'long'

				when:
	    	params.rseqc

		    input:
		    	file(bams) from rseqc_genebody_input.collect()
					file(reference) from rseqcbed_genebody_profile
		    output:
		    	file("${params.dataset_name}*") into rseqc_genebody_output
		    script:
			    """
					  geneBody_coverage.py -r $reference -i . -o ${params.dataset_name}
			    """
				}

				/*
				 * MultiQC
				 */
				process multiqc {
					tag "${params.dataset_name}"
			    cpus 1
			    module 'multiqc/1.9-foss-2018b-python-3.6.6'
					publishDir "eva/multiqc", mode: 'copy'
			    cache false

					when:
		    	params.rseqc

			    input:
						file(rseqc_clipping) from rseqc_clipping_output.collect().ifEmpty([])
						file(rseqc_mismatch) from rseqc_mismatch_output.collect().ifEmpty([])
						file(rseqc_insertion) from rseqc_insertion_output.collect().ifEmpty([])
						file(rseqc_deletion) from rseqc_deletion_output.collect().ifEmpty([])
						file(rseqc_duplication) from rseqc_read_duplication_output.collect().ifEmpty([])
						//file(rseqc_read_distribution) from rseqc_read_distribution_output.collect().ifEmpty([])
						//file(rseqc_read_quality) from rseqc_read_quality_output.collect().ifEmpty([])
						file(rseqc_junction_annotation) from rseqc_junction_annotation_output.collect().ifEmpty([])
						file(rseqc_junction_saturation) from rseqc_junction_saturation_output.collect().ifEmpty([])
						file(rseqc_genebody_profile) from rseqc_genebody_output.collect().ifEmpty([])

			    output:
			    file "*multiqc_report.html" into ch_multiqc_report
			    file "*_data"

			    script:
			    """
			      multiqc -m rseqc -f .
			    """
				}
