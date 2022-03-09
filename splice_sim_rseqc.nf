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
	}*/

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
/*process rseqc_deletion_profile {
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
 */
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
