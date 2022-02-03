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
 * evaluate_overall_performance
 */
process evaluate_bam_performance {
	tag "$name"
    cpus 1
    module 'python/3.7.2-gcccore-8.2.0:sambamba/0.6.6'
    publishDir "eva/overall_performance", mode: 'copy'
    input:
    	set name, file(bam), file(bai) from all_bams
    output:
    	file("*") into bam_performance
    	set name, file("*mismapped.bam"),file("*mismapped.bam.bai") optional true into mismapped_bams
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
    cpus 1
    time 8.h
    module 'python/3.7.2-gcccore-8.2.0:sambamba/0.6.6'
    publishDir "eva/splice_site_performance", mode: 'copy'
    //cache false
    input:
    	set name, file(bam), file(bai) from all_bams2
    output:
    	file("*") into splice_sites_performance
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
    cpus 1
    time 8.h
    label "bigmem"
    module 'python/3.7.2-gcccore-8.2.0:bedtools/2.27.1-foss-2018b:htslib/1.10.2-gcccore-7.3.0'
    publishDir "eva/splice_site_mappability", mode: 'copy'
    input:
    	set name, file(bam), file(bai) from final_bams2
    output:
    	file("*") into splice_site_mappability
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
	 * Prepare rseqc input bams
	 */
	Channel.
		fromPath( 'sim/final_bams/*bam' ).
		map { file -> tuple(file.baseName, file) }
		into { rseqc_deletion_input ; rseqc_insertion_input}

	/*
	 * Prepare rseqc input bed
	 */
	process gff_to_bed {
		tag "$name"
	    cpus 1
	    module 'bedops/2.4.35'
	    cache false

			when:
    	params.rseqc

	    input:
	    	set name, file(bam), file(bai) from mismapped_bams
	    output:
	    	file("gene_anno.bed") into rseqcbed_channel
	    script:
		    """
				  gunzip -c ${params.gene_gff} > gene_anno.gff3
				  gff2bed < gene_anno.gff3 > gene_anno.bed
		    """
		}

		/*
		 * Calculate rseqc deletion profile
		 */
		process rseqc_deletion_profile {
			tag "$name"
	    cpus 1
	    module 'rseqc/2.6.5-foss-2018b-python-2.7.15'
			publishDir "eva/rseqc", mode: 'copy'
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
