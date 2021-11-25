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
    cpus 1
    time 2.h
    module 'python/3.7.2-gcccore-8.2.0'
    publishDir "sim/reference_model", mode: 'copy'
    input:    	
    output: 
    	file("${params.dataset_name}.model") into model_file1, model_file2
    	file("${params.dataset_name}.fa") into isoform_fasta
    	file("gene_anno.gff3.gz") into gff3
    	file("gene_anno.gff3.gz.tbi") into gff3_tbi
    	file("*") into model_files
    script:
	    """
    		${params.splice_sim_cmd} build_model --config ${params.config_file} --outdir .
	    """
	} 	


/*
 * simulate reads
 */
process simulate_reads {
    module 'art/2016.06.05-gcc-7.3.0-2.30:python/3.7.2-gcccore-8.2.0:htslib/1.10.2-gcccore-7.3.0:samtools/1.10-foss-2018b'
    publishDir "sim/bams_truth", mode: 'copy'
    input: 
    	file(isoform_fasta) from isoform_fasta
    	file(model_file) from model_file1
    output: 
    	file("*.fq.gz") into fq1, fq2, fq3, fq4, fq5
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
	    	art_illumina -ss HS25 -i ${isoform_fasta} -l ${params.readlen} -f ${params.double_cov} \
	    	-na --samout -o ${params.dataset_name} ${params.seed} ${params.additional_art_params} > art.log
	    	# drop fastq - we dont need it
	    	rm ${params.dataset_name}.fq
	    	
	    	# create genome bam + fq
	    	${params.splice_sim_cmd} create_genome_bam \
	    		--model ${model_file} --art_sam ${params.dataset_name}.sam --outdir . --threads ${task.cpus}
	    	
	    	# compress fastqs
	    	for fq in ${params.dataset_name}.cr*.fq
	    	do
	    	    bgzip --threads ${task.cpus} \${fq}
	    	done
	    	
	    	# add md tags and update index
	    	for bam in *.bam
	    	do
	    		samtools calmd --threads ${task.cpus} -b \${bam} ${params.genome_fa} > tmp.bam && mv tmp.bam \${bam}
	    		samtools index \${bam}
	    	done
	    """
	} 


/*
 * STAR
 */
process map_star {
	tag "$fq" 
	module 'star/2.7.1a-foss-2018b:sambamba/0.6.6'
    publishDir "sim/bams_star", mode: 'copy'
    
    when: 
    	params.mappers.STAR
    input: 
    	file(fq) from fq1.flatten()
    output: 
    	file("${fq.getBaseName(3)}.STAR.bam") into bams_star
    	file("${fq.getBaseName(3)}.STAR.bam.bai") into bais_star
    	file("*out") into logs_star
    script:
	    """
	       # map with STAR
	       ${params.mappers.STAR.star_cmd} --readFilesCommand zcat \
		      	--outFileNamePrefix ${fq.getBaseName(3)} \
		     	--runThreadN ${task.cpus} \
		      	--genomeDir ${params.mappers.STAR.star_genome_idx} \
		      	--readFilesIn ${fq} \
		      	--outSAMattributes NH HI AS nM MD \
		      	--outReadsUnmapped Fastx \
		      	--sjdbGTFfile ${params.mappers.STAR.star_splice_gtf}

		   sambamba view -S -f bam -t ${task.cpus}  ${fq.getBaseName(3)}Aligned.out.sam -o ${fq.getBaseName(3)}.pre.bam
		   sambamba sort -t ${task.cpus}  -o ${fq.getBaseName(3)}.STAR.bam ${fq.getBaseName(3)}.pre.bam 		      
	    """
	} 

/*
 * HISAT-3N
 */
process map_hisat_3n {
	tag "$fq" 
    module 'python/3.7.2-gcccore-8.2.0:sambamba/0.6.6'
    publishDir "sim/bams_hisat3n", mode: 'copy'
    
    when: 
    	params.mappers.HISAT3N
    input: 
    	file(fq) from fq2.flatten()
    output: 
    	file("${fq.getBaseName(3)}.HISAT3N.bam") into bams_hisat3n
    	file("${fq.getBaseName(3)}.HISAT3N.bam.bai") into bais_hisat3n
    script:
	    """
	       gunzip -c ${fq} > ${fq.getBaseName(3)}.fq
	       # map with HISAT 3N
	       ${params.mappers.HISAT3N.hisat3n_cmd} \
	       		--base-change ${params.condition.ref},${params.condition.alt} \
	       		--index ${params.mappers.HISAT3N.hisat3n_idx} \
	       		-U ${fq.getBaseName(3)}.fq \
	       		--threads ${task.cpus} \
	       		--known-splicesite-infile ${params.mappers.HISAT3N.hisat3n_kss} \
	       		-S ${fq.getBaseName(3)}.sam
	      
		   sambamba view -S -f bam -t ${task.cpus}  ${fq.getBaseName(3)}.sam -o ${fq.getBaseName(3)}.pre.bam
		   sambamba sort -t ${task.cpus}  -o ${fq.getBaseName(3)}.HISAT3N.bam ${fq.getBaseName(3)}.pre.bam        
	    """
	} 


/*
 * MERANGS
 */
process map_merangs {
	tag "$fq" 
    // NOTE that merangs does not support newer versions of STAR
    module 'star/2.5.2a-foss-2018b:python/3.7.2-gcccore-8.2.0:sambamba/0.6.6:biobambam2/2.0.87-foss-2018b:samtools/1.10-foss-2018b'
    publishDir "sim/bams_merangs", mode: 'copy'
    
    when: 
    	params.mappers.MERANGS
    input: 
    	file(fq) from fq3.flatten()
    output: 
    	file("${fq.getBaseName(3)}.MERANGS.bam") into bams_merangs
    	file("${fq.getBaseName(3)}.MERANGS.bam.bai") into bais_merangs
    script:
	    """
	    	mkdir unaligned_reads
	       # map with MERANGS
	       ${params.mappers.MERANGS.merangs_cmd} align \
	       		-o . \
	       		-f ${fq} \
	       		-t ${task.cpus} \
	       		-S ${fq.getBaseName(3)}.sam \
	       		-un -ud unaligned_reads \
	       		-id ${params.mappers.MERANGS.merangs_genome_idx} \
	       		-mbp \
	       		-starcmd ${params.mappers.MERANGS.star_cmd} \
	       		-star_outSAMattributes NH HI AS nM MD \
	       		-star_sjdbGTFfile ${params.mappers.MERANGS.merangs_splice_gtf}
	      
	       # converted unaligned FASTQ to BAM
	       fastqtobam gz=1 unaligned_reads/${fq.getBaseName(3)}.truth_unmapped.fq.gz > ${fq.getBaseName(3)}_unmapped.bam
	       # merge aligned+unaligned reads
	       samtools merge ${fq.getBaseName(3)}.MERANGS.bam ${fq.getBaseName(3)}_sorted.bam ${fq.getBaseName(3)}_unmapped.bam
	       samtools index ${fq.getBaseName(3)}.MERANGS.bam
	    """
	} 


/*
 * Postprocessing
 */
process postprocess_bams {
	cpus 1
	memory '64 GB'
	time 2.h
	//cache false 
	module 'python/3.7.2-gcccore-8.2.0'
	publishDir "sim/final_bams", mode: 'copy'
	input: 
		file(bam) from bams_star.mix(bams_hisat3n, bams_merangs).ifEmpty([])
		file(bai) from bais_star.mix(bais_hisat3n, bais_merangs).ifEmpty([])
	output: 
		file("*.bam") into postprocessed_bams
		file("*.bai") into postprocessed_bais
	script:
    """
		${params.splice_sim_cmd} postfilter_bam --config ${params.config_file} --bam ${bam} --outdir .
    """
} 	



