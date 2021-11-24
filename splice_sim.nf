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
    publishDir "reference_model", mode: 'copy'
    input: 
    	
    output: 
    	file("${params.dataset_name}.model") into model_file
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
	tag "$name" // custom label for trace
    cpus 1
    memory '64 GB'
    time 2.h
    //cache false 
    module 'art/2016.06.05-gcc-7.3.0-2.30:python/3.7.2-gcccore-8.2.0:htslib/1.10.2-gcccore-7.3.0:samtools/1.10-foss-2018b'
    publishDir "bams_truth", mode: 'copy'
    input: 
    	file(isoform_fasta) from isoform_fasta
    	file(model_file) from model_file
    output: 
    	set file("${params.dataset_name}.ori.fq.gz"), file("${params.dataset_name}.mod.fq.gz") into fq1, fq2, fq3, fq4, fq5
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
	    	
	    	# create genome bam + fq
	    	${params.splice_sim_cmd} create_genome_bam \
	    		--model ${model_file} --art_sam ${params.dataset_name}.sam --outdir . --threads ${task.cpus}
	    	
	    	# compress fastq
	    	bgzip --threads ${task.cpus} ${params.dataset_name}.ori.fq 
	    	bgzip --threads ${task.cpus} ${params.dataset_name}.mod.fq 
	    	
	    	# add md tags and update index
	    	samtools calmd --threads ${task.cpus} -b ${params.dataset_name}.ori.bam ${params.genome_fa} > ${params.dataset_name}.ori.tmp.bam && mv ${params.dataset_name}.ori.tmp.bam ${params.dataset_name}.ori.bam  
	    	samtools calmd --threads ${task.cpus} -b ${params.dataset_name}.mod.bam ${params.genome_fa} > ${params.dataset_name}.mod.tmp.bam && mv ${params.dataset_name}.mod.tmp.bam ${params.dataset_name}.mod.bam 
	    	samtools index ${params.dataset_name}.ori.bam
	    	samtools index ${params.dataset_name}.mod.bam
	    """
	} 


/*
 * STAR
 */
process map_star {

	tag "$name" // custom label for trace
    cpus 1
    memory '64 GB'
    time 2.h
    //cache false 
    module 'star/2.7.1a-foss-2018b:sambamba/0.6.6'
    publishDir "bams_star", mode: 'copy'
    
    when: 
    	params.mappers.STAR
    input: 
    	file(fq) from fq1.flatten()
    output: 
    	file("${fq.getBaseName(2)}.bam") into bams_star
    	file("${fq.getBaseName(2)}.bam.bai") into bais_star
    	file("*out") into logs_star
    script:
	    """
	       # map with STAR
	       ${params.mappers.STAR.star_cmd} --readFilesCommand zcat \
		      	--outFileNamePrefix ${fq.getBaseName(2)} \
		     	--runThreadN ${task.cpus} \
		      	--genomeDir ${params.mappers.STAR.star_genome_idx} \
		      	--readFilesIn ${fq} \
		      	--outSAMattributes NH HI AS nM MD \
		      	--outReadsUnmapped Fastx \
		      	--sjdbGTFfile ${params.mappers.STAR.star_splice_gtf}

		   sambamba view -S -f bam -t ${task.cpus}  ${fq.getBaseName(2)}Aligned.out.sam -o ${fq.getBaseName(2)}.pre.bam
		   sambamba sort -t ${task.cpus}  -o ${fq.getBaseName(2)}.bam ${fq.getBaseName(2)}.pre.bam 		      
	    """
	} 

/*
 * HISAT-3N
 */
process map_hisat_3n {

	tag "$name" // custom label for trace
    cpus 1
    memory '64 GB'
    time 2.h
    //cache false 
    module 'python/3.7.2-gcccore-8.2.0:sambamba/0.6.6'
    publishDir "bams_hisat3n", mode: 'copy'
    
    when: 
    	params.mappers.HISAT3N
    input: 
    	file(fq) from fq1.flatten()
    output: 
    	file("${fq.getBaseName(2)}.bam") into bams_hisat3n
    	file("${fq.getBaseName(2)}.bam.bai") into bais_hisat3n
    script:
	    """
	       gunzip -c ${fq} > ${fq.getBaseName(2)}.fq
	       # map with HISAT 3N
	       ${params.mappers.HISAT3N.hisat3n_cmd} \
	       		--base-change ${params.condition.ref},${params.condition.alt} \
	       		--index ${params.mappers.HISAT3N.hisat3n_idx} \
	       		-U ${fq.getBaseName(2)}.fq \
	       		--threads ${task.cpus} \
	       		--known-splicesite-infile ${params.mappers.HISAT3N.hisat3n_kss} \
	       		-S ${fq.getBaseName(2)}.sam
	      
		   sambamba view -S -f bam -t ${task.cpus}  ${fq.getBaseName(2)}.sam -o ${fq.getBaseName(2)}.pre.bam
		   sambamba sort -t ${task.cpus}  -o ${fq.getBaseName(2)}.bam ${fq.getBaseName(2)}.pre.bam        
	    """
	} 


/*
 * MERANGS
 */
process map_merangs {

	tag "$name" // custom label for trace
    cpus 1
    memory '64 GB'
    time 2.h
    //cache false 
    // NOTE that merangs does not support newer versions of STAR
    module 'star/2.5.2a-foss-2018b:python/3.7.2-gcccore-8.2.0:sambamba/0.6.6:biobambam2/2.0.87-foss-2018b:samtools/1.10-foss-2018b'
    publishDir "bams_merangs", mode: 'copy'
    
    when: 
    	params.mappers.MERANGS
    input: 
    	file(fq) from fq1.flatten()
    output: 
    	file("${fq.getBaseName(2)}.bam") into bams_merangs
    	file("${fq.getBaseName(2)}.bam.bai") into bais_merangs
    	file("*out") into logs_merangs
    script:
	    """
	    	mkdir unaligned_reads
	       # map with MERANGS
	       ${params.mappers.MERANGS.merangs_cmd} align \
	       		-o . \
	       		-f ${fq} \
	       		-t ${task.cpus} \
	       		-S ${fq.getBaseName(2)}.sam \
	       		-un -ud unaligned_reads \
	       		-id ${params.mappers.MERANGS.merangs_genome_idx} \
	       		-mbp \
	       		-starcmd ${params.mappers.MERANGS.star_cmd} \
	       		-star_outSAMattributes NH HI AS nM MD \
	       		-star_sjdbGTFfile ${params.mappers.MERANGS.merangs_splice_gtf}
	      
	       # converted unaligned FASTQ to BAM
	       fastqtobam gz=1 unaligned_reads/${fq.getBaseName(2)}_unmapped.fq.gz > ${fq.getBaseName(2)}_unmapped.bam
	       # merge aligned+unaligned reads
	       samtools merge ${fq.getBaseName(2)}.bam ${fq.getBaseName(2)}_sorted.bam ${fq.getBaseName(2)}_unmapped.bam	      
	    """
	} 
