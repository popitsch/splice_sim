#!/usr/bin/bash
ml purge &> /dev/null
ml load nextflow/20.10.0

export NXF_WORK="/scratch/nextflow/"$USER"/nxf_work"
export NXF_TEMP="/scratch/nextflow/"$USER"/nxf_temp"
export NXF_ANSI_LOG=false
mkdir -p logs
nextflow -log logs/nextflow.log run splice_sim.nf -params-file splice_sim.config.json -resume
