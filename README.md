# splice_sim

`splice_sim` is a python based RNA-seq simulation and evaluation package specialized on the simulation nucleotide-conversions and splice isoforms. We use [Nextflow](https://www.nextflow.io/) and containerization via [Docker](https://hub.docker.com/repository/docker/tobneu/splice_sim) to wrap `splice_sim` into a readily executable start-to-end framework for simulating and evaluating complex experimental scenarios with mapping based processing pipelines.

Contents
========

 * [Features](#features)
 * [Installation](#installation)
 * [Usage](#usage)
 * [Configuration](#configuration)
 * [Output Structure](#output-structure)

### Features

+ Realistic Illumina short-read simulation using [ART](https://doi.org/10.1093/bioinformatics/btr708)
+ Simulation of customizable nucleotide-conversions at configurable rates
+ Simulation of isoforms at configurable splicing states per transcript
+ Mapping accuracies and performance metrics at different scopes of genomic annotation
+ Elaborate output tracks for visual inspection stratified by performance metric

### Installation
---

### Usage
---

### Configuration
---
Here's an example config :

```json
{
    "dataset_name": "test_simulation",
    "splice_sim_cmd": "python /software/splice_sim/main.py",
    "splice_eva_preprocess_cmd": "Rscript --vanilla /software/splice_sim/src/main/R/splice_sim/preprocess_results.R",
    "gene_gff": "/references/gencode.vM21.gff3.gz",
    "intron_gff": "/references/gencode.vM21.introns.sorted.gff3.gz",
    "genome_fa": "/references/Mus_musculus.GRCm38.dna.primary_assembly.fa",
    "genome_chromosome_sizes": "/references/Mus_musculus.GRCm38.dna.primary_assembly.fa.chrom.sizes",
    "genome_conservation": "/references/mm10.60way.phastCons60wayEuarchontoGlire.bw",
    "genome_mappability": "/references/mm10.k24.umap.bedgraph.gz",
    "transcript_data": "data.config.json",
    "transcript_ids": "tids.tsv",
    "isoform_mode": "1:1",
    "frac_old_mature": 0,
    "condition": {
        "ref": "T",
        "alt": "C",
        "conversion_rates": [ 0.02, 0.04 ],
        "base_coverage": 10
    },
    "mappers": {
        "STAR": {
            "star_cmd": "STAR-2.7.1a",
            "star_genome_idx": "star_2.7.1",
            "star_splice_gtf": "/indices/gencode.vM21.gtf"
            },
        "HISAT3N": {
            "hisat3n_cmd": "hisat-3n",
            "hisat3n_idx": "/indices/Mus_musculus.GRCm38.dna.primary_assembly",
            "hisat3n_kss": "/indices/gencode.vM21.gtf.hisat2_splice_sites.txt"
            },
        "MERANGS": {
            "merangs_cmd": "meRanGs",
            "star_cmd": "STAR",
            "merangs_genome_idx": "/indices/meRanTK-1.2.1b/BSgenomeIDX",
            "merangs_splice_gtf": "/indices/gencode.vM21.gtf"
            }
        },
    "create_tdf": true,
    "max_ilen": 100000,
    "min_abundance": 1,
    "random_seed": 1234,
    "readlen": 100,
    "write_reads": false,
    "write_intron_bam": false,
    "#pre_exist
}
```

### Output Structure
---
