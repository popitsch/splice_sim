# splice_sim

`splice_sim` is a python based RNA-seq simulation and evaluation package specialized on the simulation nucleotide-conversions and splice isoforms. We use [Nextflow](https://www.nextflow.io/) and containerization via [Docker](https://hub.docker.com/repository/docker/tobneu/splice_sim) to wrap `splice_sim` into a readily executable start-to-end framework for simulating and evaluating complex experimental scenarios with mapping based processing pipelines.

Contents
========

 * [Features](#features)
 * [Installation](#installation)
 * [Usage](#usage)
 * [Configuration](#configuration)
 * [Output Structure](#output-structure)

Features
========

<img src="img/splice_sim-graphical_abstract.png" width="60%" class="center">

+ Realistic Illumina short-read simulation using [ART](https://doi.org/10.1093/bioinformatics/btr708)
+ Simulation of customizable nucleotide-conversions at configurable rates
+ Simulation of isoforms at configurable splicing states per transcript
+ Mapping accuracies and performance metrics at different scopes of genomic annotation
+ Elaborate output tracks for visual inspection stratified by performance metric

Installation
============

`splice_sim` itself is a python package with several dependencies that are best simply installed in a [conda environment](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html):

```bash
conda env create -f environment.yml
```

One additional python package called `genomic_iterators` has to be installed manually into the environment:

> **Warning**
> `genomic_iterators` is not yet public, so you will have to first contact [Niko Popitsch](mailto:niko.popitsch@univie.ac.at) to obtain permissions to the repository!

```bash
source activate splice_sim
pip install git+https://github.com/popitsch/genomic_iterators.git
```

Then clone the `splice_sim` repository and now you are able to call it via the `main.py`:

```bash
git clone https://github.com/popitsch/splice_sim.git
cd splice_sim
python main.py
```

To run our full blown all-on-one Nextflow based workflow, you simply need to install [Nextflow](https://www.nextflow.io/) and [Docker](https://docs.docker.com/get-docker/) to have all depenencies available and be ready to go. To run `splice_sim` on HPC environments, most administrators prefer [Singularity](https://docs.sylabs.io/guides/3.5/user-guide/introduction.html) containers which can be seamlessly created from our [Docker container](https://hub.docker.com/repository/docker/tobneu/splice_sim).

Usage
=====

`splice_sim` itself is the Python package and includes the simulation engine to both simulate and evaluate datasets. `splice_sim` provides dedicate commands for the individual logical steps starting from creating transcript models, simulating reads from those transcripts and finally evaluating the performance of a given mapper via the produced bam file. The relevant steps for such a process are described [below](#splice_sim-engine). We have also wrapped a ready-to-use out of the box workflow that executes all step from start to end into our Nextflow workflow with a description provided in the [adjacent section](#nextflow-workflow).

## `splice_sim` engine

### build_model

The `build_model` command takes the reference and configuration provided by the user and creates the transcript model and sequence files needed that contains the composition of the transcriptome and serves and input to the read simulation step of `splice_sim`.

```shell
 python splice_sim/main.py build_model --help
usage: main.py [-h] -c config_file [-o outdir]

  Copyright (C) 2021 XXX.  All rights reserved.

  Distributed on an "AS IS" basis without warranties
  or conditions of any kind, either express or implied.

USAGE

optional arguments:
  -h, --help            show this help message and exit
  -c config_file, --config config_file
                        JSON config file
  -o outdir, --outdir outdir
                        output directory (default is current dir)
```
### create_genome_bam

The `create_genome_bam` command takes the simulated read set with your short-read simulator of choice (we use [ART](https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm)) and calculates the truth alignments that serve as the reference benchmark in `splice_sim`.

```shell
python splice_sim/main.py create_genome_bam --help
usage: main.py [-h] -m model_file -a config_file [-t threads] [-o outdir]

  Copyright (C) 2021 XXX.  All rights reserved.

  Distributed on an "AS IS" basis without warranties
  or conditions of any kind, either express or implied.

USAGE

optional arguments:
  -h, --help            show this help message and exit
  -m model_file, --model model_file
                        model file
  -a config_file, --art_sam config_file
                        ART sam file
  -t threads, --threads threads
                        threads
  -o outdir, --outdir outdir
                        output directory (default is current dir)
```

### postfilter_bam

The `postfilter_bam` command filters secondary and supplementary reads and highlights isoforms.

```bash
python splice_sim/main.py postfilter_bam --help
usage: main.py [-h] -c config_file -b bam_file -o outdir

  Copyright (C) 2021 XXX.  All rights reserved.

  Distributed on an "AS IS" basis without warranties
  or conditions of any kind, either express or implied.

USAGE

optional arguments:
  -h, --help            show this help message and exit
  -c config_file, --config config_file
                        JSON config file
  -b bam_file, --bam bam_file
                        input bam file
  -o outdir, --outdir outdir
                        output dir
```

### evaluate

The `evaluate` command runs the `splice_sim` evaluation routine on a given mapper bam file produced for a condition of the `splice_sim` simulation run.

```bash
python splice_sim/main.py evaluate --help
usage: main.py [-h] -b bam_file -c config_file -m model_file [-f filter_bed] [-t THREADS] -o outdir

  Copyright (C) 2021 XXX.  All rights reserved.

  Distributed on an "AS IS" basis without warranties
  or conditions of any kind, either express or implied.

USAGE

optional arguments:
  -h, --help            show this help message and exit
  -b bam_file, --bam_file bam_file
                        input bam file
  -c config_file, --config config_file
                        JSON config file
  -m model_file, --model model_file
                        model file
  -f filter_bed, --filter_bed filter_bed
                        Filter regions BED file
  -t THREADS, --threads THREADS
                        used threads
  -o outdir, --outdir outdir
                        output dir
```

### extract_feature_metadata

The `extract_feature_metadata` extracts comprehensive metadata that lists various characteristics of the genomic features under evaluation.

```bash
python splice_sim/main.py extract_feature_metadata --help
usage: main.py [-h] -c config_file -m model_file -o outdir

  Copyright (C) 2021 XXX.  All rights reserved.

  Distributed on an "AS IS" basis without warranties
  or conditions of any kind, either express or implied.

USAGE

optional arguments:
  -h, --help            show this help message and exit
  -c config_file, --config config_file
                        JSON config file
  -m model_file, --model model_file
                        model file
  -o outdir, --outdir outdir
                        output dir
```

# Processing results into R objects (RDS)

We provide a Rscript that takes all `splice_sim` evaluation outputs and processes them into readily usable RDS objects to be imported in R. The script is located in `splice_sim/src/main/R/splice_sim/preprocess_results.R` and all dependencies are again wrapped into our `splice_sim` [R Docker container](https://hub.docker.com/repository/docker/tobneu/splice_sim_r) `tobneu/splice_sim_r:latest`.

```bash
Rscript --vanilla splice_sim/src/main/R/splice_sim/preprocess_results.R

Error: usage: preprocess_results.R <splice_sim_config> [<outdir>]
Execution halted

```

# Nextflow workflow

Configuration
=============

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

Output Structure
===============

## Mapper count tables

These tables contain the performance metrics for a given mapper in the `count/*.counts.tsv.gz` files.

| Column            | Description                                                                                                                                               | Notes    |
|-------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------|------------|
| `mapper`          | Name of the respective mapper                                                                                                                             |          |
| `conversion_rate` | Conversion rate between 0 and 1.0                                                                                                                         |          |
| `fid`             | feature/transcript id                                                                                                                                     |          |
| `true_isoform`    | name of the isoform (as configured) the read originates from                                                                                               |          |
| `cv1`             | 1 if read contains at least one NC or 0 otherwise                                                                                                         |          |
| `cv2`             | 1 if read contains at least two NC or 0 otherwise                                                                                                         |          |
| `se1`             | 1 if read contains at least one simulated sequencing error or 0 otherwise                                                                                 |          |
| `se2`             | 1 if read contains at least two simulated sequencing errors or 0 otherwise                                                                                 |          |
| `classification`  | Read classification: TP: true positive, FN: false negative: FP_raw: false-positive/not  <br /> normalised, FP: false-positive/normalised |          |
| `count`           | read count. FP classified rows may include fractions                                                                                                       |          |
| `class_type`      | read type: acc: acceptor spanning, don: donor spanning, spl: spliced read                                                                                 |  SJ only |


## Metadata tables

These tables contain various metadata for the genomic intervals under investigation stratified at transcript level (`tx`), exon / intron feature level (`fx`) or splice-junction level (`sj`) in the `meta/*.metadata.tsv.gz` files.

| Column                 | Description                                                                                                                                  | Notes       |
|------------------------|---------------------------------------------------------------------------------------------------------------------------------------------|-------------|
| `tid`                  | Transcript ID                                                                                                                                |             |
| `fid`                  | Feature ID (intron or exon ID)                                                                                                              | FX+SJ  only |
| `ftype`                | Feature type: `tx`, `fx`, `don`, `acc` or `spl`                                                                                              |             |
| `rnk`                  | Rank. For transcripts this is the  number exons, for introns/exons it is the rank from the transcript 5'-end                                |             |
| `chromosome`           | Chromosome of the annotate feature                                                                                                          |             |
| `start` / `end`        | Genomic start/end position of the annotated feature                                                                                          |             |
| `strand`               | Strand of the annotation                                                                                                                    |             |
| `A`/`C`/`T`/`G`        | Number of A/C/T/G bases in the annotated sequence                                                                                            |             |
| `mean_map`             | Mean mappability for the annotated feature. Calculated from the configured mappability bedgraph file                                        |             |
| `tx_rnk`               | Rank in the transcript                                                                                                                      | FX+SJ only  |
| `num_exons`            | Number of exons in tx; 1,2,3,4,5,>5                                                                                                          |             |
| `tx_mappability`       | Transcript mappability, factor with levels: low, medium, high                                                                                |  FX+SJ only |
| `len`                  | Length of annotated feature                                                                                                                  |             |
| `mappability`          | Annotation mappability, factor with levels: low, medium, high                                                                                |             |
| `GC`                   | Fraction of G/C for annotated feature                                                                                                        |             |
| `frac_convertible`     | Fraction of convertible bases for annotation                                                                                                |  SJ only    |
| `convertibility`       | Convertibility, factor with levels: low, medium, high                                                                                        |             |
| `don_ex_A`/`C`/`T`/`G` | Number of A/C/T/G bases in exonic part of donor window <br /> (genomic window centred on splice donor site with size: 2xreadlen+1) |  SJ only    |
| `don_in_A`/`C`/`T`/`G` | Number of A/C/T/G bases in intronic part of donor window                                                                                    |  SJ only    |
| `don_win_map`          | Mean mappability of donor window                                                                                                            |  SJ only   |
| `don_mappability`      | Donor window mappability, factor w levels: low, medium, high                                                                                |  SJ only     |
| `don_ex_fc`            | Fraction of convertible bases in the exonic part of the donor window                                                                        |  SJ only    |
| `don_in_fc`            | Fraction of convertible bases in the intronic part of the donor window                                                                      |  SJ only    |
| `ac_*`                 | Analogous to the splice donor columns above, but for splice acceptor site                                                                    |  SJ only    |

