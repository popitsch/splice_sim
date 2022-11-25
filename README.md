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

#### Mapper count tables

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


#### Metadata tables

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

