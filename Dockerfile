FROM continuumio/miniconda3:4.10.3p0

MAINTAINER Tobias Neumann <tobias.neumann.at@gmail.com>

COPY environment.yml /tmp/environment.yml

RUN apt-get update \
    && apt-get install -y libz-dev gcc g++ make unzip \
    && apt-get clean -y \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir -p /splice_sim \
    && cd /splice_sim \
    && conda env create -f /tmp/environment.yml \
    && /opt/conda/envs/splice_sim/bin/pip install git+https://github.com/popitsch/genomic_iterators.git \
    && wget https://github.com/alexdobin/STAR/archive/refs/tags/2.7.1a.tar.gz \
    && tar xvfz 2.7.1a.tar.gz \
    && cd STAR-2.7.1a/source \
    && make STAR \
    && mv STAR STAR-2.7.1a \
    && git clone https://github.com/DaehwanKimLab/hisat2.git /splice_sim/hisat-3n \
    && cd /splice_sim/hisat-3n \
    && git checkout -b hisat-3n origin/hisat-3n \
    && make \
    && wget https://icbi.i-med.ac.at/software/meRanTK/downloads/1.2.1b/meRanTK-1.2.1b.zip -P /splice_sim \
    && cd /splice_sim \
    && unzip meRanTK-1.2.1b.zip \
    && rm -rf /opt/conda/pkgs/*

ENV PATH /splice_sim/STAR-2.7.1a/source:/splice_sim/meRanTK-1.2.1b:/splice_sim/hisat-3n:/opt/conda/envs/splice_sim/bin:$PATH
