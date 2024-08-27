FROM ubuntu:20.04

LABEL version="0.0.1"
LABEL maintainer="H A"
LABEL description="This is a base image for the CUT&TaG/ATAC-seq pipeline."

RUN echo 'debconf debconf/frontend select Noninteractive' | debconf-set-selections \
    && apt-get update && apt-get install -y \
    python3 \
    python3-pip \
    python3-dev \
    python3-setuptools \
    apt-utils \
    wget \
    bash \
    curl \
    gcc \
    make \
    zlib1g-dev \
    locales \
    git \
    tcl \
    bzip2 \
    gzip \
    build-essential \
    default-jre \
    cutadapt \
    bedtools \
    samtools \
    bowtie2 \
    fastqc \
    trim-galore \
    r-base \
    && apt-get clean

# Install Python packages
RUN pip3 install --no-cache-dir "cmake==3.26.1" \
    && pip install pandas openpyxl cykhash Cython MACS3 \
    && pip install numpy

# Install conda packages
RUN apt-get install -y wget
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
RUN bash Miniconda3-latest-Linux-x86_64.sh -b -p /miniconda3
RUN rm Miniconda3-latest-Linux-x86_64.sh
ENV PATH="/miniconda3/bin:${PATH}"
RUN conda config --add channels defaults
RUN conda config --add channels bioconda
RUN conda config --add channels conda-forge

RUN conda install --solver=classic conda-forge::conda-libmamba-solver conda-forge::libmamba conda-forge::libmambapy conda-forge::libarchive
RUN conda install -y -c bioconda -c conda-forge \
    picard=2.23.8
RUN conda install -y -c bioconda -c conda-forge \
    deeptools
RUN conda install -y -c bioconda -c conda-forge \
    multiqc    

# Install R packages
RUN R --slave -e 'install.packages("tidyverse", repos = "https://cran.irsn.fr/")' \
    && R --slave -e 'install.packages("BiocManager", repos = "https://cran.irsn.fr/")' \
    && R --slave -e 'BiocManager::install("DiffBind")' \
    && R --slave -e 'BiocManager::install("Rsamtools")'

ENV PATH="/miniconda3/bin:${PATH}"