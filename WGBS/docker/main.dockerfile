# Builder stage
FROM ubuntu:20.04 as builder

LABEL version="0.0.1"
LABEL maintainer="H A"
LABEL description="This is a base image for most of the WGBS pipeline."

RUN echo 'debconf debconf/frontend select Noninteractive' | debconf-set-selections \
    && apt-get update && apt-get install -y \
    python3 \
    python3-pip \
    python3-dev \
    python3-setuptools \
    apt-utils \
    wget \
    bash \
    gcc \
    make \
    zlib1g-dev \
    locales \
    tcl \
    bzip2 \
    gzip \
    build-essential \
    default-jre \
    cutadapt \
    bedtools \
    samtools \
    bowtie2 \
    trim-galore \
    fastqc \
    perl \
    && wget https://codeload.github.com/FelixKrueger/Bismark/tar.gz/refs/tags/v0.24.2 -O /tmp/bismark.tar.gz \
    && tar -xvf /tmp/bismark.tar.gz -C /tmp \
    && mv /tmp/Bismark-0.24.2 /Bismark-0.24.2 \
    && rm /tmp/bismark.tar.gz \
    && wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh \
    && bash /tmp/miniconda.sh -b -p /miniconda3 \
    && rm /tmp/miniconda.sh \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*
ENV PATH="/miniconda3/bin:${PATH}"
ENV PATH="/Bismark-0.24.2:${PATH}"
RUN /miniconda3/bin/conda install -c bioconda -c conda-forge picard=2.23.8 \
    && conda install -c bioconda methyldackel \
    && conda install -c conda-forge -c bioconda bsmap \
    && conda install -y -c bioconda -c conda-forge multiqc \
    && /miniconda3/bin/conda clean -afy

# Final stage
FROM ubuntu:20.04

LABEL version="0.0.1"
LABEL maintainer="H A"
LABEL description="This is a base image for most of the WGBS pipeline."

COPY --from=builder /miniconda3 /miniconda3
COPY --from=builder /Bismark-0.24.2 /Bismark-0.24.2

RUN apt-get update && apt-get install -y --no-install-recommends \
    python3 \
    python3-pip \
    bash \
    locales \
    samtools \
    bedtools \
    bowtie2 \
    trim-galore \
    fastqc \
    perl \
    && pip3 install pandas openpyxl matplotlib \
    && pip install numpy \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

ENV PATH="/Bismark-0.24.2:${PATH}"
ENV PATH="/miniconda3/bin:${PATH}"
ENV LC_ALL=en_US.UTF-8
ENV LANG=en_US.UTF-8