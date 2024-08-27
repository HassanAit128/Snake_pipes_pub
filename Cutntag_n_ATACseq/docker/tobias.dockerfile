FROM python:3.8-slim

LABEL version="0.0.1"
LABEL maintainer="H A"
LABEL description="This is a base image with TOBIAS environment."

RUN apt-get update && apt-get install -y \
    build-essential \
    python3-pip \
    r-base \
    bedtools \
    samtools \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

    RUN pip install --no-cache-dir \
    tobias uropa \
    && R --slave -e 'install.packages(c("ggplot2","devtools","gplots","gridExtra","jsonlite", "VennDiagram","snow","getopt","tidyr","UpSetR"))'

RUN pip install --no-cache-dir \
    macs2 \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*
