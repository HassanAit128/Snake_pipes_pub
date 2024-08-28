### Base Image
FROM r-base:4.2.3

LABEL version="0.0.1"
LABEL maintainer="H A"
LABEL description="This is a base image for the CUT&TaG/ATAC-seq pipeline."

LABEL maintainer="HA"

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
    && apt-get clean

RUN R --slave -e 'install.packages("tidyverse", repos = "https://cran.irsn.fr/")' \
    && R --slave -e 'install.packages("BiocManager", repos = "https://cran.irsn.fr/")' \
    && R --slave -e 'BiocManager::install(version = "3.16")' \
    && R --slave -e 'BiocManager::install("Rsamtools")' \
    && R -e "BiocManager::install('DiffBind')"