FROM debian:bullseye-slim AS builder

LABEL version="0.0.1"
LABEL maintainer="HA"
LABEL description="This is a base image for the methyldackel tool."

RUN apt-get update && apt-get install -y --no-install-recommends \
    wget \
    bzip2 \
    ca-certificates \
    && wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh \
    && bash /tmp/miniconda.sh -b -p /miniconda3 \
    && rm /tmp/miniconda.sh \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

ENV PATH="/miniconda3/bin:${PATH}"

RUN conda install -c bioconda methyldackel \
    && conda install -c conda-forge -c bioconda bsmap \
    && conda clean -afy

FROM debian:bullseye-slim

LABEL version="0.0.1"
LABEL maintainer="HA"
LABEL description="This is a base image for the methyldackel tool."

COPY --from=builder /miniconda3 /miniconda3

ENV PATH="/miniconda3/bin:${PATH}"
ENV LANG=en_US.utf-8