FROM bioconductor/bioconductor_docker

RUN apt-get update \
	&& apt-get clean \
	&& rm -rf /var/lib/apt/lists/*

RUN R --slave -e 'install.packages("tidyverse", repos = "https://cran.irsn.fr/")' \
    && R --slave -e 'BiocManager::install("methylKit")'

RUN R --slave -e 'BiocManager::install("genomation")' \
    && R --slave -e 'BiocManager::install("GenomicRanges")' \
    && rm -rf /tmp/* /var/tmp/*