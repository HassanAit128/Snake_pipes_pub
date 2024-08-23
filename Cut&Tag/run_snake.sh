#!/bin/bash
#SBATCH --cpus-per-task 48
#SBATCH --mem 100GB
#SBATCH --partition=long

module purge
module load snakemake/8.0.1
module load samtools/1.15.1
module load picard/2.23.5
module load bowtie2/2.5.1
module load deeptools/3.5.4
module load macs2/2.2.7.1
module load multiqc/1.13
module load fastqc/0.12.1
module load cutadapt/4.5
module load pigz/2.3.4
module load trim-galore/0.6.10

snakemake -s "cut&tag_part_3.smk" --cores 48 --configfile config_atac_seq.yaml