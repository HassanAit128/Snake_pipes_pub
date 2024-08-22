#!/bin/bash
#SBATCH --cpus-per-task 48
#SBATCH --mem 100GB
#SBATCH --partition=long

module purge
module load snakemake/8.0.1

snakemake -s "cut&tag_part_2.smk" --cores 48 --configfile config_atac_seq.yaml --use-singularity