## -----------------------------------------
## Snakemake workflow for WGBS data analysis
## -----------------------------------------
## ------ Part 1 - FastQC and MultiQC ------
## -----------------------------------------
## Author: Hassan AitOu
## Last updated: 2024-08-27
## Purpose: This workflow is used to perform FastQC on the raw data and then generate a MultiQC report.
## Requirements: The raw data should be in fastq format (compressed or not) and the path to the raw data should be provided in the config file.
## -----------------------------------------
## Rules:
## ------

singularity: "docker://hasba/main_d"

import glob
import os
import time
from Scripts.shared_functions import get_extension

start = time.time()

############################################# Config ##########################################################
RUN_NAME = config['name'] # Name of the run
WORKDIR = config['working_dir'] # Directory where the results will be stored
DATADIR = config['raw_data_dir'] # Directory containing the data
RESULTSDIR = WORKDIR + '/' + RUN_NAME # Directory where the results will be stored
patterns = ['/*.fq', '/*.fq.gz', '/*.fastq', '/*.fastq.gz']
FILE_LIST = []
for pattern in patterns:
    FILE_LIST.extend(glob.glob(DATADIR + pattern))
SAMPLES = [os.path.splitext(os.path.splitext(os.path.basename(x))[0])[0] for x in FILE_LIST] # List of samples
##############################################################################################################

############################################# Error Handling #################################################
if WORKDIR == "":
    raise ValueError("Please provide a working directory in the config file")
if not os.path.exists(WORKDIR):
    raise ValueError("The working directory does not exist")
if DATADIR == "":
    raise ValueError("Please provide a path to the raw data in the config file")
if not os.path.exists(DATADIR):
    raise ValueError("The path to the raw data does not exist")
if RUN_NAME == "":
    raise ValueError("Please provide a name for the project/run in the config file")
if FILE_LIST == []:
    raise ValueError("No fastq files found in the directory, please check the path to the raw data")

onsuccess:
    if not snakemake.rule == "help":
        time_to_finish = round((time.time() - start)/60,1)
        longest_line_length = max(len(f"OUTPUT folder is found at: {RESULTSDIR}"), len("Running time in minutes: %s " % time_to_finish))
        total_length = longest_line_length + 8
        hash_count = (total_length - len("Workflow finished, no error") - 4) // 2
        print("\n" + "#" * hash_count + "# " + "WORKFLOW FINISHED, NO ERROR" + " #" + "#" * hash_count)
        print(f"OUTPUT folder is found at: {RESULTSDIR}")
        print(f"Running time in minutes: {time_to_finish}")
        print("#" * total_length + "\n\n")

onerror:
    print("\n\n###################### An error occurred ######################\n")
    print("Running time in minutes: %s\n" % round((time.time() - start)/60,1))
    print("\n###############################################################\n\n")
##############################################################################################################

############################################# Workflow #######################################################
## all              : rule to run all the rules and build all final outputs
rule all:
    input:
        f"{RESULTSDIR}/multiqc_report.html"

## fastqc           : rule to run FastQC on the raw data
rule fastqc:
    input:
        lambda wildcards: f"{DATADIR}/{wildcards.sample}"+ f"{get_extension(FILE_LIST, wildcards.sample)}",
    output:
        directory("{run_dir}/FASTQC_REPORTS/{sample}")
    params:
        fastqc_args = config['fastqc_args'],
        run_dir = RESULTSDIR
    shell:
        """
        mkdir -p {params.run_dir}/FASTQC_REPORTS/{wildcards.sample}
        fastqc {input} -o {output} {params.fastqc_args}
        """

## multiqc          : rule to generate a MultiQC report
rule multiqc:
    input:
        expand("{run_di}/FASTQC_REPORTS/{sample}", sample=SAMPLES, run_di=RESULTSDIR)
    output:
        "{run_dir}/multiqc_report.html"
    params:
        run_dir = RESULTSDIR
    shell:
        """
        multiqc {params.run_dir}/FASTQC_REPORTS/ -n {output}
        """

## help             : rule to display help message
rule help:
    input: "WGBS_part1_v0.smk"
    shell:
        """
        echo -e "\\033[0;32m"
        sed -n 's/^##\\([^#]\\)/\\1/p' {input}
        echo -e "\\033[0;33m"
        echo "Example Usage: snakemake -s WGBS_part1_v0.smk --configfile Config/config.yaml --cores 4 [--use-singularity]"
        echo ""
        echo "Please refer to the README file for more information on how to use this workflow with the parameters."
        echo -e "\\033[0m"
        """
## -----------------------------------------
##############################################################################################################