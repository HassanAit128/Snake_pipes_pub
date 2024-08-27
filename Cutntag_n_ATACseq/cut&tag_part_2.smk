singularity: "docker://hasba/atac_cut:latest"

import glob
import pandas as pd
import os
import csv
import time
import re
from Scripts.shared_functions import check_if_design_matrix_is_valid, get_files_suffix

start = time.time()


############################################# Config ##########################################################
WORKDIR = config['working_directory'] # Directory where the pipeline will be run
DATADIR = config['path_to_raw_data'] # Directory containing the data
RESULTSDIR = WORKDIR + '/' + config['run_name'] # Directory where the results will be stored
patterns = ['*.fq', '*.fq.gz', '*.fastq', '*.fastq.gz']
FILE_LIST = []
for pattern in patterns:
    FILE_LIST.extend(glob.glob(os.path.join(DATADIR, '*_R1*' + pattern)))
pattern = re.compile(r'(.+)_R1.*\.(fq|fq\.gz|fastq|fastq\.gz)$')
SAMPLES = [pattern.match(os.path.basename(x)).group(1) for x in FILE_LIST if pattern.match(os.path.basename(x))]
DESIGN_MATRIX = config['design_matrix']  # Read the design matrix file into a DataFrame
RUN_NAME = config['run_name'] # Name of the run
RMDUP = config['remove_duplicates'] # Remove duplicates True/False
KEEP_UNIQ = True
KEEP_REGULAR = config['keep_only_regular_reads'] # Keep only regular reads True/False
MAPQ = config['MAPQ'] # MAPQ value to use when filtering reads
GENOME = config['genome'] # Genome used for the analysis
IDX_GNM = config["index_genome"] # Whether to index the genome or not
GENOME_PATH = config["path_to_genome_fasta"] # Path to the genome to be indexed in fasta format
GNM_IDX = config["path_to_bowtie2_genome_idx"] # Path to the genome index in bowtie2 format
DUPLI = "without" if RMDUP else "with" # Duplicates or not
ALSO_WITH_DUPLICATES = config['output_finale_file_with_duplicates'] # Also run analysis with duplicates True/False
RMV_TEMP_FILES = config['delete_tmp_files'] # Remove temporary files True/False
ATAC_SEQ_MODE = config['atac_seq'] # ATAC-seq mode True/False
TRIM = config['trim_data'] # Trim and remove adapters True/False
##############################################################################################################

############################################# Error Handling #################################################
check_if_design_matrix_is_valid(DESIGN_MATRIX, SAMPLES)

if RUN_NAME == "":
    raise ValueError("Please provide a name for the project/run in the config file")
if WORKDIR == "":
    raise ValueError("Please provide a working directory in the config file")
if not os.path.exists(WORKDIR):
    raise ValueError("Working directory not found. Please check the path to the working directory.")
if DATADIR == "":
    raise ValueError("Please provide a path to the raw data in the config file")
if not os.path.exists(DATADIR):
    raise ValueError("Data directory not found. Please check the path to the raw data.")
if SAMPLES == []:
    raise ValueError("No samples detected in the data directory. Please check the raw data directory.")
if not os.path.exists(DESIGN_MATRIX):
    raise ValueError("Design matrix file not found. Please check the path to the design matrix file.")
if not os.path.exists(config["path_to_bowtie2_genome_idx"]):
    raise ValueError("Bowtie2 genome index not found. Please check the path to the genome index. Or set 'index_genome' to True to build the index.")
if type(MAPQ) == int:
    MAPQ = [MAPQ]
if GENOME not in ["mm10", "hg38"]:
    raise ValueError("Genome not provided, or, Entered genome not supported. Please choose either 'mm10' or 'hg38'.")
if IDX_GNM is True:
    if not os.path.exists(GENOME_PATH):
        raise ValueError("Genome fasta file not found. Please check the path to the genome fasta file.")
    if GENOME_PATH is None:
        raise ValueError("Genome needed to build the index is not provided. Please provide a genome.")  
    if GENOME not in ["mm10", "hg38"]:
        raise ValueError("Genome prefix not provided, or, Entered genome not supported. Please choose either 'mm10' or 'hg38'.")
if ALSO_WITH_DUPLICATES is True and RMDUP is False:
    ALSO_WITH_DUPLICATES = False

onsuccess:
    time_to_finish = round((time.time() - start)/60,1)
    longest_line_length = max(len(f"OUTPUT folder is found at: {RESULTSDIR}"), len("Running time in minutes: %s " % time_to_finish))
    total_length = longest_line_length + 8
    suffix = create_file_path_suffix(RMDUP, KEEP_UNIQ, KEEP_REGULAR)
    files_to_keep = []
    for q in MAPQ :
        rule_all_input_temp = expand("{run_di}/BAM/{dupli}_duplicates/mapq{mapq}/{sample}.mapped.s{suffix}.s.bam", sample=SAMPLES, run_di=RESULTSDIR, suffix=create_file_path_suffix(mapq=q, rmdup=RMDUP, uniq=KEEP_UNIQ, regular=KEEP_REGULAR), mapq=q, dupli=DUPLI)
        files_to_keep.extend(rule_all_input_temp)
        rule_all_input_temp = expand("{run_di}/BAM/{dupli}_duplicates/mapq{mapq}/{sample}.mapped.s{suffix}.s.bam.bai", sample=SAMPLES, run_di=RESULTSDIR, suffix=create_file_path_suffix(mapq=q, rmdup=RMDUP, uniq=KEEP_UNIQ, regular=KEEP_REGULAR), mapq=q, dupli=DUPLI)
        files_to_keep.extend(rule_all_input_temp)
        if ALSO_WITH_DUPLICATES and RMDUP:
            for q in MAPQ :
                rule_all_input_temp = expand("{run_di}/BAM/with_duplicates/mapq{mapq}/{sample}.mapped.s{suffix}.s.bam", sample=SAMPLES, run_di=RESULTSDIR, suffix=create_file_path_suffix(mapq=q, rmdup=False, uniq=KEEP_UNIQ, regular=KEEP_REGULAR), mapq=q)
                files_to_keep.extend(rule_all_input_temp)
                rule_all_input_temp = expand("{run_di}/BAM/with_duplicates/mapq{mapq}/{sample}.mapped.s{suffix}.s.bam.bai", sample=SAMPLES, run_di=RESULTSDIR, suffix=create_file_path_suffix(mapq=q, rmdup=False, uniq=KEEP_UNIQ, regular=KEEP_REGULAR), mapq=q)
                files_to_keep.extend(rule_all_input_temp)
    if RMV_TEMP_FILES:
        delete_tmp_files(files_to_keep)

    hash_count = (total_length - len("Workflow finished, no error") - 4) // 2
    print("\n" + "#" * hash_count + "# " + "WORKFLOW FINISHED, NO ERROR" + " #" + "#" * hash_count)
    print(f"OUTPUT folder is found at: {RESULTSDIR}")
    print(f"Running time in minutes: {time_to_finish}")
    print(f"Mode: {'ATAC-seq' if ATAC_SEQ_MODE else 'CUT&Tag'}")
    print(f"Genome: {GENOME}")
    print(f"Samples detected in the data directory: {SAMPLES}")
    print(f"Raw Fastq files were {'trimmed' if TRIM else 'not trimmed'}!")
    print(f"MAPQ value used: {MAPQ}. {'' if MAPQ != '' else '(No MAPQ value provided, default value 0 was used.)'}")
    print(f"Duplicates were {'removed' if RMDUP else 'kept'}! {'(Duplicate analysis was also performed)' if ALSO_WITH_DUPLICATES else ''}")
    print(f"Only regular reads were {'kept' if KEEP_REGULAR else 'kept'}!")
    print(f"Temporary files were {'removed' if RMV_TEMP_FILES else 'kept'}!")
    print("#" * total_length + "\n\n")

onerror:
    print("\n\n###################### An error occurred ######################\n")
    print("Running time in minutes: %s\n" % round((time.time() - start)/60,1))
    print("\n###############################################################\n\n")
##############################################################################################################

############################################ Helper Functions ################################################
def delete_tmp_files(files_to_keep):
    """
    Function to delete temporary files

    Pararmeters:
    ------------
        files_to_keep : list
            List of files to keep
    
    Returns:
    --------
        None
    """
    deleted_files = []
    if RMV_TEMP_FILES:
        all_files = glob.glob(f"{RESULTSDIR}/BAM/**/*.bam", recursive=True) + glob.glob(f"{RESULTSDIR}/BAM/**/*.bam.bai", recursive=True)
        for file in all_files:
            if file not in files_to_keep and not file.endswith("mapped.bam"):
                os.remove(file)
                deleted_files.append(file)

def create_file_path_suffix(mapq, rmdup=False, uniq=False, regular=False, suffix="", step=""):
    """
    Function to create a suffix for the file path based on the parameters
    
    Parameters:
    -----------
        mapq : int
            MAPQ value
        rmdup : bool
            Remove duplicates or not
        uniq : bool
            Keep only unique reads or not
        regular : bool
            Keep only regular reads or not
        suffix : str
            Suffix to add to the file path
        step : str
            Step of the pipeline where the function is called
    
    Returns:
    --------
        suffix : str
            Suffix for the file path
    """

    if rmdup or uniq or regular:
        suffix += "."
    if rmdup:
        suffix += "rmdup"
    if uniq:
        suffix += "." if rmdup else ""
        suffix += f"mapq{mapq}.filtered"
    if regular:
        if rmdup or uniq:
            suffix += "."
        suffix += "regularChr"
    return suffix
SUFFIX = create_file_path_suffix(RMDUP, KEEP_UNIQ, KEEP_REGULAR)

def tmp_get_echo_out(rmdup=False, uniq=False, regular=False):
    """
    Function to get the parameters output
    
    Parameters:
    -----------
        rmdup : bool
            Remove duplicates or not
        uniq : bool
            Keep only unique reads or not
        regular : bool
            Keep only regular reads or not
    
    Returns:
    --------
        out : str
            Parameters output
    """
    out = ""
    if uniq:
        out += "UNIQUE " if MAPQ != 0 else "MULTIPLE "
    if regular:
        out += "REGULAR "
    if rmdup:
        out += "WITHOUT duplicates"
    else:
        out += "WITH duplicates"
    return out    
##############################################################################################################

############################################## Output ########################################################
rule_all_input = []

if SUFFIX != "":
    if KEEP_UNIQ:
        for q in MAPQ :
            rule_all_input_temp = expand("{run_di}/BAM/{dupli}_duplicates/mapq{mapq}/{sample}.mapped.s{suffix}.s.bam", sample=SAMPLES, run_di=RESULTSDIR, suffix=create_file_path_suffix(mapq=q, rmdup=RMDUP, uniq=KEEP_UNIQ, regular=KEEP_REGULAR), mapq=q, dupli=DUPLI)     
            rule_all_input.extend(rule_all_input_temp)
else :
    if KEEP_UNIQ:
        for q in MAPQ :
            rule_all_input_temp = expand("{run_di}/BAM/{dupli}_duplicates/mapq{mapq}/{sample}.mapped.s.s.bam", sample=SAMPLES, run_di=RESULTSDIR, mapq=q, dupli="with")
            rule_all_input.extend(rule_all_input_temp)

multiqc_input = []
if TRIM: 
    for sample in SAMPLES:
        rule_all_input.append(f"{RESULTSDIR}/Trimmed_Fastq/{sample}_R1{get_files_suffix(DATADIR, sample, True)[0]}val_1.fq.gz")
        rule_all_input.append(f"{RESULTSDIR}/Trimmed_Fastq/{sample}_R2{get_files_suffix(DATADIR, sample, True)[0]}val_2.fq.gz")
        multiqc_input.append(f"{RESULTSDIR}/Trimmed_Fastq/fastqc_reports/{sample}_R1{get_files_suffix(DATADIR, sample, True)[0]}val_1_fastqc.html")
        multiqc_input.append(f"{RESULTSDIR}/Trimmed_Fastq/fastqc_reports/{sample}_R2{get_files_suffix(DATADIR, sample, True)[0]}val_2_fastqc.html")
    rule_all_input.extend(multiqc_input)
    rule_all_input.append(f"{RESULTSDIR}/Trimmed_Fastq/multiqc_report.html")

rule_all_input.append(f"{RESULTSDIR}/{RUN_NAME}_parameters_summary_part_2.txt")

summary_stats = expand("{run_dir}/stats/Summary/mapq{q}_{dupli}_duplicates.xlsx", q=MAPQ, run_dir=RESULTSDIR, dupli=DUPLI)
rule_all_input.extend(summary_stats)

if IDX_GNM:
    rule_all_input.append(f"{RESULTSDIR}/bt2_idx/{GENOME}/{GENOME}.1.bt2")
    rule_all_input.append(f"{RESULTSDIR}/bt2_idx/{GENOME}/{GENOME}.2.bt2")
    rule_all_input.append(f"{RESULTSDIR}/bt2_idx/{GENOME}/{GENOME}.3.bt2")
    rule_all_input.append(f"{RESULTSDIR}/bt2_idx/{GENOME}/{GENOME}.4.bt2")
    rule_all_input.append(f"{RESULTSDIR}/bt2_idx/{GENOME}/{GENOME}.rev.1.bt2")
    rule_all_input.append(f"{RESULTSDIR}/bt2_idx/{GENOME}/{GENOME}.rev.2.bt2")
##############################################################################################################

############################################## Rules #########################################################
rule all:
    input:
        rule_all_input

rule trim_galore:
    input: 
        R1 = lambda wildcards: f"{DATADIR}/{wildcards.sample}_R1{get_files_suffix(DATADIR, wildcards.sample)[0]}" + f"{get_files_suffix(DATADIR, wildcards.sample)[1]}",
        R2 = lambda wildcards: f"{DATADIR}/{wildcards.sample}_R2{get_files_suffix(DATADIR, wildcards.sample)[0]}" + f"{get_files_suffix(DATADIR, wildcards.sample)[1]}"
    output:
        R1 = "{run_dir}/Trimmed_Fastq/{sample}_R1{suffix}val_1.fq.gz",
        R2 = "{run_dir}/Trimmed_Fastq/{sample}_R2{suffix}val_2.fq.gz",
        fasqc_R1 = "{run_dir}/Trimmed_Fastq/fastqc_reports/{sample}_R1{suffix}val_1_fastqc.html",
        fasqc_R2 = "{run_dir}/Trimmed_Fastq/fastqc_reports/{sample}_R2{suffix}val_2_fastqc.html"
    params:
        run_dir = RESULTSDIR,
        adapter = config['adapter'],
        suffix = lambda wildcards: get_files_suffix(DATADIR, wildcards.sample, True)[0]
    threads: 4 if config['nbCores_trim_galore'] > 4 else config['nbCores_trim_galore']
    shell:
        """
        mkdir -p {params.run_dir}/Trimmed_Fastq/
        mkdir -p {params.run_dir}/Trimmed_Fastq/fastqc_reports/
        if [{params.adapter} == ""]
        then
            trim_galore --paired --cores {threads} --output_dir {params.run_dir}/Trimmed_Fastq/ --fastqc --fastqc_args "--outdir {params.run_dir}/Trimmed_Fastq/fastqc_reports/" {input.R1} {input.R2}
        else
            trim_galore --paired --cores {threads} --gzip --output_dir {params.run_dir}/Trimmed_Fastq/ --fastqc --fastqc_args "--outdir {params.run_dir}/Trimmed_Fastq/fastqc_reports/" --{params.adapter} {input.R1} {input.R2}
        fi
        """

rule multiqc:
    input:
        multiqc_input  
    output:
        "{run_dir}/Trimmed_Fastq/multiqc_report.html"
    params:
        run_dir = RESULTSDIR
    shell:
        """
        multiqc {params.run_dir}/Trimmed_Fastq/fastqc_reports/ -n {output}
        """

rule index_genome:
    input:
        genome = GENOME_PATH
    output:
        multiext(
            f"{RESULTSDIR}/bt2_idx/{GENOME}/{GENOME}",
            ".1.bt2",
            ".2.bt2",
            ".3.bt2",
            ".4.bt2",
            ".rev.1.bt2",
            ".rev.2.bt2",
        )
    params:
        run_dir = RESULTSDIR,
        genome_prefix = GENOME,
    threads: config['nbCores_bowtie2']
    shell:
        """
        mkdir -p {params.run_dir}/bt2_idx/{params.genome_prefix}
        bowtie2-build --threads {threads} {input.genome} {params.run_dir}/bt2_idx/{params.genome_prefix}/{params.genome_prefix}
        """

rule run_bowtie2:
    input:
        R1 = lambda wildcards: f"{DATADIR}/{wildcards.sample}_R1{get_files_suffix(DATADIR, wildcards.sample)[0]}" + f"{get_files_suffix(DATADIR, wildcards.sample)[1]}" if not TRIM else f"{RESULTSDIR}/Trimmed_Fastq/{wildcards.sample}_R1{get_files_suffix(DATADIR, wildcards.sample, True)[0]}val_1.fq.gz",
        R2 = lambda wildcards: f"{DATADIR}/{wildcards.sample}_R2{get_files_suffix(DATADIR, wildcards.sample)[0]}" + f"{get_files_suffix(DATADIR, wildcards.sample)[1]}" if not TRIM else f"{RESULTSDIR}/Trimmed_Fastq/{wildcards.sample}_R2{get_files_suffix(DATADIR, wildcards.sample, True)[0]}val_2.fq.gz",
        genomeIndex = GNM_IDX if not IDX_GNM else [f"{RESULTSDIR}/bt2_idx/{GENOME}/{GENOME}.1.bt2", f"{RESULTSDIR}/bt2_idx/{GENOME}/{GENOME}.2.bt2", f"{RESULTSDIR}/bt2_idx/{GENOME}/{GENOME}.3.bt2", f"{RESULTSDIR}/bt2_idx/{GENOME}/{GENOME}.4.bt2", f"{RESULTSDIR}/bt2_idx/{GENOME}/{GENOME}.rev.1.bt2", f"{RESULTSDIR}/bt2_idx/{GENOME}/{GENOME}.rev.2.bt2"]
    output:
        temp("{run_dir}/BAM/{sample}.bef.sam")
    log: 
        "{run_dir}/stats/{sample}_bowtie2Reports.txt"
    params:
        genome = config['genome'],
        genome_prefix = GENOME,
        run_dir = RESULTSDIR,
        txt = "{run_dir}/stats/{sample}_stats.txt",
        mode = "--end-to-end --very-sensitive --no-mixed --no-discordant --phred33 -X 2000" if ATAC_SEQ_MODE else "--end-to-end --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700"
    threads: config['nbCores_bowtie2']
    shell:
        """
        mkdir -p {params.run_dir}
        mkdir -p {params.run_dir}/BAM
        mkdir -p {params.run_dir}/stats
        mkdir -p {params.run_dir}/stats/Summary
        bowtie2 {params.mode} \
            -p {threads} -x "{input.genomeIndex}/{params.genome}" \
            -S {output} \
            -q \
            -1 {input.R1} \
            -2 {input.R2} &> {log}
        echo "Stats for : {wildcards.sample}" >> {params.txt}
        echo "" >> {params.txt}
        echo "### Bowtie 2 ###" >> {params.txt}
        cat {params.run_dir}/stats/{wildcards.sample}_bowtie2Reports.txt >> {params.txt}
        """

rule sam_to_bam:
    input:
        "{run_dir}/BAM/{sample}.bef.sam"
    output:
        bam = protected("{run_dir}/BAM/{sample}.aft.bam"),
    threads: config['nbCores_samtools']
    params:
        run_dir = RESULTSDIR
    log:
        "{run_dir}/LOGS/samtools/{sample}_sam_to_bam.txt"
    shell:
        """
        samtools view --threads {threads} -bo {output.bam} {input} 2> {log}
        read_count=$(samtools view -c {output.bam} --threads {threads})
        echo "" >> {params.run_dir}/stats/{wildcards.sample}_stats.txt
        echo "### Proccessing 1 ###" >> {params.run_dir}/stats/{wildcards.sample}_stats.txt
        echo "Number of reads after conversion from sam to bam (_aft.bam): $read_count" >> {params.run_dir}/stats/{wildcards.sample}_stats.txt
        """
            
rule map_bam:
    input:
        "{run_dir}/BAM/{sample}.aft.bam"
    output:
        bam = "{run_dir}/BAM/{sample}.mapped.bam"
    threads: config['nbCores_samtools']
    params:
        run_dir = RESULTSDIR
    log:
        "{run_dir}/LOGS/samtools/{sample}_map_bam.txt"
    shell:
        """
        samtools view --threads {threads} -b -F 4 {input} -o {output.bam} 2> {log}
        read_count=$(samtools view -F 4 -c {output.bam} --threads {threads})
        echo "Mapped reads after mapping bam using -F = 4 (_mapped.bam): $read_count" >> {params.run_dir}/stats/{wildcards.sample}_stats.txt
        """
        
rule sort_bam:
    input:
        "{run_dir}/BAM/{sample}.mapped.bam"
    output:
        bam = "{run_dir}/BAM/{sample}.mapped.s.bam"
    threads: config['nbCores_samtools']
    params:
        run_dir = RESULTSDIR
    log:
        "{run_dir}/LOGS/samtools/{sample}_sort_bam.txt"
    shell:
        """
        samtools sort --threads {threads} {input} -o {output.bam} 2> {log}
        read_count=$(samtools view -F 4 -c {output.bam} --threads {threads})
        echo "Sorted reads: $read_count" >> {params.run_dir}/stats/{wildcards.sample}_stats.txt
        """

rule remove_duplicates:
    input:
        "{run_dir}/BAM/{sample}.mapped.s.bam"
    output:
        r1 = "{run_dir}/BAM/{dupli}_duplicates/{sample}.mapped.s.rmdup.bam",
    threads: config['nbCores_picard']
    params:
        run_dir = RESULTSDIR,
        duplicates = DUPLI
    log:
        "{run_dir}/LOGS/picard/{dupli}_duplicates/{sample}_remove_duplicates.txt"
    shell:
        """
        picard -Xmx8g MarkDuplicates \
            -I {input} \
            -O {output.r1} \
            --REMOVE_DUPLICATES true \
            --METRICS_FILE {params.run_dir}/stats/{wildcards.sample}_mapped_s_picard_rmdup.txt \
            --VALIDATION_STRINGENCY LENIENT > {log} 2>&1
        mkdir -p {wildcards.run_dir}/stats
        read_count=$(samtools view -F 4 -c {output.r1} --threads {threads})
        cp {params.run_dir}/stats/{wildcards.sample}_stats.txt {params.run_dir}/stats/{wildcards.sample}_{params.duplicates}_duplicates_stats.txt
        echo "" >> {params.run_dir}/stats/{wildcards.sample}_{params.duplicates}_duplicates_stats.txt
        echo "### Picard ###" >> {params.run_dir}/stats/{wildcards.sample}_{params.duplicates}_duplicates_stats.txt
        awk '/## METRICS CLASS	picard.sam.DuplicationMetrics/{{flag=1;next}} /##/{{flag=0}} flag {{print}}' {params.run_dir}/stats/{wildcards.sample}_mapped_s_picard_rmdup.txt >> {params.run_dir}/stats/{wildcards.sample}_{params.duplicates}_duplicates_stats.txt
        echo "### Proccessing 2 ###" >> {params.run_dir}/stats/{wildcards.sample}_{params.duplicates}_duplicates_stats.txt
        echo "Reads WITHOUT duplicates (_mapped_s_rmdup.bam): $read_count" >> {params.run_dir}/stats/{wildcards.sample}_{params.duplicates}_duplicates_stats.txt
        """

rule filter_reads:
    input:
        "{run_dir}/BAM/{dupli}_duplicates/{sample}.mapped.s.rmdup.bam" if RMDUP else "{run_dir}/BAM/{sample}.mapped.s.bam"
    output:
        bam = "{run_dir}/BAM/{dupli}_duplicates/mapq{mapq}/{sample}.mapped.s{suffix}.mapq{mapq2}.filtered.bam" if RMDUP else "{run_dir}/BAM/{dupli}_duplicates/mapq{mapq}/{sample}.mapped.s.mapq{mapq}.filtered.bam"
    params:
        quality = config["MAPQ"],
        run_dir = RESULTSDIR,
        mapq2 = lambda wildcards: wildcards.mapq,
        suffix = lambda wildcards: create_file_path_suffix(mapq=wildcards.mapq, rmdup=RMDUP, step="filter_reads"),
        tem_echo_out = tmp_get_echo_out(RMDUP, True),
        also_with_duplicates = ALSO_WITH_DUPLICATES,
        duplicates = DUPLI
    threads: config['nbCores_samtools']
    log:
        "{run_dir}/LOGS/samtools/{dupli}_duplicates/mapq{mapq}/{sample}_{suffix}_{mapq2}_filter_reads.txt" if RMDUP else "{run_dir}/LOGS/samtools/{dupli}_duplicates/{sample}_{mapq}_filter_reads.txt"
    shell:
        """
        mkdir -p {params.run_dir}/stats/Summary
        samtools view --threads {threads} -q {wildcards.mapq} {input} -o {output.bam} 2> {log}
        read_count=$(samtools view -F 4 -c {output.bam} --threads {threads})
        if [ {params.duplicates} = without ]
        then
            cp {params.run_dir}/stats/{wildcards.sample}_{params.duplicates}_duplicates_stats.txt {params.run_dir}/stats/{wildcards.sample}_{params.duplicates}_duplicates_mapq{wildcards.mapq}_stats.txt
        else
            cp {params.run_dir}/stats/{wildcards.sample}_stats.txt {params.run_dir}/stats/{wildcards.sample}_with_duplicates_mapq{wildcards.mapq}_stats.txt
        fi
        echo "{params.tem_echo_out} reads with  Q>={wildcards.mapq} (_mapped_s{params.suffix}.bam): $read_count" >> {params.run_dir}/stats/{wildcards.sample}_{params.duplicates}_duplicates_mapq{wildcards.mapq}_stats.txt
        if [ {params.also_with_duplicates} = True ]
        then
            mkdir -p {params.run_dir}/BAM/with_duplicates
            mkdir -p {params.run_dir}/BAM/with_duplicates/mapq{wildcards.mapq}
            samtools view --threads {threads} -q {wildcards.mapq} {params.run_dir}/BAM/{wildcards.sample}.mapped.s.bam -o {params.run_dir}/BAM/with_duplicates/mapq{wildcards.mapq}/{wildcards.sample}.mapped.s.mapq{wildcards.mapq}.filtered.bam
            read_count=$(samtools view -F 4 -c {params.run_dir}/BAM/with_duplicates/mapq{wildcards.mapq}/{wildcards.sample}.mapped.s.mapq{wildcards.mapq}.filtered.bam --threads {threads})
            cp {params.run_dir}/stats/{wildcards.sample}_stats.txt {params.run_dir}/stats/{wildcards.sample}_with_duplicates_mapq{wildcards.mapq}_stats.txt
            if [ {wildcards.mapq} = 0 ]
            then
                multi_or_unique="MULTIPLE"
            else
                multi_or_unique="UNIQUE"
            fi
            echo "$multi_or_unique Reads WITH duplicates with Q>={wildcards.mapq} (_mapped_s.mapq{wildcards.mapq}.filtered.bam): $read_count" >> {params.run_dir}/stats/{wildcards.sample}_with_duplicates_mapq{wildcards.mapq}_stats.txt
        fi
        """     

rule keep_regular_chr:
    input:
        "{run_dir}/BAM/{dupli}_duplicates/mapq{mapq}/{sample}.mapped.s{suffix}.bam" if RMDUP or KEEP_UNIQ else "{run_dir}/BAM/{dupli}_duplicates/mapq{mapq}/{sample}.mapped.s.bam"
    output:
        "{run_dir}/BAM/{dupli}_duplicates/mapq{mapq}/{sample}.mapped.s{suffix}.regularChr.bam" if RMDUP or KEEP_UNIQ else "{run_dir}/BAM/{dupli}_duplicates/mapq{mapq}/{sample}.mapped.s.regularChr.bam"
    threads: config['nbCores_samtools']
    log:
        "{run_dir}/LOGS/samtools/{dupli}_duplicates/mapq{mapq}/{sample}_{suffix}_keep_regular_chr.txt" if RMDUP or KEEP_UNIQ else "{run_dir}/LOGS/samtools/{dupli}_duplicates/mapq{mapq}/{sample}_keep_regular_chr.txt"
    params:
        quality = config["MAPQ"],
        run_dir = RESULTSDIR,
        suffix = lambda wildcards: create_file_path_suffix(mapq=wildcards.mapq, rmdup=RMDUP, uniq=KEEP_UNIQ, step="keep_regular_chr"),
        suffix2 = lambda wildcards: create_file_path_suffix(mapq=wildcards.mapq, rmdup=RMDUP, uniq=KEEP_UNIQ, regular=True, step="keep_regular_chr"),
        tem_echo_out = tmp_get_echo_out(RMDUP, KEEP_UNIQ, True),
        genome = config['genome'],
        also_with_duplicates = ALSO_WITH_DUPLICATES,
        duplicates = DUPLI
    shell:
        """
        samtools index -@ {threads} {input}
        if [ {params.genome} = hg38 ]
        then
            samtools view --threads {threads} -hb {input} chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY > {output} 2> {log}
        elif [ {params.genome} = mm10 ]
        then
            samtools view --threads {threads} -hb {input} chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX chrY > {output} 2> {log}
        fi
        read_count=$(samtools view -c {output} --threads {threads})
        echo "{params.tem_echo_out} reads with Q>={wildcards.mapq} (_mapped_s{params.suffix2}.bam): $read_count" >> {params.run_dir}/stats/{wildcards.sample}_{params.duplicates}_duplicates_mapq{wildcards.mapq}_stats.txt
        if [ {params.also_with_duplicates} = True ]
        then
            samtools index -@ {threads} {params.run_dir}/BAM/with_duplicates/mapq{wildcards.mapq}/{wildcards.sample}.mapped.s.mapq{wildcards.mapq}.filtered.bam
            if [ {params.genome} = hg38 ]
            then 
                samtools view --threads {threads} -hb {params.run_dir}/BAM/with_duplicates/mapq{wildcards.mapq}/{wildcards.sample}.mapped.s.mapq{wildcards.mapq}.filtered.bam chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY > {params.run_dir}/BAM/with_duplicates/mapq{wildcards.mapq}/{wildcards.sample}.mapped.s.mapq{wildcards.mapq}.filtered.regularChr.bam
            elif [ {params.genome} = mm10 ]
            then
                samtools view --threads {threads} -hb {params.run_dir}/BAM/with_duplicates/mapq{wildcards.mapq}/{wildcards.sample}.mapped.s.mapq{wildcards.mapq}.filtered.bam chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX chrY > {params.run_dir}/BAM/with_duplicates/mapq{wildcards.mapq}/{wildcards.sample}.mapped.s.mapq{wildcards.mapq}.filtered.regularChr.bam
            fi
            read_count=$(samtools view -c {params.run_dir}/BAM/with_duplicates/mapq{wildcards.mapq}/{wildcards.sample}.mapped.s.mapq{wildcards.mapq}.filtered.regularChr.bam --threads {threads})
            if [ {wildcards.mapq} = 0 ]
            then
                multi_or_unique="MULTIPLE"
            else
                multi_or_unique="UNIQUE"
            fi
            echo "$multi_or_unique Reads WITH duplicates with Q>={wildcards.mapq} (_mapped_s.mapq{wildcards.mapq}.filtered.regularChr.bam): $read_count" >> {params.run_dir}/stats/{wildcards.sample}_with_duplicates_mapq{wildcards.mapq}_stats.txt
        fi
        """    

rule final_sort_count:
    input:
        "{run_dir}/BAM/{dupli}_duplicates/mapq{mapq}/{sample}.mapped.s{suffix}.bam" if SUFFIX != "" else "{run_dir}/BAM/{dupli}_duplicates/mapq{mapq}/{sample}.mapped.s.bam"
    output:
        "{run_dir}/BAM/{dupli}_duplicates/mapq{mapq}/{sample}.mapped.s{suffix}.s.bam" if SUFFIX != "" else "{run_dir}/BAM/{dupli}_duplicates/mapq{mapq}/{sample}.mapped.s.s.bam"
    threads: config['nbCores_samtools']
    params:
        quality = MAPQ,
        run_dir = RESULTSDIR,
        suffix = lambda wildcards: create_file_path_suffix(mapq=wildcards.mapq, rmdup=RMDUP, uniq=KEEP_UNIQ, regular=KEEP_REGULAR, step="final_sort_count"),
        tem_echo_out = tmp_get_echo_out(RMDUP, KEEP_UNIQ, KEEP_REGULAR),
        duplicates = DUPLI,
        also_with_duplicates = ALSO_WITH_DUPLICATES
    log:
        "{run_dir}/LOGS/samtools/{dupli}_duplicates/mapq{mapq}/{sample}_{suffix}_final_sort_count.txt" if SUFFIX != "" else "{run_dir}/LOGS/samtools/{dupli}_duplicates/mapq{mapq}/{sample}_final_sort_count.txt"
    shell:
        """
        samtools sort --threads {threads} {input} -o {output} 2> {log}
        samtools index -@ {threads} {output} 2> {log}
        read_count=$(samtools view -c {output} --threads {threads})
        echo "FINALE COUNT with {params.tem_echo_out} sorted reads and Q>={wildcards.mapq} (_Final__Reads): $read_count" >> {params.run_dir}/stats/{wildcards.sample}_{params.duplicates}_duplicates_mapq{wildcards.mapq}_stats.txt
        if [ {params.also_with_duplicates} = True ]
        then
            samtools sort --threads {threads} {params.run_dir}/BAM/with_duplicates/mapq{wildcards.mapq}/{wildcards.sample}.mapped.s.mapq{wildcards.mapq}.filtered.regularChr.bam -o {params.run_dir}/BAM/with_duplicates/mapq{wildcards.mapq}/{wildcards.sample}.mapped.s.mapq{wildcards.mapq}.filtered.regularChr.s.bam
            samtools index -@ {threads} {params.run_dir}/BAM/with_duplicates/mapq{wildcards.mapq}/{wildcards.sample}.mapped.s.mapq{wildcards.mapq}.filtered.regularChr.s.bam
            read_count=$(samtools view -c {params.run_dir}/BAM/with_duplicates/mapq{wildcards.mapq}/{wildcards.sample}.mapped.s.mapq{wildcards.mapq}.filtered.regularChr.s.bam --threads {threads})
            if [ {wildcards.mapq} = 0 ]
            then
                multi_or_unique="MULTIPLE"
            else
                multi_or_unique="UNIQUE"
            fi
            echo "FINALE COUNT with $multi_or_unique sorted reads and Q>={wildcards.mapq} (_Final__Reads): $read_count" >> {params.run_dir}/stats/{wildcards.sample}_with_duplicates_mapq{wildcards.mapq}_stats.txt
        fi
        """
        
rule gather_stats:
    input:
        expand("{run_dir}/BAM/{dupli}_duplicates/mapq{mapq}/{sample}.mapped.s{rmdup}.mapq{mapq}.filtered{regular}.s.bam", sample=SAMPLES, run_dir=RESULTSDIR, mapq=MAPQ, rmdup=".rmdup" if RMDUP else "", regular=".regularChr" if KEEP_REGULAR else "", dupli=DUPLI)
    output:
        expand("{run_dir}/stats/Summary/mapq{q}_{dupli}_duplicates.xlsx", q=MAPQ, run_dir=RESULTSDIR, dupli=DUPLI)
    params:
        run_dir = RESULTSDIR,
        design_matrix = DESIGN_MATRIX,
        rmdup = RMDUP,
        del_tmp_files = config['delete_tmp_files'],
        mapq = MAPQ,
        samples = SAMPLES
    singularity: "docker://hasba/small_stats:latest"
    script:
        "Scripts/gather_stats.py"
                
rule write_parameters:
    output:
        "{run_dir}/{RUN_NAME}_parameters_summary_part_2.txt"
    params:
        RUN_NAME = RUN_NAME,
        run_dir = RESULTSDIR,
        genome = GENOME
    shell:
        """
        echo "##### Pipeline Part 2 #####" > {output}
        echo "" >> {output}
        echo "- Run name: {RUN_NAME}" >> {output}
        echo "" >> {output}
        if [ {ATAC_SEQ_MODE} = True ]
        then
            echo "- Mode: ATAC-seq" >> {output}
        else
            echo "- Mode: CUT&Tag" >> {output}
        fi
        echo "" >> {output}
        echo "- Path to Results Directory": >> {output}
        echo "  -> Results are found in: {RESULTSDIR}" >> {output}
        echo "     --> Part 2 of the pipeline generates the following folders:" >> {output}
        echo "         - BAM: contains the BAM files from Bowtie2 and downstream processing" >> {output}
        echo "         - BAM/with[out]_duplicates/: contains the BAM files with or without duplicates" >> {output}
        echo "         - BAM/with[out]_duplicates/mapq[MAPQ]/: contains the BAM files filtered by MAPQ" >> {output}
        echo "         - stats: contains the stats files for each sample" >> {output}
        echo "         - stats/Summary: contains the summary stats for all samples by condition" >> {output}
        echo "" >> {output}
        echo "- Detected Samples in the following data directory: {DATADIR}" >> {output}
        for sample in {SAMPLES}; do
            echo "  ->" $sample >> {output}
        done 
        echo "" >> {output} 
        echo "- Steps:" >> {output}
        echo "  -> Trimming: {TRIM}" >> {output}
        echo "  -> Indexing genome: {IDX_GNM}" >> {output}
        if [ {IDX_GNM} = True ]
        then
            echo "  -> Indexing genome: {params.genome}" >> {output}
        fi
        echo "  -> Mapping: bowtie2, run with parameters:" >> {output}
        if [ {ATAC_SEQ_MODE} = True ]
        then
            echo "      --> End-to-end, very-sensitive, no-mixed, no-discordant, phred33, X 2000" >> {output}
        else
            echo "      --> End-to-end, very-sensitive, no-mixed, no-discordant, phred33, I 10, X 700" >> {output}
        fi
        echo "      --> Using genome : {params.genome}" >> {output}
        echo "  -> Remove duplicates: {RMDUP}" >> {output}
        if [ {RMDUP} = True ]
        then
            echo "      --> Picard MarkDuplicates, run with parameters: REMOVE_DUPLICATES true, VALIDATION_STRINGENCY LENIENT" >> {output}
        fi        
        if [ {ALSO_WITH_DUPLICATES} = True ]
        then
            echo "      --> Also ouput BAM files with duplicates" >> {output}
        fi
        if [ {KEEP_UNIQ} = True ]
        then
            echo "  -> Filter reads by MAPQ: {KEEP_UNIQ}" >> {output}
            echo "      --> Filtered using MAPQ values = {MAPQ}" >> {output}
        else
            echo "  -> Kept all reads (did not filter using MAPQ)" >> {output}
        fi
        echo "  -> Keep only regular reads: {KEEP_REGULAR}" >> {output}
        if [ {KEEP_REGULAR} = True ]
        then
            if [ {params.genome} = hg38 ]
            then
                echo "      --> Using genome: hg38" >> {output}
                echo "          ---> Kept only reads from chromosomes: chr1, chr2, chr3, chr4, chr5, chr6, chr7, chr8, chr9, chr10, chr11, chr12, chr13, chr14, chr15, chr16, chr17, chr18, chr19, chr20, chr21, chr22, chrX, chrY" >> {output}
            elif [ {params.genome} = mm10 ]
            then
                echo "      --> Using genome: mm10" >> {output}
                echo "         ---> Kept only reads from chromosomes: chr1, chr2, chr3, chr4, chr5, chr6, chr7, chr8, chr9, chr10, chr11, chr12, chr13, chr14, chr15, chr16, chr17, chr18, chr19, chrX, chrY" >> {output}
            fi
        fi
        """
##############################################################################################################