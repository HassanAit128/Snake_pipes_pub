singularity: "docker://hasba/atac_cut:latest"

import glob
import pandas as pd
import os
import re
import csv
import time
from Scripts.shared_functions import check_if_design_matrix_is_valid, get_groups

start = time.time()

############################################# CONFIG ##############################################
RUN_NAME = config['run_name'] # Name of the run/experiment
DOWNSAMPLE = config['downsample'] # Whether to downsample or not
DOWNSAMPLE_GRP = config['downsample_specific_group'] # Group to downsample
DOWNSAMPLE_COL = config['downsample_column'] # Column to downsample
DOWNSAMPLE_TYPE = config['downsample_type'] # Type of downsampling
BIGWIG = config['convert_to_BigWig'] # Whether to convert to BigWig or not
PEAK_CALLING = config['peak_calling'] # Whether to perform peak calling or not
COMPUTE_MATRIX = config['compute_matrix'] # Whether to compute the matrix or not
REGIONS = config['path_to_regions'] # Path to the regions file
PLOT_HEATMAP = config['plot_heatmap'] # Whether to plot the heatmap or not
PLOT_PROFILE = config['plot_profile'] # Whether to plot the profile or not
CUSTOM_REGIONS = config['custom_regions'] # Path to the custom regions file for the matrix
DIRECTORY = config['working_directory'] + "/" + RUN_NAME # Directory containing the data
DESIGN_MATRIX_FILE = config['design_matrix'] # Path to the design matix file
DESIGN_MATRIX = pd.read_csv(DESIGN_MATRIX_FILE)  # Read the design matrix file into a DataFrame
FILE_LIST = glob.glob(DIRECTORY + "/BAM/*/*/*.s.bam") # List of BAM files
DIFFBIND = config['diffbind'] # Whether to perform diffbind or not
USER_PEAKS = config['user_peaks'] # Path to the user peaks file
USER_PEAKS_CALLER = config['user_peak_caller'] # Peak caller used for the user peaks file
MAPQ = config['MAPQ'] # Mapping quality
if MAPQ == "":
    MAPQ = [0]
if type(MAPQ) == int:
    MAPQ = [MAPQ]
RMDUP = config['remove_duplicates'] # Whether to remove duplicates or not
DUPLI = "without" if RMDUP else "with" # Duplicates or not
ALSO_WITH_DUPLICATES = config['output_finale_file_with_duplicates'] # Also run analysis with duplicates True/False
KEEP_UNIQ = True
KEEP_REGULAR = config['keep_only_regular_reads'] # Keep only regular reads True/False
ATAC_SEQ_MODE = config['atac_seq'] # ATAC-seq mode True/False
MACS2_PARAMS = config['macs2_parameters'] # MACS2 parameters
TOBIAS = config['run_tobias'] # Whether to run TOBIAS or not
#####################################################################################################

############################################# ERROR CHECK ###########################################
samples_check = [os.path.basename(f).replace(".bam", "").split(".")[0] for f in FILE_LIST]
check_if_design_matrix_is_valid(DESIGN_MATRIX_FILE, samples_check)

if RUN_NAME == "":
    raise Exception("Please provide a name for the project/run in the config file.")
if not FILE_LIST:
    raise Exception("No files found. Please run the previous part of the pipeline to generate the BAM files.")
if DOWNSAMPLE == "Group" and (DOWNSAMPLE_GRP == None or DOWNSAMPLE_GRP == "" or DOWNSAMPLE_COL == None or DOWNSAMPLE_COL == ""):
    raise Exception("Please specify the group/column to downsample in the config file.")
if DOWNSAMPLE_COL != "" and DOWNSAMPLE_COL not in DESIGN_MATRIX.columns:
    raise Exception("The column specified for downsampling is not present in the design matrix file. Please check the column name in the config file.")
if DIFFBIND :
    if not PEAK_CALLING and (USER_PEAKS == None or USER_PEAKS == ""):
        raise Exception("Cannot perform DiffBind without peak calling or peaks file. Please set PEAK_CALLING to True, OR, add a path to custom peaks file in the config file.")         
    if USER_PEAKS != "" and USER_PEAKS_CALLER == "":
        raise Exception("Please specify the peak caller used for the custom peaks file in the config file.") 
if USER_PEAKS != "" and not os.path.exists(USER_PEAKS):
    raise Exception("The path to the custom peaks file is incorrect. Please check the path in the config file.")
if not os.path.exists(DESIGN_MATRIX_FILE):
    raise Exception("The path to the design matrix file is incorrect. Please check the path in the config file.")
if MACS2_PARAMS == "" and PEAK_CALLING:
    print("No MACS2 parameters specified. Default parameters will be used.") 


if (PLOT_HEATMAP or PLOT_PROFILE) and not COMPUTE_MATRIX and not os.path.exists(REGIONS): 
    raise Exception("Cannot plot the heatmap or profile without computing the matrix. Please set COMPUTE_MATRIX to True, OR, add a path to the regions file in the config file.")

onsuccess:
    time_to_finish = round((time.time() - start)/60,1)
    longest_line_length = max(len(f"OUTPUT folder is found at: {DIRECTORY}"), len("Running time in minutes: %s " % time_to_finish))
    total_length = longest_line_length + 8
    hash_count = (total_length - len("Workflow finished, no error") - 4) // 2
    print("\n" + "#" * hash_count + "# " + "WORKFLOW FINISHED, NO ERROR" + " #" + "#" * hash_count)
    print(f"OUTPUT folder is found at: {DIRECTORY}")
    print(f"Running time in minutes: {time_to_finish}")
    print(f"Mode: {'ATAC-seq' if ATAC_SEQ_MODE else 'CUT&Tag'}")
    print(f"Downsample: {DOWNSAMPLE}")
    if DOWNSAMPLE:
        print(f"Downsample type: {DOWNSAMPLE_TYPE}")
        if DOWNSAMPLE_TYPE != "ALL":
            print(f"    Downsample column: {DOWNSAMPLE_COL}")
            print(f"    Downsample by group: {DOWNSAMPLE_GRP if DOWNSAMPLE_GRP != '' else 'None'}")
            print(f"    Number of downsampled files: {len(files_to_downsample)}")
            print(f"    Number of non-downsampled files: {len(files_not_to_downsample)}")
    print(f"Convert to BigWig: {BIGWIG}")
    print(f"Peak calling: {PEAK_CALLING}")
    if PEAK_CALLING:
        print(f"    Peak caller: MACS2")
        print(f"    MACS2 parameters: {MACS2_PARAMS if MACS2_PARAMS != '' else 'Default parameters used (-f BAMPE --max-gap 2000 --min-length 200 -B --broad --broad-cutoff 0.1)'}")
    print(f"DiffBind: {DIFFBIND}")
    if DIFFBIND:
        print(f"    Peak caller: {'MACS2' if USER_PEAKS == '' else USER_PEAKS_CALLER}")
        print(f"    Custom peaks file: {USER_PEAKS if USER_PEAKS != '' else 'None'}")
    if TOBIAS:
        print(f"TOBIAS: {TOBIAS}")
    print("#" * total_length)

onerror:
    print("\n\n###################### An error occurred ######################\n")
    print("Running time in minutes: %s\n" % round((time.time() - start)/60,1))
    print("\n###############################################################\n\n")
###############################################################################################

######################################### HELPER FUNCTIONS ####################################
def get_files_to_downsample():
    """
    Get the files to downsample and the files not to downsample
    
    Parameters:
    -----------
        None

    Returns:
    --------
        files_to_downsample: list
            List of files to downsample
        files_not_to_downsample: list
            List of files not to downsample
    """
    if DOWNSAMPLE:
        if DOWNSAMPLE_TYPE == "ALL":
            files_to_downsample = FILE_LIST
            files_not_to_downsample = []
        else:
            if DOWNSAMPLE_GRP != "":
                samples_to_downsample = DESIGN_MATRIX[DESIGN_MATRIX[DOWNSAMPLE_COL] == DOWNSAMPLE_GRP]["SampleID"].tolist()
                files_to_downsample = [f for f in FILE_LIST if os.path.basename(f).replace(".bam", "").split(".")[0] in samples_to_downsample]
                files_not_to_downsample = [f for f in FILE_LIST if f not in files_to_downsample]
            else:
                files_to_downsample = FILE_LIST
                files_not_to_downsample = []
    else :
        files_to_downsample = []
        files_not_to_downsample = FILE_LIST
    return files_to_downsample, files_not_to_downsample

def generate_rule_all_input(files):
    """
    Generate the input for the rule all

    Parameters:
    -----------
        files: list
            List of files
    
    Returns:
    --------
        rule_all_input: list
            List of files to be used as input for the rule all
        downsampled_samples: list
            List of downsampled samples
        non_downsampled_samples: list
            List of non-downsampled samples
    """
    rule_all_input = []
    downsampled_samples = []
    non_downsampled_samples = []
    for f in FILE_LIST:
        base_path = DIRECTORY + "/BAM"
        sample_name = os.path.relpath(f, base_path)
        sample_name = sample_name.replace(".bam", "")
        if f in files:
            rule_all_input.append(f"{DIRECTORY}/DS/{sample_name}.downsampled.bam")
            downsampled_samples.append(sample_name)
        else:
            non_downsampled_samples.append(sample_name)
        if BIGWIG:
            rule_all_input.append(f"{DIRECTORY}/bigwig/{sample_name}.bw")
        if PEAK_CALLING:
            if MACS2_PARAMS != "":
                if "broad" in MACS2_PARAMS:
                    rule_all_input.append(f"{DIRECTORY}/peak_calling/{sample_name}_peaks.broadPeak")
                else:
                    rule_all_input.append(f"{DIRECTORY}/peak_calling/{sample_name}_peaks.narrowPeak")
            else:
                rule_all_input.append(f"{DIRECTORY}/peak_calling/{sample_name}_peaks.broadPeak")
        if DIFFBIND:
            diffbinds = expand("{run_dir}/DiffBind/diffbind_mapq{q}_{dupli}_duplicates.csv", run_dir = DIRECTORY, q = MAPQ, dupli = ['without', 'with'] if RMDUP and ALSO_WITH_DUPLICATES else ['without'] if RMDUP else ['with'])
            rule_all_input.extend(diffbinds)
            diffbinds_rdata = expand("{run_dir}/DiffBind/diffbind_mapq{q}_{dupli}_duplicates.Rdata", run_dir = DIRECTORY, q = MAPQ, dupli = ['without', 'with'] if RMDUP and ALSO_WITH_DUPLICATES else ['without'] if RMDUP else ['with'])
            rule_all_input.extend(diffbinds_rdata)
            plots = expand("{run_dir}/DiffBind/diffbind_mapq{q}_{dupli}_duplicates_diffbind_annotated_TE_summary.png", run_dir = DIRECTORY, q = MAPQ, dupli = ['without', 'with'] if RMDUP and ALSO_WITH_DUPLICATES else ['without'] if RMDUP else ['with'])
            rule_all_input.extend(plots)
        if COMPUTE_MATRIX:
            matrices = expand("{run_dir}/Matrix/mapq{q}_{dupli}_duplicates_matrix.gz", run_dir = DIRECTORY, q = MAPQ, dupli = ['without', 'with'] if RMDUP and ALSO_WITH_DUPLICATES else ['without'] if RMDUP else ['with'])
            rule_all_input.extend(matrices)
        if TOBIAS:
            rule_all_input.append(f"{DIRECTORY}/TOBIAS/BINDetect_output/bindetect_results.txt")
    rule_all_input.append(f"{DIRECTORY}/{RUN_NAME}_parameters_summary_part_3.txt")
    return rule_all_input, downsampled_samples, non_downsampled_samples

files_to_downsample, files_not_to_downsample = get_files_to_downsample()
rule_all_input, downsampled_samples, non_downsampled_samples = generate_rule_all_input(files_to_downsample)
################################################################################################

############################################# RULES ############################################
rule all:
    input:
        rule_all_input,
        f"{DIRECTORY}/{RUN_NAME}_report.html"
    
rule report:
    input:
        rule_all_input
    output:
        "{run_dir}/{run}_report.html"
    params:
        run_dir = DIRECTORY,
        run = RUN_NAME,
        summaries = glob.glob(f"{DIRECTORY}/stats/Summary/*"),
        dir = DIRECTORY,
        total_samples = len(downsampled_samples) + len(non_downsampled_samples),
        multiqc_path = f"{DIRECTORY}/multiqc_report.html",
        multiqc_trimmed_path = f"{DIRECTORY}/Trimmed_Fastq/multiqc_report.html",
        Mode = "CUT&Tag" if not ATAC_SEQ_MODE else "ATAC-seq"
    singularity: "docker://hasba/pyrpsta:latest"
    shell:
        """
        python Scripts/report.py --run {params.run} --summaries {params.summaries} --dir {params.run_dir} --total_samples {params.total_samples} --multiqc_path {params.multiqc_path} --multiqc_trimmed_path {params.multiqc_trimmed_path} --Mode "{params.Mode}"
        """

rule find_read_count:
    input:
        lambda wildcards: f"{DIRECTORY}/BAM/{wildcards.sample}.bam"
    output:
        temp("{run_dir}/BAM/{sample}_count.txt")
    threads: config['nbCores_samtools']
    params:
        run_dir = DIRECTORY
    shell:
        """
        mkdir -p {params.run_dir}
        mkdir -p {params.run_dir}/DS
        samtools view --threads {threads} -c {input} > {output}
        """
        
rule create_downsample_cmds:
    input:
        files = expand("{run_dir}/BAM/{sample}.bam", sample=[f for f in downsampled_samples], run_dir = DIRECTORY),
        counts = expand("{run_dir}/BAM/{sample}_count.txt", sample=[f for f in downsampled_samples], run_dir = DIRECTORY)
    output:
        dss = expand("{run_dir}/DS/{sample}.downsampled.bam", sample=[f for f in downsampled_samples], run_dir = DIRECTORY)
    params:
        ds_type = DOWNSAMPLE_TYPE,
        ds_col = DOWNSAMPLE_COL,
        ds_grp = DOWNSAMPLE_GRP,
        des_mat = DESIGN_MATRIX_FILE,
        run_dir = DIRECTORY,
        samples = downsampled_samples,
        mapq = MAPQ,
        memory = config['available_memory'],
        cores = config['available_cores'],
    script:
        "Scripts/downsample.py"

rule convert_to_BigWig:
    input:
        lambda wildcards: f"{DIRECTORY}/DS/{wildcards.sample}.downsampled.bam" if wildcards.sample in downsampled_samples else f"{DIRECTORY}/BAM/{wildcards.sample}.bam"
    output:
        "{run_dir}/bigwig/{sample}.bw"
    threads: config['nbCores_deepTools']
    params: 
        nbCores_samtools = config['nbCores_samtools'],
        nbCores_deepTools = config['nbCores_deepTools'],
        run_dir = DIRECTORY
    log:
        "{run_dir}/LOGS/{sample}_bamCoverage.log"
    shell:
        """
        mkdir -p {params.run_dir}/bigwig
        samtools index -@ {params.nbCores_samtools} {input}
        bamCoverage -b {input} -o {output} -bs 5 --normalizeUsing BPM -p {threads}
        """
        
rule peak_calling:
    input:
        lambda wildcards: f"{DIRECTORY}/DS/{wildcards.sample}.downsampled.bam" if wildcards.sample in downsampled_samples else f"{DIRECTORY}/BAM/{wildcards.sample}.bam"
    output:
        "{run_dir}/peak_calling/{sample}_peaks.broadPeak" if "broad" in MACS2_PARAMS else "{run_dir}/peak_calling/{sample}_peaks.narrowPeak"
    params:
        run_dir = DIRECTORY,
        mode_params = MACS2_PARAMS,
        genome = "mm" if config['genome'] == "mm10" else "hs"
    log:
        "{run_dir}/LOGS/MACS2/{sample}_MACS2.log"
    singularity: "docker://fooliu/macs2"
    shell:
        """
        mkdir -p {params.run_dir}/peak_calling
        macs2 callpeak -t {input} --outdir {DIRECTORY}/peak_calling -n {wildcards.sample} -g {params.genome} {params.mode_params} &> {log}
        """

rule create_diffbind_input:
    input:
        files = expand("{run_dir}/peak_calling/{sample}_peaks.{prefix}Peak", sample=[f for f in downsampled_samples] + [f for f in non_downsampled_samples], run_dir = DIRECTORY, prefix = "broad" if "broad" in MACS2_PARAMS else "narrow") if USER_PEAKS == "" else expand("{run_dir}/DS/{sample}.downsampled.bam", sample=[f for f in downsampled_samples], run_dir = DIRECTORY) + expand("{run_dir}/BAM/{sample}.bam", sample=[f for f in non_downsampled_samples], run_dir = DIRECTORY)
    output:
        expand("{run_dir}/DiffBind/diffbind_mapq{q}_{dupli}_duplicates.csv", run_dir = DIRECTORY, q = MAPQ, dupli = ['without', 'with'] if RMDUP and ALSO_WITH_DUPLICATES else ['without'] if RMDUP else ['with'])
    params:
        run_dir = DIRECTORY,
        design_matrix = DESIGN_MATRIX_FILE,
        samples = FILE_LIST,
        downsampled_samples = downsampled_samples,
        peak_called = PEAK_CALLING,
        peaks_path = USER_PEAKS,
        peaks_path_caller = USER_PEAKS_CALLER,
        mapq = MAPQ
    script:
        "Scripts/create_diffbind_input.py"

rule run_diffbind:
    input:
        matrices = expand("{run_dir}/DiffBind/diffbind_mapq{q}_{dupli}_duplicates.csv", run_dir = DIRECTORY, q = MAPQ, dupli = ['without', 'with'] if RMDUP and ALSO_WITH_DUPLICATES else ['without'] if RMDUP else ['with'])
    output:
        expand("{run_dir}/DiffBind/diffbind_mapq{q}_{dupli}_duplicates.Rdata", run_dir = DIRECTORY, q = MAPQ, dupli = ['without', 'with'] if RMDUP and ALSO_WITH_DUPLICATES else ['without'] if RMDUP else ['with'])
    params:
        run_dir = DIRECTORY,
        pvalue = config['diffbind_pvalue'],
        mode = 'cutntag' if not ATAC_SEQ_MODE else 'atac-seq',
        flip = config['diffbind_flip']
    script:
        "Scripts/Diffbind.R"

rule annotate_peaks:
    input:
        rdata = expand("{run_dir}/DiffBind/diffbind_mapq{q}_{dupli}_duplicates.Rdata", run_dir = DIRECTORY, q = MAPQ, dupli = ['without', 'with'] if RMDUP and ALSO_WITH_DUPLICATES else ['without'] if RMDUP else ['with'])    
    output:
        annotated = expand("{run_dir}/DiffBind/diffbind_mapq{q}_{dupli}_duplicates_diffbind_annotated.bed", run_dir = DIRECTORY, q = MAPQ, dupli = ['without', 'with'] if RMDUP and ALSO_WITH_DUPLICATES else ['without'] if RMDUP else ['with'])
    params:
        run_dir = DIRECTORY,
        overlaps_param = config['overlaps_parameters'],
        peaks = expand("{run_dir}/DiffBind/diffbind_mapq{q}_{dupli}_duplicates_diffbind.bed", run_dir = DIRECTORY, q = MAPQ, dupli = ['without', 'with'] if RMDUP and ALSO_WITH_DUPLICATES else ['without'] if RMDUP else ['with']),
        regions = REGIONS
    shell:
        """
        for peak_file in {params.peaks}; do
            echo "Processing $peak_file"
            intersectBed -a $peak_file -b {params.regions} -wao {params.overlaps_param} > {params.run_dir}/DiffBind/$(basename $peak_file .bed)_annotated.bed
        done
        """

rule plot_annotation:
    input:
        annotated_peaks = expand("{run_dir}/DiffBind/diffbind_mapq{q}_{dupli}_duplicates_diffbind_annotated.bed", run_dir = DIRECTORY, q = MAPQ, dupli = ['without', 'with'] if RMDUP and ALSO_WITH_DUPLICATES else ['without'] if RMDUP else ['with'])
    output:
        expand("{run_dir}/DiffBind/diffbind_mapq{q}_{dupli}_duplicates_diffbind_annotated_TE_summary.png", run_dir = DIRECTORY, q = MAPQ, dupli = ['without', 'with'] if RMDUP and ALSO_WITH_DUPLICATES else ['without'] if RMDUP else ['with'])
    params:
        run_dir = DIRECTORY
    script:
        "Scripts/a_plots.R"      

        
rule compute_matrix:
    input:
        bw = expand("{run_dir}/bigwig/{sample}.bw", sample=[f for f in downsampled_samples] + [f for f in non_downsampled_samples], run_dir = DIRECTORY),
        annotated_peaks = expand("{run_dir}/DiffBind/diffbind_mapq{q}_{dupli}_duplicates_diffbind_annotated.bed", run_dir = DIRECTORY, q = MAPQ, dupli = ['without', 'with'] if RMDUP and ALSO_WITH_DUPLICATES else ['without'] if RMDUP else ['with'])
    output:
        matrix = expand("{run_dir}/Matrix/mapq{q}_{dupli}_duplicates_matrix.gz", run_dir = DIRECTORY, q = MAPQ, dupli = ['without', 'with'] if RMDUP and ALSO_WITH_DUPLICATES else ['without'] if RMDUP else ['with'])
    threads: config['nbCores_deepTools'] if config['nbCores_deepTools'] > 0 else 1
    params:
        run_dir = DIRECTORY,
        path_to_regions = REGIONS,
        addit_params = config['compute_matrix_parameters'] if config['compute_matrix_parameters'] != "" else "None",
        type_a = config['compute_matrix_type'],
        specific = config['specific_regions'] if config['specific_regions'] != "" else "None",
        design_matrix = DESIGN_MATRIX_FILE,
        plot_heatmap = "True" if PLOT_HEATMAP else "False",
        ph_addit_params = config['plot_heatmap_parameters'] if config['plot_heatmap_parameters'] != "" else "None",
        plot_profile = "True" if PLOT_PROFILE else "False",
        pp_addit_params = config['plot_profile_parameters'] if config['plot_profile_parameters'] != "" else "None",
        custom_regions = CUSTOM_REGIONS if CUSTOM_REGIONS != "" else "None",
        upstream = config['upstream'],
        downstream = config['downstream']
    shell:
        """
        mkdir -p {params.run_dir}/Matrix
        python Scripts/compute_matrix.py -b {input.bw} -a {input.annotated_peaks} -sr {params.specific} -d {params.design_matrix} -r {params.run_dir} -t {params.type_a} -ao "{params.addit_params}" -pp {params.plot_profile} -ppa "{params.pp_addit_params}" -ph {params.plot_heatmap} -pha "{params.ph_addit_params}" -cr {params.custom_regions} -p {threads} -up {params.upstream} -down {params.downstream}
        """

unique_conditions, condition1_samples, condition2_samples = get_groups(DESIGN_MATRIX_FILE)
bam1, bam2 = "", ""
if len(unique_conditions) == 2:
    bam1 = f"{DIRECTORY}/TOBIAS/{unique_conditions[0]}.bam"
    bam2 = f"{DIRECTORY}/TOBIAS/{unique_conditions[1]}.bam"

rule tobias_prep:
    input:
        matrix = expand("{run_dir}/DS/{sample}.downsampled.bam", sample=[f for f in downsampled_samples], run_dir = DIRECTORY),
    output:
        expand("{run_dir}/TOBIAS/{sample}.bam", sample=[f for f in [bam1.split('/')[-1].split(".")[0], bam2.split('/')[-1].split(".")[0]]], run_dir = DIRECTORY)
    params:
        run_dir = DIRECTORY,
        des_mat = DESIGN_MATRIX_FILE,
        genome = config['path_to_genome_fasta'],
        macs2_params = MACS2_PARAMS,
        uropa_config = config['uropa_config'],
        motifs = config['motifs'],
        blacklist = config['blacklist'],
        mapq = MAPQ,
        genomesize = "mm" if config['genome'] == "mm10" else "hs"
    threads: config['nbCores_TOBIAS']
    singularity: "docker://hasba/toburo"
    shell:
        """
        mkdir -p {params.run_dir}/TOBIAS
        python Scripts/tobias.py -q {params.mapq} -d {params.run_dir} -dm {params.des_mat} -b {input} -g {params.genome} -p "{params.macs2_params}" -c {params.uropa_config} -m {params.motifs} -bl {params.blacklist} --cores {threads} -gs {params.genomesize}
        """    
    
rule tobias_atac_correct:
    input:
        bamfiles = [bam1, bam2]
    output:
        expand("{run_dir}/TOBIAS/ATACorrect/{sample}_corrected.bw", sample=[b for b in [bam1.split('/')[-1].split('.')[0], bam2.split('/')[-1].split('.')[0]]], run_dir = DIRECTORY)
    params:
        run_dir = DIRECTORY,
        genome = config['path_to_genome_fasta'],
        peaks = REGIONS,
        blacklist = config['blacklist'],
    threads: config['nbCores_TOBIAS']
    singularity: "docker://hasba/toburo"
    shell:
        """
        mkdir -p {params.run_dir}/TOBIAS/ATACorrect
        for bam in {input.bamfiles}; do
            echo "---> Processing $bam"
            TOBIAS ATACorrect --bam $bam --genome {params.genome} --peaks {params.peaks} --blacklist {params.blacklist} --outdir {params.run_dir}/TOBIAS/ATACorrect --cores {threads}
        done
        """

rule tobias_footprint_scores:
    input:
        bws = expand("{run_dir}/TOBIAS/ATACorrect/{sample}_corrected.bw", sample=[b for b in [bam1.split('/')[-1].split('.')[0], bam2.split('/')[-1].split('.')[0]]], run_dir = DIRECTORY)
    output:
        expand("{run_dir}/TOBIAS/FootprintScores/{sample}_footprints.bw", sample=[b for b in [bam1.split('/')[-1].split('.')[0], bam2.split('/')[-1].split('.')[0]]], run_dir = DIRECTORY)
    params:
        run_dir = DIRECTORY,
        cores = config['nbCores_TOBIAS'],
        peaks = REGIONS
    singularity: "docker://hasba/toburo"
    shell:
        """
        mkdir -p {params.run_dir}/TOBIAS/FootprintScores
        for sample in {params.samples}; do
            bw={params.run_dir}/TOBIAS/ATACorrect/${{sample}}_corrected.bw
            echo "---> Processing $bw"
            TOBIAS FootprintScores --signal $bw --regions {params.peaks} --output {params.run_dir}/TOBIAS/FootprintScores/$sample_footprints.bw --cores {params.cores}
        done
        """

rule tobias_bindetect:
    input:
        expand("{run_dir}/TOBIAS/FootprintScores/{sample}_footprints.bw", sample=[b for b in [bam1.split('/')[-1].split('.')[0], bam2.split('/')[-1].split('.')[0]]], run_dir = DIRECTORY)
    output:
        "{run_dir}/TOBIAS/BINDetect_output/bindetect_results.txt"
    params:
        run_dir = DIRECTORY,
        motifs = config['motifs'],
        genome = config['path_to_genome_fasta'],
        peaks = REGIONS,
    singularity: "docker://hasba/toburo"
    threads: config['nbCores_TOBIAS']
    shell:
        """
        mkdir -p {params.run_dir}/TOBIAS/BINDetect_output
        TOBIAS BINDetect --motifs {params.motifs} --signals {input} --genome {params.genome} --peaks {params.peaks} --peak_header {params.run_dir}/TOBIAS/merged_peaks_annotated_header.txt --outdir {params.run_dir}/TOBIAS/BINDetect_output --cores {threads}
        """

rule write_parameters:
    output:
        "{run_dir}/{run_name}_parameters_summary_part_3.txt"
    params:
        run_dir = DIRECTORY,
        run_name = RUN_NAME
    shell:
        """
        echo "##### Pipeline Part 3 #####" > {output}
        echo "" >> {output}
        echo "- Run name: {RUN_NAME}" >> {output}
        echo "" >> {output}
        echo "- Path to Results Directory": >> {output}
        echo "  -> Results are found in: {DIRECTORY}" >> {output}
        echo "     --> Part 3 of the pipeline generates the following folders:" >> {output}
        if [ "{DOWNSAMPLE}" = "True" ]
        then
            echo "         - DS: contains the downsampled BAM files" >> {output}
            echo "         - Downsampling statistics can be found in the stats/Summary folder" >> {output}
        fi
        if [ "{BIGWIG}" = "True" ]
        then
            echo "         - bigwig: contains the BigWig files" >> {output}
        fi
        if [ "{PEAK_CALLING}" = "True" ]
        then
            echo "         - peak_calling: contains the peak calling results" >> {output}
        fi
        if [ "{DIFFBIND}" = "True" ]
        then
            echo "         - DiffBind: contains the DiffBind matrix for each condition" >> {output}
        fi
        echo "" >> {output}
        echo "- Detected BAM files": >> {output}  
        for file in {FILE_LIST}; do
            echo "  ->" $(basename $file) >> {output}
        done 
        echo "" >> {output} 
        echo "- Steps:" >> {output}
        echo "  -> Downsample: {DOWNSAMPLE}" >> {output}
        if [ {DOWNSAMPLE} = True ]
        then
            echo "      --> Downsample type: {DOWNSAMPLE_TYPE}" >> {output}
        fi
        if [ {DOWNSAMPLE_TYPE} = "Group" ]
        then
            echo "      --> Downsample column: {DOWNSAMPLE_COL}" >> {output}
            echo "      --> Downsample by group: {DOWNSAMPLE_GRP}" >> {output}
            echo "      --> Files to downsample:" >> {output}
            if [ -z "{files_to_downsample}" ]
            then
                echo "          ---> None" >> {output}
            else
                for file in {files_to_downsample}; do
                    echo "          --->" $(basename $file) >> {output}
                done
            fi
            echo "      --> Files not to downsample:" >> {output}
            if [ -z "{files_not_to_downsample}" ]
            then
                echo "          ---> None" >> {output}
            else
                for file in {files_not_to_downsample}; do
                    echo "          --->" $(basename $file) >> {output}
                done
            fi
        fi
        echo "  -> Convert to BigWig: {BIGWIG}" >> {output}
        if [ "{BIGWIG}" = "True" ]
        then
            echo "      --> Method: bamCoverage, run with parameters:" >> {output}
            echo "          ---> -bs 5 --normalizeUsing BPM" >> {output}
        fi
        echo "  -> Peak calling: {PEAK_CALLING}" >> {output}
        if [ "{PEAK_CALLING}" = "True" ]
        then
            echo "      --> Method: MACS2, run with parameters:" >> {output}
            if [ "{MACS2_PARAMS}" = "" ]
            then
                echo "          ---> Default parameters used : -f BAMPE --max-gap 2000 --min-length 200 -B --broad --broad-cutoff 0.1" >> {output}
            else
                echo "          ---> {MACS2_PARAMS}" >> {output}
            fi
        fi
        echo "  -> Run DiffBind: {DIFFBIND}" >> {output}
        if [ "{DIFFBIND}" = "True" ]
        then
            if [ "{USER_PEAKS}" = "" ]
            then
                echo "      --> Peak caller: MACS2" >> {output}
            else
                echo "      --> Used custom peaks file: {USER_PEAKS}" >> {output}
                echo "      --> Peak caller: {USER_PEAKS_CALLER}" >> {output}
            fi
        fi
        """