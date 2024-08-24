include: "common.smk"
path_to_genome_folder = os.path.dirname(GENOME_PATH)

rule bismark_index:
    input:
        genome = path_to_genome_folder
    output:
        idx = idx_output
    params:
        run_dir = PROJECT_DIR
    log:
        log = "{run_dir}/LOGS/IDX_LOGS/bismark_idx.log"
    threads: config['bismark']['bismark_cores']
    shell:
        """
        bismark_genome_preparation --parallel {threads} {input.genome} > {log} 2>&1
        """

rule align:
    input:
        R1 = lambda wildcards: f"{PROJECT_DIR}/TRIMMED_READS/{wildcards.sample}_R1{get_suffixes(wildcards.sample, RAW_DATA_DIR, trim = True)[0]}val_1.fq.gz" if TRIM else f"{RAW_DATA_DIR}/{wildcards.sample}_R1{get_suffixes(wildcards.sample, RAW_DATA_DIR)[0]}" + get_suffixes(wildcards.sample, RAW_DATA_DIR)[1],
        R2 = lambda wildcards: f"{PROJECT_DIR}/TRIMMED_READS/{wildcards.sample}_R2{get_suffixes(wildcards.sample, RAW_DATA_DIR, trim = True)[0]}val_2.fq.gz" if TRIM else f"{RAW_DATA_DIR}/{wildcards.sample}_R2{get_suffixes(wildcards.sample, RAW_DATA_DIR)[0]}" + get_suffixes(wildcards.sample, RAW_DATA_DIR)[1]
    output:
        bam = "{run_dir}/ALIGN_BAM/{sample}_aligned_pe.bam"
    params:
        genome = path_to_genome_folder,
        run_dir = PROJECT_DIR,
        additional_params = ALIGN_ARGS,
    threads: ALIGN_CORES
    log:
        log = "{run_dir}/LOGS/ALIGN_LOGS/{sample}_bismark.log"
    shell:
        """
        bismark -q --gzip --parallel {threads} --genome {params.genome} {params.additional_params} -1 {input.R1} -2 {input.R2} -o {params.run_dir}/ALIGN_BAM > {log} 2>&1
        mv {params.run_dir}/ALIGN_BAM/{wildcards.sample}*_report.txt {params.run_dir}/REPORTS/ALIGN_REPORTS
        mv {params.run_dir}/ALIGN_BAM/{wildcards.sample}*_bismark_bt2_pe.bam {output.bam}
        """

rule dedup:
    input:
        bam = lambda wildcards: f"{PROJECT_DIR}/ALIGN_BAM/{wildcards.sample}_aligned_pe.bam" 
    output:
        dedup_bam = "{run_dir}/BAMs/{sample}.deduplicated.bam"
    params:
        run_dir = PROJECT_DIR,
        additional_params = DEDUP_ARGS,
        trim = TRIM
    log:
        log = "{run_dir}/LOGS/DEDUP_LOGS/{sample}_dedup.log"
    shell:
        """
        deduplicate_bismark --paired {params.additional_params} --output_dir {params.run_dir}/BAMs {input.bam} > {log} 2>&1
        mv {params.run_dir}/BAMs/{wildcards.sample}_*deduplication_report.txt {params.run_dir}/REPORTS/DEDUP_REPORTS
        if [ "{params.trim}" = "True" ]; then
            cp {params.run_dir}/BAMs/{wildcards.sample}_aligned_pe.deduplicated.bam {output.dedup_bam}
            rm {params.run_dir}/BAMs/{wildcards.sample}_aligned_pe.deduplicated.bam
        else
            cp {params.run_dir}/BAMs/{wildcards.sample}_aligned_pe.deduplicated.bam {output.dedup_bam}
            rm {params.run_dir}/BAMs/{wildcards.sample}_aligned_pe.deduplicated.bam
        fi         
        """ 

rule mbias_only:
    input:
        bam = lambda wildcards: f"{PROJECT_DIR}/BAMs/{wildcards.sample}.deduplicated.regular.sorted.bam" if not DOWNSAMPLE else f"{PROJECT_DIR}/BAMs/{wildcards.sample}.deduplicated.regular.ds.sorted.bam"
    output:
        mbias = "{run_dir}/METHYLATION/METHYLATION_BIAS/{sample}.deduplicated.regular.sorted.M-bias.txt" if not DOWNSAMPLE else "{run_dir}/METHYLATION/METHYLATION_BIAS/{sample}.deduplicated.regular.ds.sorted.M-bias.txt"
    params:
        run_dir = PROJECT_DIR,
        samtools_cores = SAMTOOLS_CORES,
        additional_params = CALL_METH_ARGS,
        suffix = lambda wildcards: f"{wildcards.sample}.deduplicated.regular.sorted" if not DOWNSAMPLE else f"{wildcards.sample}.deduplicated.regular.ds.sorted"
    threads: CALL_METH_CORES
    log:
        log = "{run_dir}/LOGS/METHYLATION_LOGS/METHYLATION_BIAS_LOGS/{sample}_mbias.log"
    shell:
        """
        mkdir -p {params.run_dir}/TEMP/MBIAS_ONLY
        samtools sort -@ {params.samtools_cores} -n -o {params.run_dir}/TEMP/MBIAS_ONLY/{params.suffix}.bam {input.bam} 
        bismark_methylation_extractor {params.additional_params} --mbias_only --parallel {threads} --output_dir {params.run_dir}/METHYLATION/METHYLATION_BIAS {params.run_dir}/TEMP/MBIAS_ONLY/{params.suffix}.bam > {log} 2>&1
        mv {params.run_dir}/METHYLATION/METHYLATION_BIAS/{wildcards.sample}*_report.txt {params.run_dir}/REPORTS/METHYLATION_REPORTS/METHYLATION_BIAS_REPORTS
        """

rule draw_mbias_plot:
    input:
        mbias = lambda wildcards: f"{PROJECT_DIR}/METHYLATION/METHYLATION_BIAS/{wildcards.sample}.deduplicated.regular.sorted.M-bias.txt" if not DOWNSAMPLE else f"{PROJECT_DIR}/METHYLATION/METHYLATION_BIAS/{wildcards.sample}.deduplicated.regular.ds.sorted.M-bias.txt"
    output:
        plot = "{run_dir}/METHYLATION/METHYLATION_BIAS/{sample}.deduplicated.regular.sorted.M-bias_R1.png" if not DOWNSAMPLE else "{run_dir}/METHYLATION/METHYLATION_BIAS/{sample}.deduplicated.regular.ds.sorted.M-bias_R1.png",
        plot2 = "{run_dir}/METHYLATION/METHYLATION_BIAS/{sample}.deduplicated.regular.sorted.M-bias_R2.png" if not DOWNSAMPLE else "{run_dir}/METHYLATION/METHYLATION_BIAS/{sample}.deduplicated.regular.ds.sorted.M-bias_R2.png"
    params:
        run_dir = PROJECT_DIR
    log:    
        log = "{run_dir}/LOGS/METHYLATION_LOGS/METHYLATION_BIAS_LOGS/{sample}_mbias_plot.log"
    shell:
        """
        python Scripts/plot_mbias.py -p {input.mbias}
        """

rule methylation_calling:
    input:
        bam = lambda wildcards: f"{PROJECT_DIR}/BAMs/{wildcards.sample}.deduplicated.regular.sorted.bam" if not DOWNSAMPLE else f"{PROJECT_DIR}/BAMs/{wildcards.sample}.deduplicated.regular.ds.sorted.bam"
    output:
        methylation = "{run_dir}/METHYLATION/{sample}.deduplicated.regular.sorted.bismark.cov.gz"  if not DOWNSAMPLE else "{run_dir}/METHYLATION/{sample}.deduplicated.regular.ds.sorted.bismark.cov.gz",
        cov = "{run_dir}/METHYLATION/{sample}.deduplicated.regular.sorted.bismark.cov" if not DOWNSAMPLE else "{run_dir}/METHYLATION/{sample}.deduplicated.regular.ds.sorted.bismark.cov"
    params:
        run_dir = PROJECT_DIR,
        R1_5prime_ignore = f"--ignore {config['bismark_call']['R1_5prime_ignore']}" if config['bismark_call']['R1_5prime_ignore'] > 0 else "",
        R1_3prime_ignore = f"--ignore_3prime {config['bismark_call']['R1_3prime_ignore']}" if config['bismark_call']['R1_3prime_ignore'] > 0 else "",
        R2_5prime_ignore = f"--ignore_r2 {config['bismark_call']['R2_5prime_ignore']}" if config['bismark_call']['R2_5prime_ignore'] > 0 else "",
        R2_3prime_ignore = f"--ignore_3prime_r2 {config['bismark_call']['R2_3prime_ignore']}" if config['bismark_call']['R2_3prime_ignore'] > 0 else "",
        additional_params = CALL_METH_ARGS,
        samtools_cores = SAMTOOLS_CORES,
        suffix = lambda wildcards: f"{wildcards.sample}.deduplicated.regular.sorted" if not DOWNSAMPLE else f"{wildcards.sample}.deduplicated.regular.ds.sorted"
    threads: CALL_METH_CORES
    log:
        log = "{run_dir}/LOGS/METHYLATION_LOGS/{sample}_methylation.log"
    shell:
        """
        mkdir -p {params.run_dir}/TEMP/METHCALL
        samtools sort -@ {params.samtools_cores} -n -o {params.run_dir}/TEMP/METHCALL/{params.suffix}.bam {input.bam}
        bismark_methylation_extractor {params.additional_params} --paired-end --parallel {threads} --output_dir {params.run_dir}/METHYLATION --gzip --comprehensive --bedGraph {params.R1_5prime_ignore} {params.R1_3prime_ignore} {params.R2_5prime_ignore} {params.R2_3prime_ignore} {params.additional_params} {params.run_dir}/TEMP/METHCALL/{params.suffix}.bam > {log} 2>&1
        mv {params.run_dir}/METHYLATION/{wildcards.sample}*_report.txt {params.run_dir}/REPORTS/METHYLATION_REPORTS
        """
