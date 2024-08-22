include: "common.smk"
path_to_genome_folder = os.path.dirname(GENOME_PATH)

rule align:
    input:
        R1 = lambda wildcards: f"{PROJECT_DIR}/TRIMMED_READS/{wildcards.sample}_R1{get_suffixes(wildcards.sample, RAW_DATA_DIR, trim = True)[0]}val_1.fq.gz" if TRIM else f"{RAW_DATA_DIR}/{wildcards.sample}_R1{get_suffixes(wildcards.sample, RAW_DATA_DIR)[0]}" + get_suffixes(wildcards.sample, RAW_DATA_DIR)[1],
        R2 = lambda wildcards: f"{PROJECT_DIR}/TRIMMED_READS/{wildcards.sample}_R2{get_suffixes(wildcards.sample, RAW_DATA_DIR, trim = True)[0]}val_2.fq.gz" if TRIM else f"{RAW_DATA_DIR}/{wildcards.sample}_R2{get_suffixes(wildcards.sample, RAW_DATA_DIR)[0]}" + get_suffixes(wildcards.sample, RAW_DATA_DIR)[1]
    output:
        bam = "{run_dir}/ALIGN_BAM/{sample}_aligned_pe.bam"
    params:
        genome = path_to_genome_folder,
        genome_fasta = GENOME_PATH,
        run_dir = PROJECT_DIR,
        additional_params = ALIGN_ARGS,
        mode = 0 if ALIGN_MODE == 'unique' else 1,
        samtools_cores = SAMTOOLS_CORES
    threads: ALIGN_CORES
    log:
        log = "{run_dir}/LOGS/ALIGN_LOGS/{sample}_bsmap.log"
    shell:
        """
        bsmap -p {threads} -r {params.mode} {params.additional_params} -a {input.R1} -b {input.R2} -d {params.genome_fasta} -o {params.run_dir}/TEMP/{wildcards.sample}_aligned_not_filtered.bam > {log} 2>&1
        samtools view -@ {params.samtools_cores} -b -f 2 -F 4 -o {output.bam} {params.run_dir}/TEMP/{wildcards.sample}_aligned_not_filtered.bam
        """
