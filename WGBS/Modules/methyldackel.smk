include: "common.smk"


rule mbias_only:
    input:
        bam = lambda wildcards: f"{PROJECT_DIR}/BAMs/{wildcards.sample}.deduplicated.regular.sorted.bam" if not DOWNSAMPLE else f"{PROJECT_DIR}/BAMs/{wildcards.sample}.deduplicated.regular.ds.sorted.bam"
    output:
        mbias = "{run_dir}/METHYLATION/METHYLATION_BIAS/{sample}_OT.svg",
        mbias2 = "{run_dir}/METHYLATION/METHYLATION_BIAS/{sample}_OB.svg",
        log = "{run_dir}/LOGS/METHYLATION_LOGS/METHYLATION_BIAS_LOGS/{sample}_mbias.log"
    params:
        run_dir = PROJECT_DIR,
        additional_params = CALL_METH_ARGS,
        genome = GENOME_PATH,
        samtools_cores = SAMTOOLS_CORES,
        prefix = lambda wildcards: f"{wildcards.sample}.deduplicated.regular.sorted" if not DOWNSAMPLE else f"{wildcards.sample}.deduplicated.regular.ds.sorted"
    threads: CALL_MBIAS_ARGS
    shell:
        """
        MethylDackel mbias -@ {threads} {params.additional_params} {params.genome} {input.bam} {params.run_dir}/METHYLATION/METHYLATION_BIAS/{params.prefix} > {output.log} 2>&1
        """


def input_meth(sample):
    input = []
    if DOWNSAMPLE:
        input.append(f"{PROJECT_DIR}/BAMs/{sample}.deduplicated.regular.ds.sorted.bam")
    else:
        input.append(f"{PROJECT_DIR}/BAMs/{sample}.deduplicated.regular.sorted.bam")
    if AUTO_IGNORE:
        input.append(f"{PROJECT_DIR}/LOGS/METHYLATION_LOGS/METHYLATION_BIAS_LOGS/{sample}_mbias.log")
    
    
rule methylation_calling:
    input:
        bam = lambda wildcards: f"{PROJECT_DIR}/BAMs/{wildcards.sample}.deduplicated.regular.sorted.bam" if not DOWNSAMPLE else f"{PROJECT_DIR}/BAMs/{wildcards.sample}.deduplicated.regular.ds.sorted.bam",
    output:
        methylation = "{run_dir}/METHYLATION/{sample}.deduplicated.regular.sorted_CpG.bedGraph" if not DOWNSAMPLE else "{run_dir}/METHYLATION/{sample}.deduplicated.regular.ds.sorted_CpG.bedGraph",
        methylkit = "{run_dir}/METHYLATION/{sample}.deduplicated.regular.sorted_CpG.methylKit" if ALSO_METHYLKIT and not DOWNSAMPLE else "{run_dir}/METHYLATION/{sample}.deduplicated.regular.ds.sorted_CpG.methylKit" if ALSO_METHYLKIT and DOWNSAMPLE else ""
    params:
        run_dir = PROJECT_DIR,
        additional_params = CALL_METH_ARGS,
        genome = GENOME_PATH,
        prefix = lambda wildcards: f"{wildcards.sample}.deduplicated.regular.sorted" if not DOWNSAMPLE else f"{wildcards.sample}.deduplicated.regular.ds.sorted",
        ignore = f"--OT {config['Methyldackel']['OT_ignore']} --OB {config['Methyldackel']['OB_ignore']}" if not AUTO_IGNORE else "",
        also_methylkit = ALSO_METHYLKIT
    threads: CALL_METH_CORES
    log:
        log = "{run_dir}/LOGS/METHYLATION_LOGS/{sample}_methylation.log"
    shell:
        """
        if [ "{params.ignore}" == "" ]; then
            line=$(cat {params.run_dir}/LOGS/METHYLATION_LOGS/METHYLATION_BIAS_LOGS/{wildcards.sample}_mbias.log | grep "Suggested inclusion options")
            options=$(echo "$line" | awk -F': ' '{{print $2}}')
        else
            options="{params.ignore}"
        fi
        MethylDackel extract -@ {threads} {params.additional_params} $options -o {params.run_dir}/METHYLATION/{params.prefix} {params.genome} {input.bam} > {log} 2>&1
        if [ "{params.also_methylkit}" == "True" ]; then
            MethylDackel extract -@ {threads} --methylKit {params.additional_params} $options -o {params.run_dir}/METHYLATION/{params.prefix} {params.genome} {input.bam} > {log} 2>&1
        fi
        """