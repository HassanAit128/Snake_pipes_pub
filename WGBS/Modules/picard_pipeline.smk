include: "common.smk"


rule dedup:
    input:
        bam = lambda wildcards: f"{PROJECT_DIR}/ALIGN_BAM/{wildcards.sample}_aligned_pe.bam"
    output:
        dedup_bam = "{run_dir}/BAMs/{sample}.picard.dedup.bam"
    params:
        run_dir = PROJECT_DIR,
        additional_params = DEDUP_ARGS,
        trim = TRIM
    log:
        log = "{run_dir}/LOGS/DEDUP_LOGS/{sample}_dedup.log"
    threads: SAMTOOLS_CORES
    shell:
        """
        samtools sort -@ {threads} -o {params.run_dir}/TEMP/{wildcards.sample}.sorted.bam {input}
        picard -Xmx8g MarkDuplicates \
            -I {params.run_dir}/TEMP/{wildcards.sample}.sorted.bam \
            -O {output.dedup_bam} \
            --REMOVE_DUPLICATES true \
            --METRICS_FILE {params.run_dir}/REPORTS/DEDUP_REPORTS/{wildcards.sample}.picard.dedup.txt > {log} 2>&1
        """