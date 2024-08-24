include: "Modules/common.smk"

path_to_genome_folder = os.path.dirname(GENOME_PATH)

############################################# HELPER FUNCTIONS #############################################
def set_up_directories():
    """
    Create the directories needed for the pipeline
    """
    if ALIGN:
        if not os.path.exists(PROJECT_DIR + '/ALIGN_BAM'):
            os.makedirs(PROJECT_DIR + '/ALIGN_BAM')
            os.makedirs(PROJECT_DIR + '/REPORTS/ALIGN_REPORTS', exist_ok=True)
            os.makedirs(PROJECT_DIR + '/LOGS/ALIGN_LOGS', exist_ok=True)
        if ALIGN_TOOL == 'bismark':
            if IDX_GENOME:
                if not os.path.exists(PROJECT_DIR + '/IDX_GENOME'):
                    os.makedirs(PROJECT_DIR + '/IDX_GENOME')
                    os.makedirs(PROJECT_DIR + '/LOGS/IDX_LOGS', exist_ok=True)
    if TRIM:
        if not os.path.exists(PROJECT_DIR + '/TRIMMED_READS'):
            os.makedirs(PROJECT_DIR + '/TRIMMED_READS', exist_ok=True)
            os.makedirs(PROJECT_DIR + '/LOGS/TRIM_LOGS', exist_ok=True)
            os.makedirs(PROJECT_DIR + '/REPORTS/TRIM_REPORTS', exist_ok=True)
    if DEDUP:
        if not os.path.exists(PROJECT_DIR + '/BAMs'):
            os.makedirs(PROJECT_DIR + '/BAMs', exist_ok=True)
            os.makedirs(PROJECT_DIR + '/REPORTS/DEDUP_REPORTS', exist_ok=True)
            os.makedirs(PROJECT_DIR + '/LOGS/DEDUP_LOGS', exist_ok=True)
    if CALL_METHYLATION:
        if not os.path.exists(PROJECT_DIR + '/METHYLATION'):
            os.makedirs(PROJECT_DIR + '/METHYLATION', exist_ok=True)
            os.makedirs(PROJECT_DIR + '/REPORTS/METHYLATION_REPORTS', exist_ok=True)
            os.makedirs(PROJECT_DIR + '/LOGS/METHYLATION_LOGS', exist_ok=True)
    if CALL_METHYLATION_BIAS:
        if not os.path.exists(PROJECT_DIR + '/METHYLATION/METHYLATION_BIAS'):
            os.makedirs(PROJECT_DIR + '/METHYLATION/METHYLATION_BIAS', exist_ok=True)
            os.makedirs(PROJECT_DIR + '/REPORTS/METHYLATION_REPORTS/METHYLATION_BIAS_REPORTS', exist_ok=True)
            os.makedirs(PROJECT_DIR + '/LOGS/METHYLATION_LOGS/METHYLATION_BIAS_LOGS', exist_ok=True)
        if DMR_analysis:
            if not os.path.exists(PROJECT_DIR + '/METHYLATION/DMR'):
                os.makedirs(PROJECT_DIR + '/METHYLATION/DMR', exist_ok=True)
                os.makedirs(PROJECT_DIR + '/REPORTS/METHYLATION_REPORTS/DMR_REPORTS', exist_ok=True)
                os.makedirs(PROJECT_DIR + '/LOGS/METHYLATION_LOGS/DMR_LOGS', exist_ok=True)
    if GENERATE_REPORT:
        if not os.path.exists(PROJECT_DIR + '/FINAL_REPORT'):
            os.makedirs(PROJECT_DIR + '/FINAL_REPORT', exist_ok=True)
    if not os.path.exists(PROJECT_DIR + '/TEMP'):
        os.makedirs(PROJECT_DIR + '/TEMP', exist_ok=True)

idx_output = []

def set_up_output_files():
    all_output = []
    if TRIM:
        for sample in SAMPLES:
            all_output.append(PROJECT_DIR + '/TRIMMED_READS/' + sample + "_R1" + get_suffixes(sample, RAW_DATA_DIR, trim = True)[0] + 'val_1.fq.gz')
            all_output.append(PROJECT_DIR + '/TRIMMED_READS/' + sample + "_R2" + get_suffixes(sample, RAW_DATA_DIR, trim = True)[0] + 'val_2.fq.gz')
    if ALIGN:
        if ALIGN_TOOL == "bismark":
            if IDX_GENOME:
                idx_output = []
                for conv in ["CT", "GA"]:
                    idx_output.append(os.path.dirname(GENOME_PATH) + '/Bisulfite_Genome/' + conv + '_conversion/' + 'BS_' + conv + '.1.bt2')
                    idx_output.append(os.path.dirname(GENOME_PATH) + '/Bisulfite_Genome/' + conv + '_conversion/' + 'BS_' + conv + '.2.bt2')
                    idx_output.append(os.path.dirname(GENOME_PATH) + '/Bisulfite_Genome/' + conv + '_conversion/' + 'BS_' + conv + '.3.bt2')
                    idx_output.append(os.path.dirname(GENOME_PATH) + '/Bisulfite_Genome/' + conv + '_conversion/' + 'BS_' + conv + '.4.bt2')
                    idx_output.append(os.path.dirname(GENOME_PATH) + '/Bisulfite_Genome/' + conv + '_conversion/' + 'BS_' + conv + '.rev.1.bt2')
                    idx_output.append(os.path.dirname(GENOME_PATH) + '/Bisulfite_Genome/' + conv + '_conversion/' + 'BS_' + conv + '.rev.2.bt2')
                all_output.extend(idx_output)
            align_output = [PROJECT_DIR + '/ALIGN_BAM/' + sample + '_aligned_pe.bam' for sample in SAMPLES]
            all_output.extend(align_output)
        if ALIGN_TOOL == 'bsmap':
            pass
    if DEDUP:
        dedup_output = [PROJECT_DIR + '/BAMs/' + sample + '.deduplicated.bam' for sample in SAMPLES]
        all_output.extend(dedup_output)
    if KEEP_REGULAR:
        regular_output = [PROJECT_DIR + '/BAMs/' + sample + '.deduplicated.regular.bam' for sample in SAMPLES]
        all_output.extend(regular_output)
    if DOWNSAMPLE:
        ds_output = [PROJECT_DIR + '/BAMs/' + sample + '.deduplicated.regular.ds.sorted.bam' for sample in SAMPLES]
        all_output.extend(ds_output)
    prefix = ".deduplicated.regular.sorted" if not DOWNSAMPLE else ".deduplicated.regular.ds.sorted"
    if CALL_METHYLATION_BIAS:
        suffix = ".M-bias.txt" if CALL_METHYLATION_TOOL == 'bismark' else "_OT.svg"
        mbias_output = [PROJECT_DIR + '/METHYLATION/METHYLATION_BIAS/' + sample + prefix + suffix for sample in SAMPLES]
        if CALL_METHYLATION_TOOL == 'bismark':
            mbias_output.extend([PROJECT_DIR + '/METHYLATION/METHYLATION_BIAS/' + sample + prefix + ".M-bias_R1.png" for sample in SAMPLES])
        all_output.extend(mbias_output)
        if STOP_ALL and not AUTO_IGNORE:
            return all_output
    if CALL_METHYLATION:
        suffix = ".bedGraph.gz" if CALL_METHYLATION_TOOL == 'bismark' else "_CpG.bedGraph"
        methylation_output = [PROJECT_DIR + '/METHYLATION/' + sample + prefix + suffix for sample in SAMPLES]
        all_output.extend(methylation_output)
        if DMR_analysis:
            if DMR_analysis_tool.lower() == 'metilene':
                metilene_output = f"{PROJECT_DIR}/METHYLATION/DMR/metilene_DMRs.txt"
                all_output.append(metilene_output)
            elif DMR_analysis_tool.lower() == 'methylkit':
                methylkit_output = f"{PROJECT_DIR}/METHYLATION/DMR/{config['name']}_methylkit_analysis.RData"
                all_output.append(methylkit_output)
    return all_output
############################################################################################################

set_up_directories()
output_all = set_up_output_files()

############################################# RULES ########################################################
rule ALL:
    input:
        output_all

rule trim_galore:
    input:
        R1 = lambda wildcards: f"{RAW_DATA_DIR}/{wildcards.sample}_R1{get_suffixes(wildcards.sample, RAW_DATA_DIR)[0]}" + get_suffixes(wildcards.sample, RAW_DATA_DIR)[1],
        R2 = lambda wildcards: f"{RAW_DATA_DIR}/{wildcards.sample}_R2{get_suffixes(wildcards.sample, RAW_DATA_DIR)[0]}" + get_suffixes(wildcards.sample, RAW_DATA_DIR)[1]
    output:
        R1 = "{run_dir}/TRIMMED_READS/{sample}_R1{suffix}val_1.fq.gz",
        R2 = "{run_dir}/TRIMMED_READS/{sample}_R2{suffix}val_2.fq.gz",
        fasqc_R1 = "{run_dir}/REPORTS/TRIM_REPORTS/fastqc_reports/{sample}_R1{suffix}val_1_fastqc.html",
        fasqc_R2 = "{run_dir}/REPORTS/TRIM_REPORTS/fastqc_reports/{sample}_R2{suffix}val_2_fastqc.html"
    params:
        run_dir = PROJECT_DIR,
        suffix = lambda wildcards: get_suffixes(wildcards.sample, RAW_DATA_DIR, trim = True)[0],
        additional_params = TRIM_ARGS
    log:
        Log = "{run_dir}/LOGS/TRIM_LOGS/{sample}_{suffix}.log"
    threads: TRIM_CORES if TRIM_CORES <= 4 else 4
    singularity: 
        "docker://hasba/main_d"
    shell:
        """
        trim_galore --paired {params.additional_params} --gzip --cores {threads} --output_dir {params.run_dir}/TRIMMED_READS/ --fastqc --fastqc_args "--outdir {params.run_dir}/REPORTS/TRIM_REPORTS/fastqc_reports/" {input.R1} {input.R2} > {log.Log} 2>&1
        mv {params.run_dir}/TRIMMED_READS/{wildcards.sample}*_report.txt {params.run_dir}/REPORTS/TRIM_REPORTS
        """

module alignment_module:
    snakefile:
        "Modules/bismark_pipeline.smk" if ALIGN_TOOL == 'bismark' else "Modules/bsmap_pipeline.smk"
    config:
        config

use rule align from alignment_module as align_rule with:
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
        mode = 0 if ALIGN_TOOL == 'bsmap' and ALIGN_MODE == 'unique' else 1,
        samtools_cores = SAMTOOLS_CORES
    threads: ALIGN_CORES
    singularity: 
        "docker://hasba/main_d"
    log:
        log = "{run_dir}/LOGS/ALIGN_LOGS/{sample}_bismark.log"

module deduplication_module:
    snakefile:
        "Modules/bismark_pipeline.smk" if DEDUP_TOOL == 'bismark' else "Modules/picard_pipeline.smk"
    config:
        config

use rule dedup from deduplication_module as dedup_rule with:
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
    singularity: 
        "docker://hasba/main_d"
    threads: SAMTOOLS_CORES
    
rule keep_regular:
    input:
        bam = lambda wildcards: f"{PROJECT_DIR}/BAMs/{wildcards.sample}.deduplicated.bam" 
    output:
        bam = "{run_dir}/BAMs/{sample}.deduplicated.regular.bam"
    params:
        run_dir = PROJECT_DIR,
        genome = GENOME
    log:
        log = "{run_dir}/LOGS/SAMTOOLS/{sample}_keep_regular.log"
    threads: SAMTOOLS_CORES
    singularity: 
        "docker://hasba/main_d"
    shell:
        """
        samtools sort -@ {threads} -o {params.run_dir}/TEMP/{wildcards.sample}.deduplicated.sorted.bam {input}
        samtools index -@ {threads} {params.run_dir}/TEMP/{wildcards.sample}.deduplicated.sorted.bam
        if [ {params.genome} = hg38 ]
        then
            samtools view --threads {threads} -hb {params.run_dir}/TEMP/{wildcards.sample}.deduplicated.sorted.bam chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY > {params.run_dir}/TEMP/{wildcards.sample}.deduplicated.regular.bam 2> {log}
        elif [ {params.genome} = mm10 ]
        then
            samtools view --threads {threads} -hb {params.run_dir}/TEMP/{wildcards.sample}.deduplicated.sorted.bam chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX chrY > {params.run_dir}/TEMP/{wildcards.sample}.deduplicated.regular.bam 2> {log}
            samtools sort -n -@ {threads} -o {output} {params.run_dir}/TEMP/{wildcards.sample}.deduplicated.regular.bam
        fi
        """

rule sort:
    input:
        bam = lambda wildcards: f"{PROJECT_DIR}/BAMs/{wildcards.sample}.deduplicated.regular.bam"
    output:
        bam = "{run_dir}/BAMs/{sample}.deduplicated.regular.sorted.bam"
    params:
        run_dir = PROJECT_DIR,
        sort_type = "-n" if not DOWNSAMPLE else ""
    log:
        log = "{run_dir}/LOGS/SAMTOOLS/{sample}_sort.log"
    threads: SAMTOOLS_CORES
    singularity: 
        "docker://hasba/main_d"
    shell:
        """
        samtools sort {params.sort_type} -@ {threads} -o {output.bam} {input.bam} 2> {log}
        """
        
rule find_read_count:
    input:
        lambda wildcards: f"{PROJECT_DIR}/BAMs/{wildcards.sample}.deduplicated.regular.sorted.bam"
    output:
        "{run_dir}/TEMP/{sample}.read_count.txt"
    threads: SAMTOOLS_CORES
    params:
        run_dir = PROJECT_DIR
    singularity: 
        "docker://hasba/main_d"
    shell:
        """
        samtools view --threads {threads} -c {input} > {output}
        """

rule get_lowest_count:
    input:
        read_counts = expand("{run_dir}/TEMP/{sample}.read_count.txt", run_dir = PROJECT_DIR, sample = SAMPLES)
    output:
        "{run_dir}/TEMP/lowest_count.txt"
    shell:
        """
        cat {input} | sort -n | head -n 1 > {output}
        """

rule downsample:
    input:
        bam = lambda wildcards: f"{PROJECT_DIR}/BAMs/{wildcards.sample}.deduplicated.regular.sorted.bam",
        lowest_count = "{run_dir}/TEMP/lowest_count.txt"
    output:
        bam = "{run_dir}/BAMs/{sample}.deduplicated.regular.ds.sorted.bam"
    params:
        run_dir = PROJECT_DIR,
    threads: SAMTOOLS_CORES
    log:
        log = "{run_dir}/LOGS/SAMTOOLS/{sample}_downsample.log"
    singularity: 
        "docker://hasba/main_d"
    shell:
        """
        samtools index -@ {threads} {input.bam} 2> {log}
        lowest_count=$(cat {input.lowest_count})
        frac=$( samtools idxstats {input.bam} | cut -f3 | awk -v lowest_count=$lowest_count 'BEGIN {{total=0}} {{total += $1}} END {{frac=lowest_count/total; if (frac > 1) {{print 1}} else {{print frac}}}}' )
        if [ $(echo "$frac 1.0" | awk '{{if ($1 >= $2) print 1; else print 0}}') -eq 1 ]; then
            cp {input.bam} "{params.run_dir}/TEMP/{wildcards.sample}.deduplicated.regular.ds.bam" 2> {log}
        else
            samtools view -s $frac -b {input.bam} -@ 12 > "{params.run_dir}/TEMP/{wildcards.sample}.deduplicated.regular.ds.bam" 2> {log}
        fi
        samtools sort -@ {threads} -o {output.bam} "{params.run_dir}/TEMP/{wildcards.sample}.deduplicated.regular.ds.bam" 2> {log}
        """

module methylation_module:
    snakefile:
        "Modules/bismark_pipeline.smk" if CALL_METHYLATION_TOOL == 'bismark' else "Modules/methyldackel.smk"
    config:
        config

if CALL_METHYLATION_TOOL == "bismark":
    use rule mbias_only from methylation_module as mbias_only_rule with:
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
        singularity: 
            "docker://hasba/main_d"
        log:
            log = "{run_dir}/LOGS/METHYLATION_LOGS/METHYLATION_BIAS_LOGS/{sample}_mbias.log"

    use rule methylation_calling from methylation_module as methylation_calling_rule with:
        input:
            bam = lambda wildcards: f"{PROJECT_DIR}/BAMs/{wildcards.sample}.deduplicated.regular.sorted.bam" if not DOWNSAMPLE else f"{PROJECT_DIR}/BAMs/{wildcards.sample}.deduplicated.regular.ds.sorted.bam"
        output:
            methylation = "{run_dir}/METHYLATION/{sample}.deduplicated.regular.sorted.bedGraph.gz"  if not DOWNSAMPLE else "{run_dir}/METHYLATION/{sample}.deduplicated.regular.ds.sorted.bedGraph.gz",
            cov = "{run_dir}/METHYLATION/{sample}.deduplicated.regular.sorted.bismark.cov.gz" if not DOWNSAMPLE else "{run_dir}/METHYLATION/{sample}.deduplicated.regular.ds.sorted.bismark.cov.gz"
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
        singularity: 
            "docker://hasba/main_d"
        log:
            log = "{run_dir}/LOGS/METHYLATION_LOGS/{sample}_methylation.log"

    use rule draw_mbias_plot from methylation_module as draw_mbias_plot_rule with:
        input:
            mbias = lambda wildcards: f"{PROJECT_DIR}/METHYLATION/METHYLATION_BIAS/{wildcards.sample}.deduplicated.regular.sorted.M-bias.txt" if not DOWNSAMPLE else f"{PROJECT_DIR}/METHYLATION/METHYLATION_BIAS/{wildcards.sample}.deduplicated.regular.ds.sorted.M-bias.txt"
        output:
            plot = "{run_dir}/METHYLATION/METHYLATION_BIAS/{sample}.deduplicated.regular.sorted.M-bias_R1.png" if not DOWNSAMPLE else "{run_dir}/METHYLATION/METHYLATION_BIAS/{sample}.deduplicated.regular.ds.sorted.M-bias_R1.png",
            plot2 = "{run_dir}/METHYLATION/METHYLATION_BIAS/{sample}.deduplicated.regular.sorted.M-bias_R2.png" if not DOWNSAMPLE else "{run_dir}/METHYLATION/METHYLATION_BIAS/{sample}.deduplicated.regular.ds.sorted.M-bias_R2.png"
        params:
            run_dir = PROJECT_DIR
        log:    
            log = "{run_dir}/LOGS/METHYLATION_LOGS/METHYLATION_BIAS_LOGS/{sample}_mbias_plot.log"

def input_meth(sample):
    input = []
    if DOWNSAMPLE:
        input.append(f"{PROJECT_DIR}/BAMs/{sample}.deduplicated.regular.ds.sorted.bam")
    else:
        input.append(f"{PROJECT_DIR}/BAMs/{sample}.deduplicated.regular.sorted.bam")
    if AUTO_IGNORE:
        input.append(f"{PROJECT_DIR}/LOGS/METHYLATION_LOGS/METHYLATION_BIAS_LOGS/{sample}_mbias.log")
    return input

if CALL_METHYLATION_TOOL == "methyldackel":
    use rule mbias_only from methylation_module as mbias_only_rule with:
        input:
            bam = lambda wildcards: f"{PROJECT_DIR}/BAMs/{wildcards.sample}.deduplicated.regular.sorted.bam" if not DOWNSAMPLE else f"{PROJECT_DIR}/BAMs/{wildcards.sample}.deduplicated.regular.ds.sorted.bam"
        output:
            mbias = "{run_dir}/METHYLATION/METHYLATION_BIAS/{sample}.deduplicated.regular.sorted_OT.svg" if not DOWNSAMPLE else "{run_dir}/METHYLATION/METHYLATION_BIAS/{sample}.deduplicated.regular.ds.sorted_OT.svg",
            log = "{run_dir}/LOGS/METHYLATION_LOGS/METHYLATION_BIAS_LOGS/{sample}_mbias.log"
        params:
            run_dir = PROJECT_DIR,
            additional_params = CALL_METH_ARGS,
            genome = GENOME_PATH,
            samtools_cores = SAMTOOLS_CORES,
            prefix = lambda wildcards: f"{wildcards.sample}.deduplicated.regular.sorted" if not DOWNSAMPLE else f"{wildcards.sample}.deduplicated.regular.ds.sorted"
        singularity: 
            "docker://hasba/main_d"
        threads: CALL_METH_CORES

    use rule methylation_calling from methylation_module as methylation_calling_rule with:
        input:
            bam = lambda wildcards: input_meth(wildcards.sample)
        output:
            methylation = ["{run_dir}/METHYLATION/{sample}.deduplicated.regular.sorted_CpG.bedGraph" if not DOWNSAMPLE else "{run_dir}/METHYLATION/{sample}.deduplicated.regular.ds.sorted_CpG.bedGraph"] + ["{run_dir}/METHYLATION/{sample}.deduplicated.regular.sorted_CpG.methylKit" if not DOWNSAMPLE and ALSO_METHYLKIT else "{run_dir}/METHYLATION/{sample}.deduplicated.regular.ds.sorted_CpG.methylKit" if DOWNSAMPLE and ALSO_METHYLKIT else ""]
        params:
            run_dir = PROJECT_DIR,
            additional_params = CALL_METH_ARGS,
            genome = GENOME_PATH,
            prefix = lambda wildcards: f"{wildcards.sample}.deduplicated.regular.sorted" if not DOWNSAMPLE else f"{wildcards.sample}.deduplicated.regular.ds.sorted",
            also_methylkit = ALSO_METHYLKIT,
            ignore = f"--OT {config['Methyldackel']['OT_ignore']} --OB {config['Methyldackel']['OB_ignore']}" if not AUTO_IGNORE else ""
        threads: CALL_METH_CORES
        singularity: 
            "docker://hasba/main_d"
        log:
            log = "{run_dir}/LOGS/METHYLATION_LOGS/{sample}_methylation.log"

DS = 'ds.' if DOWNSAMPLE else ''

rule sort_for_metilene:
    input: 
        call = lambda wildcards: f"{PROJECT_DIR}/METHYLATION/{wildcards.sample}.deduplicated.regular.{DS}sorted.bedGraph.gz" if CALL_METHYLATION_TOOL == 'bismark' else f"{PROJECT_DIR}/METHYLATION/{wildcards.sample}.deduplicated.regular.{DS}sorted_CpG.bedGraph"
    output:
        "{run_dir}/METHYLATION/{sample}.sorted_mod.bedGraph"
    params:
        run_dir = PROJECT_DIR,
        coverage = config["DMR_analysis"]["minimum_coverage"]
    shell:
        """
        inputfile={input}
        if [[ "{input}" == *.gz ]]; then
            gzip -d -k {input}
            inputfile=$(echo {input} | sed 's/.gz//')
        fi
        sort -i $inputfile > {params.run_dir}/TEMP/{wildcards.sample}.sorted.bedGraph
        awk '{{print $1,$2,$3,$4}}' {params.run_dir}/TEMP/{wildcards.sample}.sorted.bedGraph > {params.run_dir}/TEMP/{wildcards.sample}.sorted_m.bedGraph
        cat {params.run_dir}/TEMP/{wildcards.sample}.sorted_m.bedGraph | tr ' ' '\t' > {output}
        """ 

rule run_metilene:
    input:
        expand("{PROJECT_DIR}/METHYLATION/{sample}.sorted_mod.bedGraph", PROJECT_DIR = PROJECT_DIR, sample = SAMPLES)
    output:
        "{run_dir}/METHYLATION/DMR/metilene_DMRs.txt"
    params:
        run_dir = PROJECT_DIR,
        path_to_metilene = config['metilene']['path_to_metilene'],
        matrix = DESIGN_MATRIX,
        samples = SAMPLES,
        maxdist = config['metilene']['maxdist'],
        mincpgs = config['metilene']['mincpgs'],
        minMethDiff = config['metilene']['minMethDiff'],
        additional_params = config['metilene']['metilene_args'] if config['metilene']['metilene_args'] != "" else "None"
    threads: config['metilene']['metilene_cores']
    shell:
        """
        python Scripts/metilene_prep.py -m {params.matrix} -d {params.run_dir} -mp {params.path_to_metilene} -s {params.samples} -sp {input} -o {output} -c {threads} -md {params.maxdist} -mc {params.mincpgs} -mdf {params.minMethDiff} -aa {params.additional_params}
        """

rule run_methylkit:
    input:
        calls = expand("{run_dir}/METHYLATION/{sample}.deduplicated.regular.{ds}sorted.bismark.cov.gz", run_dir = PROJECT_DIR, sample = SAMPLES, ds = DS) if CALL_METHYLATION_TOOL == 'bismark' else expand("{run_dir}/METHYLATION/{sample}.deduplicated.regular.{ds}sorted_CpG.methylKit", run_dir = PROJECT_DIR, sample = SAMPLES, ds = DS)
    output:
        "{run_dir}/METHYLATION/DMR/{run}_methylkit_analysis.RData"
    params:
        run_dir = PROJECT_DIR,
        genome = GENOME,
        pp = "bismark" if CALL_METHYLATION_TOOL == 'bismark' else "amp",
        mincov = config['DMR_analysis']['minimum_coverage'],
        min_methylation_diff = config['methylkit']['min_methylation_diff'],
        pvalue = config['methylkit']['pvalue'],
        min_DMR_length = config['methylkit']['min_DMR_length'],
        regions = config['methylkit']['regions'],
        design_matrix = DESIGN_MATRIX,
        run = config['name']
    singularity: 
        "docker://hasba/mekit"
    script:
        "Scripts/methylkit.R"

############################################################################################################