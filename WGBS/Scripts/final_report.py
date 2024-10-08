import os
import subprocess
import re
import glob
import pandas as pd

def extract_bsmap_stats(file):
    """
    Get alignement stats from BSMap log file.
    
    Parameters
    ----------
    file : str
        Path to the BSMap log file.
        
    Returns
    -------
    stats : dict
        Dictionary containing the alignment stats.
    """
    stats = {}
    with open(file, "r") as f :
        data = f.read()
    stats['SampleID'] = re.search(r'\/([^\/]+?)_[^\/]+\.bam', data).group(1)
    stats['Total reads'] = int(re.search(r'total read pairs: (\d+)', data).group(1))*2
    stats['Aligned reads'] = int(re.search(r'aligned pairs: (\d+)', data).group(1))*2
    stats['Aligned reads %'] = float(re.search(r'aligned pairs: \d+ \((\d+\.\d+)%\)', data).group(1)) 
    stats['Unique reads'] = int(re.search(r'unique pairs: (\d+)', data).group(1))*2
    stats['Non unique reads'] = int(re.search(r'non-unique pairs: (\d+)', data).group(1))*2
    stats['Total time consumed (h)'] = int(re.search(r'total time consumed:  (\d+)', data).group(1))/3600
    return stats

def extract_bismark_stats(file):
    """
    Extract alignment stats from Bismark report file.

    Parameters
    ----------
    file : str
        Path to the Bismark report file.
        
    Returns
    -------
    stats : dict
        Dictionary containing the alignment stats.
    """
    stats = {}
    with open(file, "r") as f :
        data = f.read()
    stats['SampleID'] = re.search(r'\/(SRR\d+)_', data).group(1)
    stats['Total read pairs'] = int(re.search(r'Sequence pairs analysed in total:\t(\d+)', data).group(1))
    stats['Unique pairs'] = int(re.search(r'Number of paired-end alignments with a unique best hit:\t(\d+)', data).group(1))
    stats['Unique pairs %'] = float(re.search(r'Mapping efficiency:\t(\d+\.\d+)%', data).group(1))
    stats['Non unique pairs'] = int(re.search(r'Sequence pairs did not map uniquely:\t(\d+)', data).group(1))
    stats['Total time consumed (h)'] = int(re.search(r'Bismark completed in (\d+)d (\d+)h (\d+)m (\d+)s', data).group(1)) * 24 + int(re.search(r'Bismark completed in (\d+)d (\d+)h (\d+)m (\d+)s', data).group(2)) + int(re.search(r'Bismark completed in (\d+)d (\d+)h (\d+)m (\d+)s', data).group(3))/60 + int(re.search(r'Bismark completed in (\d+)d (\d+)h (\d+)m (\d+)s', data).group(4))/3600
    return stats

def create_align_reports(run_dir, run_name, tool):
    """
    Create alignment reports for each sample by calling the appropriate function depending on the alignment tool used.
    
    Parameters
    ----------
    run_dir : str
        Path to the run directory.
    run_name : str
        Name of the run.
    tool : str
        Alignment tool used (BSMap or Bismark).
        
    Returns
    -------
    stats : list
        List of dictionaries containing the alignment stats for each sample.
    """
    stats = []
    if tool.lower() == "bsmap":
        print("--- Extracting BSMap stats ---")
        align_reports = glob.glob(f"{run_dir}/{run_name}/LOGS/ALIGN_LOGS/*.log")
        for file in align_reports:
            stats.append(extract_bsmap_stats(file))
    elif tool.lower() == "bismark":
        print("--- Extracting Bismark stats ---")
        align_reports = glob.glob(f"{run_dir}/{run_name}/REPORTS/ALIGN_REPORTS/*.txt")
        for file in align_reports:
            stats.append(extract_bismark_stats(file))
    return stats
        
def samtools_count(bam_file, cores):
    """
    Count the number of reads in a BAM file using samtools.
    
    Parameters
    ----------
    bam_file : str
        Path to the BAM file.
    cores : int
        Number of cores to use.
    
    Returns
    -------
    count : int
        Number of reads in the BAM file.
    """
    cmd = f"samtools view -@ {cores} -c {bam_file}"
    count = subprocess.check_output(cmd, shell=True).decode().strip()
    return count 

def get_bam_stats(samples, run_dir, run_name, cores, downsampling):
    """
    Count the number of reads in the BAM files for each sample and for each step of the pipeline.
    
    Parameters
    ----------
    samples : list
        List of sample IDs.
    run_dir : str
        Path to the run directory.
    run_name : str
        Name of the run.
    cores : int
        Number of cores to use when calling samtools.
    downsampling : bool
        Whether downsampling was performed or not.
        
    Returns
    -------
    stats : list
        List of dictionaries containing the alignment stats for each sample.
    """
    stats = []
    for sample in samples:
        print(f"--- Getting stats for {sample} ---")
        dedplicated_count = samtools_count(f"{run_dir}/{run_name}/BAMs/{sample}.deduplicated.bam", cores)
        regular_count = samtools_count(f"{run_dir}/{run_name}/BAMs/{sample}.deduplicated.regular.bam", cores)
        final_count = samtools_count(f"{run_dir}/{run_name}/BAMs/{sample}.deduplicated.regular.sorted.bam", cores)
        if downsampling:
            ds_count = samtools_count(f"{run_dir}/{run_name}/BAMs/{sample}.deduplicated.regular.ds.sorted.bam", cores)
            stats.append({"SampleID": sample, "Deduplicated": dedplicated_count, "Regular Chr": regular_count, "Final_Count": final_count, "Final Count (million)": int(final_count)/1000000, "Downsampled": ds_count})
        else:
            stats.append({"SampleID": sample, "Deduplicated": dedplicated_count, "Regular Chr": regular_count, "Final_Count": final_count, "Final Count (million)": int(final_count)/1000000})
        print(f"--- Done for {sample} ---")
    return stats

def get_merge_and_save_stats(design_matrix, run_dir, run_name, cores, tool, downsampling):
    """
    'Main' function to get alignment stats and BAM stats, merge them with the design matrix and save the final report.
    
    Parameters
    ----------
    design_matrix : str
        Path to the design matrix.
    run_dir : str
        Path to the run directory.
    run_name : str
        Name of the run.
    cores : int
        Number of cores to use when calling samtools.
    tool : str
        Alignment tool used (BSMap or Bismark).
    downsampling : bool
        Whether downsampling was performed or not.
    
    Returns
    -------
    merged_df : pandas.DataFrame
        DataFrame containing the final report.
    """
    print("--- Reading design matrix ---")
    design_matrix = pd.read_csv(design_matrix)
    print("--- OK ---")
    print("--- Creating alignment reports ---")
    align_stats = create_align_reports(run_dir, run_name, tool)
    print("--- OK ---")
    print("--- Getting BAM stats ---")
    samples = [sample["SampleID"] for sample in align_stats]
    print(samples)
    bam_stats = get_bam_stats(samples, run_dir, run_name, cores, downsampling)
    print("--- OK ---")
    print("--- Merging and saving stats ---")
    align_stats_df = pd.DataFrame(align_stats)
    bam_stats_df = pd.DataFrame(bam_stats)
    merged_df = pd.merge(design_matrix, align_stats_df, on="SampleID")
    merged_df = pd.merge(merged_df, bam_stats_df, on="SampleID")
    merged_df.to_excel(f"{run_dir}/{run_name}/FINAL_REPORT/Stats.xlsx", index=False)
    return merged_df