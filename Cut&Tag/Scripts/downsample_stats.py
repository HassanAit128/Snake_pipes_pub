import csv
import re
import pandas as pd
import glob
import os
from shared_functions import open_user_design_matrix
import subprocess
import time

def wait_for_files(files):
    """
    Wait for the files to be created
    
    Parameters
    ----------
    files : list
        List of files to wait for
    
    Returns
    -------
    None
    """
    all_files_exist = False
    while not all_files_exist:
        all_files_exist = all([os.path.exists(x) for x in files])
        time.sleep(1)
    return      

def count_reads(file, threads=1):
    """
    Count the number of reads in a BAM file
    
    Parameters
    ----------
    file : str
        Path to the BAM file
    
    Returns
    -------
    count : int
        Number of reads
    """
    command = f"samtools view -@ {threads} -c {file}"
    count = subprocess.run(command, capture_output=True, shell=True, text=True)
    return int(count.stdout.strip())

def create_downsample_stats_excel(BAM_files, directory, duplicates, mapq):
    dict = {}
    all_outputs = []
    for d in duplicates:
        for m in mapq:
            for f in BAM_files:
                if f"{d}_duplicates" in f and f"mapq{m}" in f:
                    dict[os.path.basename(f).split(".")[0]] = count_reads(f)
            df = pd.DataFrame.from_dict(dict, orient='index', columns=['Reads in downsampled BAM'])
            df.index.name = 'Sample'
            output = f"{directory}/stats/Summary/mapq{m}_{d}_duplicates_Downsample_stats.xlsx"
            all_outputs.append(output)
            df.to_excel(output)
    wait_for_files(all_outputs)
                    
            
create_downsample_stats_excel(snakemake.input, snakemake.params.run_dir, snakemake.params.dupli, snakemake.params.mapq)
