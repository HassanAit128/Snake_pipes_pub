import csv
import re
import pandas as pd
import glob
import os
import time
from shared_functions import open_user_design_matrix
import subprocess
import concurrent.futures


def create_cmd_for_downsampling(files, dir, threads, memory):
    """
    Create a list containing the commands to downsample the BAM files
    
    Parameters
    ----------
        files : list
            List of BAM files to downsample
        quality : int
            Mapq value
        dir : str
            Path to the results directory
        threads : int
            Number of threads
        memory : int
            Amount of memory
        rmdup : bool
            If the duplicates were removed
        
    Returns
    -------
        ouptut_batch_file : list
            list of commands
    """
    commands = []
    counts = {}
    for file in files:
        count_file = f"{dir}/BAM/{file}_count.txt"
        with open(count_file, 'r') as f:
            count = f.read()
        counts[file] = int(count)
    min_count = min(counts.values())
    for file in counts.keys():
        fraction = min_count/counts[file]
        if fraction != 1:
            commands.append(f"samtools view -s {fraction} -b {dir}/BAM/{file}.bam -@ {threads} > {dir}/DS/{file}.downsampled.bam")
        else:
            commands.append(f"cp {dir}/BAM/{file}.bam {dir}/DS/{file}.downsampled.bam")
    return commands

def downsample_BAM_files(type, samples, mapq, dir, threads, memory, Des_Mat="", DS_Col="", DS_Grp=""):
    """
    Downsample BAM files

    Parameters
    ----------
        type : str
            Type of downsampling
        samples : list
            List of sample files
        mapq : list
            List of mapq values
        dir : str
            Directory path
        threads : int
            Number of threads
        memory : int
            Amount of memory
        Des_Mat : str, optional
            Path to the design matrix
        DS_Col : str, optional
            Column in the design matrix for downsampling
        DS_Grp : str, optional
            Group in the design matrix for downsampling

    Returns
    -------
        list
            List of all commands
    """
    all_commands = []
    wanted_output = []
    def process_files(files_q_rmdup, files_q_no_rmdup):
        nonlocal all_commands, wanted_output
        if files_q_rmdup:
            commands = create_cmd_for_downsampling(files_q_rmdup, dir, threads, memory)
            all_commands.extend(commands)
            wanted_output.extend([f"{dir}/DS/{x}.downsampled.bam" for x in files_q_rmdup])
        if files_q_no_rmdup:
            commands = create_cmd_for_downsampling(files_q_no_rmdup, dir, threads, memory)
            all_commands.extend(commands)
            wanted_output.extend([f"{dir}/DS/{x}.downsampled.bam" for x in files_q_no_rmdup])
    if type == "ALL":
        for q in mapq:
            files_q = [x for x in samples if f"mapq{q}" in x]
            files_q_rmdup = [x for x in files_q if "without_duplicates" in x]
            files_q_no_rmdup = [x for x in files_q if "with_duplicates" in x]
            process_files(files_q_rmdup, files_q_no_rmdup)
    else:
        design_matrix = open_user_design_matrix(Des_Mat)
        if DS_Grp != "": # downsample only the samples in the group specified
            for q in mapq:
                samples_to_downsample = [x for x in samples if design_matrix[x.split("/")[-1].split(".")[0]][DS_Col] == DS_Grp]
                files_q = [x for x in samples_to_downsample if f"mapq{q}" in x]
                files_q_rmdup = [x for x in files_q if "without_duplicates" in x]
                files_q_no_rmdup = [x for x in files_q if "with_duplicates" in x]
                process_files(files_q_rmdup, files_q_no_rmdup)
        else: # downsample all the samples by group of samples
            unique_grps = list({design_matrix[key][DS_Col] for key in design_matrix})
            for group in unique_grps:
                for q in mapq:
                    samples_to_downsample = [x for x in samples if design_matrix[x.split("/")[-1].split(".")[0]][DS_Col] == group]
                    files_q = [x for x in samples_to_downsample if f"mapq{q}" in x]
                    files_q_rmdup = [x for x in files_q if "without_duplicates" in x]
                    files_q_no_rmdup = [x for x in files_q if "with_duplicates" in x]
                    process_files(files_q_rmdup, files_q_no_rmdup)

    return all_commands

def execute_command(command):
    """Execute a single command using subprocess."""
    try:
        result = subprocess.run(command, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        return result.stdout
    except subprocess.CalledProcessError as e:
        return e.stderr

def calculate_optimal_resources(total_cores, total_memory, num_commands):
    max_workers = min(num_commands, total_cores)
    memory_per_worker = total_memory // max_workers
    cores_per_worker = max(1, total_cores // max_workers)    
    if cores_per_worker * max_workers > total_cores:
        max_workers = total_cores
        cores_per_worker = 1 
    return max_workers, memory_per_worker, cores_per_worker

def run_commands_in_parallel(num_threads, memory):
    """Run commands in parallel based on the specified number of threads."""
    max_workers, max_memory_per_worker, max_threads_per_worker = calculate_optimal_resources(num_threads, memory, len(snakemake.params.samples))
    all_commands = downsample_BAM_files(snakemake.params.ds_type, snakemake.params.samples, snakemake.params.mapq, snakemake.params.run_dir, max_threads_per_worker, max_memory_per_worker, snakemake.params.des_mat, snakemake.params.ds_col, snakemake.params.ds_grp)
    os.makedirs(f"{snakemake.params.run_dir}/DS", exist_ok=True)
    os.makedirs(f"{snakemake.params.run_dir}/LOGS", exist_ok=True)
    os.makedirs(f"{snakemake.params.run_dir}/LOGS/DS", exist_ok=True)
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
        results = executor.map(execute_command, all_commands)
        with open(f"{snakemake.params.run_dir}/LOGS/DS/downsample.log", 'w') as f:
            f.write("Running the following commands and their results:\n\n")
            for command, result in zip(all_commands, results):
                f.write(f"Command: {command}\n")
                f.write(f"Errors: {result.decode('utf-8')}\n" if result else "Errors: None\n\n")
                f.write("\n")
                        

run_commands_in_parallel(snakemake.params.cores, snakemake.params.memory)