import json
import sys
import argparse
import csv
import subprocess
from shared_functions import get_groups

def merge_bam_files(design_matrix, bam_files, output):
    unique_conditions, condition1_samples, condition2_samples = get_groups(design_matrix)
    print("unique_conditions", unique_conditions)
    print("condition1_samples", condition1_samples)
    print("condition2_samples", condition2_samples)
    condition1_bam_files = [bam_files[i] for i in range(len(bam_files)) if bam_files[i].split("/")[-1].split(".")[0] in condition1_samples]
    condition2_bam_files = [bam_files[i] for i in range(len(bam_files)) if bam_files[i].split("/")[-1].split(".")[0] in condition2_samples]
    condition1_bam_files = " ".join(condition1_bam_files)
    condition2_bam_files = " ".join(condition2_bam_files)
    result1 = subprocess.run(f"samtools merge -f -@ 8 {output}/TOBIAS/{unique_conditions[0]}.bam {condition1_bam_files}", shell=True)
    result2 = subprocess.run(f"samtools merge -f -@ 8 {output}/TOBIAS/{unique_conditions[1]}.bam {condition2_bam_files}", shell=True)
    if result1.returncode != 0 or result2.returncode != 0:
        sys.exit('\033[91m' + f"Error merging BAM files with return code {result1.returncode} and {result2.returncode}" + '\033[0m')
    print(f"Successfully merged BAM files for {unique_conditions[0]}: {condition1_bam_files}")
    print(f"Successfully merged BAM files for {unique_conditions[1]}: {condition2_bam_files}")
    return f"{output}/TOBIAS/{unique_conditions[0]}.bam", f"{output}/TOBIAS/{unique_conditions[1]}.bam"

def run_macs2(run_dir, merged_bam_files, genome, params):
    for file in merged_bam_files:
        subprocess.run(f"mkdir -p {run_dir}/TOBIAS/peak_calling", shell=True)
        result = subprocess.run(f"macs2 callpeak -t {file} --outdir {run_dir}/TOBIAS/peak_calling -n {file.split('/')[-1].split('.')[0]} -g {genome} {params}", shell=True)                                
        if result.returncode != 0:
            sys.exit('\033[91m' + f"Error running MACS2 with return code {result.returncode}" + '\033[0m')
        print(f"Successfully ran MACS2 for {file}")
    suffix = "broadPeak" if "broad" in params else "narrowPeak"
    return f"{run_dir}/TOBIAS/peak_calling/{merged_bam_files[0].split('/')[-1].split('.')[0]}_peaks.{suffix}", f"{run_dir}/TOBIAS/peak_calling/{merged_bam_files[1].split('/')[-1].split('.')[0]}_peaks.{suffix}"

def create_merge_bed_peaks(run_dir, bed_files):
    bed_files = " ".join(bed_files)
    result = subprocess.run(f"cat {bed_files} | bedtools sort | bedtools merge > {run_dir}/TOBIAS/merged_peaks.bed", shell=True)
    if result.returncode != 0:
        sys.exit('\033[91m' + f"Error merging BED files with return code {result.returncode}" + '\033[0m')
    print(f"Successfully merged BED files")

def create_uropa_input(run_dir, config_file):
    input_file = config_file
    output_file = f"{run_dir}/TOBIAS/uropa_config.json"
    data_to_add = {
        "bed": f"{run_dir}/TOBIAS/merged_peaks.bed",
        "prefix": "merged_peaks_annotated",
        "outdir": f"{run_dir}/TOBIAS",
        "output_by_query": "False",
        "show_attributes": ["all", ""],
        "priority": "False",
    }
    with open(input_file, 'r') as f:
        data = json.load(f)
    keys_to_remove = ["bed", "prefix", "outdir", "output_by_query", "show_attributes", "priority"]
    for key in keys_to_remove:
        data.pop(key, None)
    data.update(data_to_add)
    with open(output_file, 'w') as f:
        json.dump(data, f)
    print(f"Successfully created UROPA config file")
    
def run_uropa(run_dir, config):
    result = subprocess.run(f"uropa -i {config}", shell=True)
    if result.returncode != 0:
        sys.exit('\033[91m' + f"Error running UROPA with return code {result.returncode}" + '\033[0m')
    header_file = f"{run_dir}/TOBIAS/merged_peaks_annotated_header.txt"
    result2 = subprocess.run(f"head -n 1 {run_dir}/TOBIAS/merged_peaks_annotated_finalhits.txt > {header_file}", shell=True)
    if result2.returncode != 0:
        sys.exit('\033[91m' + f"Error creating header file with return code {result2.returncode}" + '\033[0m')
    print(f"Successfully ran UROPA")
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Create TOBIAS input files')
    parser.add_argument('-q', '--mapq', help='MAPQ values', required=True, nargs='+')
    parser.add_argument('-d', '--run_dir', help='Path to the run directory', required=True)
    parser.add_argument('-dm', '--design_matrix', help='Path to the design matrix', required=True)
    parser.add_argument('-b', '--bam_files', help='Path to the BAM files', required=True, nargs='+')
    parser.add_argument('-g', '--genome', help='Genome size', required=True)
    parser.add_argument('-p', '--params', help='MACS2 parameters', required=True)
    parser.add_argument('-c', '--config', help='UROPA config file', required=True)
    parser.add_argument('-m', '--motifs', help='Path to the motif file', required=True)
    parser.add_argument('-bl', '--blacklist', help='Path to the blacklist file', required=True)
    parser.add_argument('--cores', help='Number of cores to use', required=True)
    parser.add_argument('-gs', '--genome_size', help='Genome size, mm or hs', required=True)
    args = parser.parse_args()
    
    for q in args.mapq:
        files_q = [x for x in args.bam_files if f"mapq{q}" in x]
        files_q_rmdup = [x for x in files_q if "without_duplicates" in x]
        files_q_no_rmdup = [x for x in files_q if "with_duplicates" in x]  
        if files_q_rmdup != []:
            merged_bam_files = merge_bam_files(args.design_matrix, files_q_rmdup, args.run_dir)
            bed_files = run_macs2(args.run_dir, merged_bam_files, args.genome_size, args.params)
            create_merge_bed_peaks(args.run_dir, bed_files)
            create_uropa_input(args.run_dir, args.config)
            run_uropa(args.run_dir, f"{args.run_dir}/TOBIAS/uropa_config.json")
            merged_bam_files_list = [merged_bam_files[0], merged_bam_files[1]]
        if files_q_no_rmdup != []:
            merged_bam_files = merge_bam_files(args.design_matrix, files_q_no_rmdup, args.run_dir)
            bed_files = run_macs2(args.run_dir, merged_bam_files, args.genome_size, args.params)
            create_merge_bed_peaks(args.run_dir, bed_files)
            create_uropa_input(args.run_dir, args.config)
            run_uropa(args.run_dir, f"{args.run_dir}/TOBIAS/uropa_config.json")
            merged_bam_files_list = [merged_bam_files[0], merged_bam_files[1]]
    

