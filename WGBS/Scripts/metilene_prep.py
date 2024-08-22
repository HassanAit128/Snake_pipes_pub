import csv
import argparse
import subprocess
import sys

def get_groups(design_matrix):
    with open(design_matrix, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        conditions = []
        samples = []
        for row in reader:
            conditions.append(row["Condition"])
            samples.append(row["SampleID"])
        unique_conditions = list(set(conditions))
        if "wt" in unique_conditions[0].lower():
            unique_conditions = unique_conditions[::-1]
        condition1_samples = [samples[i] for i in range(len(samples)) if conditions[i] == unique_conditions[0]]
        condition2_samples = [samples[i] for i in range(len(samples)) if conditions[i] == unique_conditions[1]]
    return unique_conditions, condition1_samples, condition2_samples

def create_input_for_metilene(run_dir, path_to_metilene, conditions, condition1_samples, condition2_samples, samples, paths):
    out = f"{run_dir}/METHYLATION/DMR/input_for_metilene.txt"
    condition1_paths = [paths[samples.index(sample)] for sample in condition1_samples]
    condition2_paths = [paths[samples.index(sample)] for sample in condition2_samples]
    path_to_metilene = path_to_metilene + "/metilene_input.pl" if path_to_metilene != "" else "metilene_input.pl"
    command = f"perl {path_to_metilene} --in1 {','.join(condition1_paths)} --in2 {','.join(condition2_paths)} --out {out} --h1 {conditions[0]} --h2 {conditions[1]}"
    result = subprocess.run(command, shell=True)
    if result.returncode != 0:
        sys.exit('\033[91m' + f"Error creating input for Metilene with return code {result.returncode}" + '\033[0m')

def create_metilene_output(path_to_metilene, conditions, input, output, cores, maxdist=300, mincpgs=10, minMethDiff=0.1, additional_args=""):
    mode = 1
    threads = cores if cores != "" or cores != None or cores > 0 else 1
    path_to_metilene = path_to_metilene + "/metilene" if path_to_metilene != "" else "metilene"
    command = f"{path_to_metilene} -a {conditions[0]} -b {conditions[1]} -M {maxdist} -m {mincpgs} -d {minMethDiff} -t {threads} -f {mode} {additional_args} {input} > {output}"
    result = subprocess.run(command, shell=True)
    if result.returncode != 0:
        sys.exit('\033[91m' + f"Error running Metilene with return code {result.returncode}" + '\033[0m')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Prepare data for Metilene.')
    parser.add_argument('-m', '--matrix', required=True, help='Matrix file')

    parser.add_argument('-d', '--run_dir', required=True, help='Run directory')
    parser.add_argument('-mp', '--path_to_metilene', required=True, help='Path to Metilene')
    parser.add_argument('-s', '--samples', nargs='+', required=True, help='List of samples')
    parser.add_argument('-sp', '--samples_paths', nargs='+', required=True, help='List of all files')
    
    parser.add_argument('-o', '--output', required=True, help='Output file')
    parser.add_argument('-c', '--cores', required=False, help='Number of cores')
    parser.add_argument('-md', '--maxdist', required=False, help='Max distance')
    parser.add_argument('-mc', '--mincpgs', required=False, help='Min CpGs')
    parser.add_argument('-mdf', '--minMethDiff', required=False, help='Min Methhylation difference')
    parser.add_argument('-aa', '--additional_args', required=False, help='Additional arguments')
    args = parser.parse_args()

    print("Creating input for Metilene...")
    conditions, condition1_samples, condition2_samples = get_groups(args.matrix)
    create_input_for_metilene(args.run_dir, args.path_to_metilene, conditions, condition1_samples, condition2_samples, args.samples, args.samples_paths)
    print("Input for Metilene created successfully")
    print("Running Metilene...")
    create_metilene_output(args.path_to_metilene, conditions, f"{args.run_dir}/METHYLATION/DMR/input_for_metilene.txt", args.output, args.cores, args.maxdist, args.mincpgs, args.minMethDiff)
    print("Metilene finished successfully") 