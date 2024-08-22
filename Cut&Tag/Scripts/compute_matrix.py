import csv
import subprocess
import os
import argparse
import sys


def computeMatrix(bigwigs, annotated_regions, specific_region, design_matrix, run_dir, type, additional_options, plotprofile, plotprofile_additional_options, plotheatmap, plotheatmap_additional_options, custom_regions,cores, up, down, how_to_compute="all"):
    # seprate the bigwigs into groups for each condition in the design matrix
    with open(design_matrix, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        des_mat = list(reader)
    # Extract unique conditions
    conditions = list(set(row["Condition"] for row in des_mat))
    if "wt" in conditions[1].lower():
        conditions = conditions[::-1]
    # Extract sampleIDs for each condition
    condition1_samples = [row["SampleID"] for row in des_mat if row["Condition"] == conditions[0]]
    condition2_samples = [row["SampleID"] for row in des_mat if row["Condition"] == conditions[1]]
    # Separate the bigwigs into groups for each condition
    condition1_bigwigs = [bw for bw in bigwigs if any(sample in bw for sample in condition1_samples)]
    condition2_bigwigs = [bw for bw in bigwigs if any(sample in bw for sample in condition2_samples)]
    # compute the matrix for each annotated region
    print(f"Computing matrix for custom regions {custom_regions}")
    for an_region in annotated_regions:
        if custom_regions != "None":
            region = custom_regions
            print(f"Region: {region}")
            outFileSortedRegions = f"{run_dir}/Matrix/{conditions[0]}_{conditions[1]}_{custom_regions.split('/')[-1].split('.')[0]}.bed"
        elif specific_region != "None":
            # run a command to grep the specific region from the annotated regions
            region = f"{run_dir}/Matrix/{specific_region}_regions.bed"
            outFileSortedRegions = f"{run_dir}/Matrix/{specific_region}_{conditions[0]}_{conditions[1]}_{an_region.split('/')[-1].split('.')[0]}.bed"
            subprocess.run(f"grep {specific_region} {an_region} > {run_dir}/Matrix/{specific_region}_regions.bed", shell=True)
        else:
            region = an_region
            outFileSortedRegions = f"{run_dir}/Matrix/{conditions[0]}_{conditions[1]}_{an_region.split('/')[-1].split('.')[0]}.bed"
        additional_options = "" if additional_options == "None" else additional_options
        base_name = os.path.basename(an_region)        
        parts = base_name.split('_')
        name = '_'.join(parts[1:4])
        output_name = f"{run_dir}/Matrix/{name}_matrix.gz"
        if how_to_compute.lower() == "all":
            print(f"Computing matrix for {name}")
            command = f"computeMatrix {type} {additional_options} -a {up} -b {down} -S {' '.join(condition1_bigwigs) + ' ' + ' '.join(condition2_bigwigs)} -R {region} -o {output_name} -p {cores} --outFileSortedRegions {outFileSortedRegions} --samplesLabel {' '.join(condition1_samples) + ' ' + ' '.join(condition2_samples)}"
            print(f"Running command: {command}")
            subprocess.run(command, shell=True)
            # check if no error occured
            if not os.path.exists(output_name):
                sys.exit(f"Error: Matrix computation failed for {name}")
            print(f"Matrix created for {name}")
            if plotheatmap == "True":
                print(f"Plotting heatmap for {name}")
                plotHeatmap(output_name, run_dir, plotheatmap_additional_options)
            if plotprofile == "True":
                print(f"Plotting profile for {name}")
                plotProfile(output_name, run_dir, plotprofile_additional_options)

def plotHeatmap(matrix_file, run_dir, additional_options):
    additional_options = "" if additional_options == "None" else additional_options
    base_name = os.path.basename(matrix_file)        
    parts = base_name.split('_')
    name = '_'.join(parts[0:3])
    output_name = f"{run_dir}/Matrix/{name}_heatmap.png"
    command = f"plotHeatmap -m {matrix_file} -o {output_name} {additional_options}"
    print(f"Running command: {command}")
    subprocess.run(command, shell=True)

def plotProfile(matrix_file, run_dir, additional_options):
    additional_options = "" if additional_options == "None" else additional_options
    base_name = os.path.basename(matrix_file)        
    parts = base_name.split('_')
    name = '_'.join(parts[0:3])
    output_name = f"{run_dir}/Matrix/{name}_profile.png"
    command = f"plotProfile -m {matrix_file} -o {output_name} {additional_options}"
    print(f"Running command: {command}")
    subprocess.run(command, shell=True)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Compute matrix for regions of interest.')
    parser.add_argument('-b', '--bigwigs', nargs='+', required=True, help='List of bigwigs')
    parser.add_argument('-a', '--annotated_regions', nargs='+', required=True, help='List of annotated regions')
    parser.add_argument('-sr', '--specific_region', required=True, help='Specific regions to compute matrix for')
    parser.add_argument('-d', '--design_matrix', required=True, help='Design matrix file')
    parser.add_argument('-r', '--run_dir', required=True, help='Run directory')
    parser.add_argument('-t', '--type', required=True, help='Type of matrix to compute')
    parser.add_argument('-ao', '--additional_options', required=False, help='Additional options for computeMatrix')
    parser.add_argument('-htc', '--how_to_compute', required=False, help='How to compute the matrix')
    parser.add_argument('-pp', '--plotprofile', required=False, help='Plot profile')
    parser.add_argument('-ph', '--plotheatmap', required=False, help='Plot heatmap')
    parser.add_argument('-ppa', '--plotprofile_additional_options', required=False, help='Additional options for plotProfile')
    parser.add_argument('-pha', '--plotheatmap_additional_options', required=False, help='Additional options for plotHeatmap')
    parser.add_argument('-cr', '--custom_regions', required=False, help='Custom regions to compute matrix for')
    parser.add_argument('-p', '--cores', required=False, help='Number of cores to use')
    parser.add_argument('-up', '--upstream', required=False, help='Upstream distance from the reference')
    parser.add_argument('-down', '--downstream', required=False, help='Downstream distance from the reference')
    args = parser.parse_args()
    
    computeMatrix(args.bigwigs, args.annotated_regions, args.specific_region, args.design_matrix, args.run_dir, args.type, args.additional_options, args.plotprofile, args.plotprofile_additional_options, args.plotheatmap, args.plotheatmap_additional_options, args.custom_regions, args.cores, args.upstream, args.downstream, "all")