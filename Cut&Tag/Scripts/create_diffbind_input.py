import csv
import re
import pandas as pd
import glob
import os


def open_user_design_matrix(DESIGN_MATRIX):
    design_matrix = {}
    with open(DESIGN_MATRIX, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            sample = row["SampleID"]
            design_matrix[sample] = row
    return design_matrix

def create_diffbind_matrix(dir, des_mat, paths, mapq, downsampled_samples, peak_called, peaks_path, peaks_path_caller):
    os.makedirs(f"{dir}/DiffBind", exist_ok=True)
    meta_data = open_user_design_matrix(des_mat)
    base_path = dir + "/BAM/"
    for q in mapq:
        files_q = [x for x in paths if f"mapq{q}" in x]
        files_q_rmdup = [x for x in files_q if "without_duplicates" in x]
        files_q_no_rmdup = [x for x in files_q if "with_duplicates" in x]            
        if files_q_rmdup != []:
            for f in files_q_rmdup:
                full_sample_name, file_extension = os.path.splitext(os.path.basename(f))
                sample = full_sample_name.split(".")[0]
                condition_sample_path = os.path.relpath(f, base_path).replace(".bam", "")
                meta_data[sample]['bamReads'] = f"{dir}/DS/{condition_sample_path}.downsampled.bam" if condition_sample_path in downsampled_samples else f"{dir}/BAM/{condition_sample_path}.bam"
                if peak_called == True and (peaks_path == None or peaks_path == ""):
                    meta_data[sample]['Peaks'] = f"{dir}/peak_calling/{condition_sample_path}_peaks.broadPeak"
                    meta_data[sample]['PeakCaller'] = "bed"
                elif peaks_path != None or peaks_path != "":
                    meta_data[sample]['Peaks'] = peaks_path
                    meta_data[sample]['PeakCaller'] = peaks_path_caller
            output = f"{dir}/DiffBind/diffbind_mapq{q}_without_duplicates.csv"
            with open(output, 'w', newline='') as f:
                writer = csv.DictWriter(f, fieldnames=meta_data[sample].keys())
                writer.writeheader()
                for f in files_q_rmdup:
                    full_sample_name, file_extension = os.path.splitext(os.path.basename(f))
                    sample = full_sample_name.split(".")[0]
                    writer.writerow(meta_data[sample])
        if files_q_no_rmdup != []:
            for f in files_q_no_rmdup:
                full_sample_name, file_extension = os.path.splitext(os.path.basename(f))
                sample = full_sample_name.split(".")[0]
                condition_sample_path = os.path.relpath(f, base_path).replace(".bam", "")
                meta_data[sample]['bamReads'] = f"{dir}/DS/{condition_sample_path}.downsampled.bam" if condition_sample_path in downsampled_samples else f"{dir}/BAM/{condition_sample_path}.bam"
                if peak_called == True and (peaks_path == None or peaks_path == ""):
                    meta_data[sample]['Peaks'] = f"{dir}/peak_calling/{condition_sample_path}_peaks.broadPeak"
                    meta_data[sample]['PeakCaller'] = "bed"
                elif peaks_path != None or peaks_path != "":
                    meta_data[sample]['Peaks'] = peaks_path
                    meta_data[sample]['PeakCaller'] = peaks_path_caller
            output = f"{dir}/DiffBind/diffbind_mapq{q}_with_duplicates.csv"
            with open(output, 'w', newline='') as f:
                writer = csv.DictWriter(f, fieldnames=meta_data[sample].keys())
                writer.writeheader()
                for f in files_q_no_rmdup:
                    full_sample_name, file_extension = os.path.splitext(os.path.basename(f))
                    sample = full_sample_name.split(".")[0]
                    writer.writerow(meta_data[sample])
                    
create_diffbind_matrix(dir = snakemake.params.run_dir, des_mat = snakemake.params.design_matrix, paths = snakemake.params.samples, mapq = snakemake.params.mapq, downsampled_samples = snakemake.params.downsampled_samples, peak_called = snakemake.params.peak_called, peaks_path = snakemake.params.peaks_path, peaks_path_caller = snakemake.params.peaks_path_caller)