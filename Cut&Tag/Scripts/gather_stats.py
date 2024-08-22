import csv
import re
import pandas as pd
import glob
import os
from shared_functions import open_user_design_matrix

def transform_picard_stats_to_dict(data):
    """
    Extract the Picard stats from the txt file and transform it into a dictionary
    
    Parameters
    ----------
    data : str
        Content of the txt file
    
    Returns
    -------
    stats : dict
        Dictionary with the Picard stats
    """
    stats = {}
    picard_stats = re.search(r'### Picard ###\n(.*?)(?=###|$)', data, re.DOTALL).group(1)
    picard_stats = picard_stats.split('\n')
    picard_stats = [x.split('\t') for x in picard_stats if x != '']
    stats["PCR_Duplicates__Reads"] = int(picard_stats[1][6])
    stats["PCR_duplicates__percent"] = float(picard_stats[1][8])
    return stats

def transform_processing_stats_to_dict(data, rmdup):
    """
    Extract the processing stats from the txt file and transform it into a dictionary
    
    Parameters
    ----------
    data : str
        Content of the txt file
    rmdup : bool
        If the duplicates were removed
    
    Returns
    -------
    stats : dict
        Dictionary with the processing stats
    """
    def process_section(variable_processing_stats, stats):
        """
        Inner function to process the stats
        
        Parameters
        ----------
        variable_processing_stats : str
            Section of the stats
        stats : dict
            Dictionary with the stats
        
        Returns
        -------
        None
        """
        for line in variable_processing_stats.split('\n'):
            if line != '':
                key, value = line.split(': ')
                column = re.search(r'\((.*?)\)', key)
                if column is not None:
                    column = column.group(1)
                    stats[column] = int(value)

    stats = {}
    processing_1_stats = re.search(r'### Proccessing 1 ###\n(.*?)(?=###|$)', data, re.DOTALL).group(1)
    process_section(processing_1_stats, stats)
    
    if rmdup:
        processing_2_stats = re.search(r'### Proccessing 2 ###\n(.*?)(?=###|$)', data, re.DOTALL).group(1)
        process_section(processing_2_stats, stats)
    
    stats["_Final__Million"] = stats["_Final__Reads"] / 1000000
    return stats
    
def transform_bowtie2_stats_to_dict(data):
    """
    Extract the Bowtie2 stats from the txt file and transform it into a dictionary
    
    Parameters
    ----------
    data : str
        Content of the txt file
    
    Returns
    -------
    stats : dict
        Dictionary with the Bowtie2 stats
    """
    stats = {}
    stats['Raw_reads'] = int(re.search(r'(\d+) reads; of these:', data).group(1))*2
    stats['Uniq_mapped__Reads'] = int(re.search(r'(\d+) \((\d+\.\d+)%\) aligned concordantly exactly 1 time', data).group(1))*2
    stats['Uniq_mapped__percent (%)'] = float(re.search(r'(\d+\.\d+)%\) aligned concordantly exactly 1 time', data).group(1))
    stats['MutiMapped__Reads'] = int(re.search(r'(\d+) \((\d+\.\d+)%\) aligned concordantly >1 times', data).group(1))*2
    stats['MutiMapped__percent'] = float(re.search(r'(\d+\.\d+)%\) aligned concordantly >1 times', data).group(1))
    stats['Unmapped__Reads'] = int(re.search(r'(\d+) \((\d+\.\d+)%\) aligned concordantly 0 times', data).group(1))*2
    stats['Unmapped__Percent'] = float(re.search(r'(\d+\.\d+)%\) aligned concordantly 0 times', data).group(1))
    stats['Uniq_and_MutiMapped__Reads'] = int(stats['Uniq_mapped__Reads']) + int(stats['MutiMapped__Reads'])
    stats['Uniq_and_MutiMapped__percent'] = float(re.search(r'(\d+\.\d+)% overall alignment rate', data).group(1))
    return stats

def extract_stats_from_txt_file(txt_file, rmdup):
    """
    Get the stats from the txt file
    
    Parameters
    ----------
    txt_file : str
        Path to the txt file
    rmdup : bool
        If the duplicates were removed
    
    Returns
    -------
    stats : dict
        Dictionary with the stats
    """
    with open(txt_file, 'r') as f:
        data = f.read()
    stats = {}
    stats['SampleID'] = re.search(r'Stats for : ([^\n]+)', data).group(1)
    bowtie2_stats = transform_bowtie2_stats_to_dict(data)
    stats.update(bowtie2_stats)
    if rmdup:
        picard_stats = transform_picard_stats_to_dict(data)
        stats.update(picard_stats)
    samtools_stats = transform_processing_stats_to_dict(data, rmdup)
    stats.update(samtools_stats)
    return stats
    
def create_stats_excel_for_different_mapq(design_matrix, mapq, dir, samples):
    """
    Create an Excel file with the stats for each sample and mapq value
    
    Parameters
    ----------
    design_matrix : str
        Path to the user design matrix
    mapq : list
        List of mapq values
    dir : str
        Path to the run directory
    samples : list
        List of samples
    
    Returns
    -------
    None
    """
    def save_to_excel(files, output, meta_data, rmdup):
        """
        Inner function to save the stats to an Excel file
        
        Parameters
        ----------
        files : list
            List of files
        output : str
            Path to the output Excel file
        meta_data : dict
            Dictionary with the meta data
        rmdup : bool
            If the duplicates were removed
        
        Returns
        -------
        None
        """
        stats = []
        outputs = []
        for file in files:
            stats.append(extract_stats_from_txt_file(file, rmdup))
        for stat in stats:
            out = {}
            sample = stat['SampleID']
            out.update(meta_data[sample])
            out.update(stat)
            outputs.append(out)
        df = pd.DataFrame(outputs)
        df.to_excel(output, index=False)
        
    for sample in samples:
        os.remove(f"{dir}/stats/{sample}_stats.txt")
    files = glob.glob(f"{dir}/stats/*_stats.txt")
    meta_data = open_user_design_matrix(design_matrix)
    for q in mapq:
        files_q = [x for x in files if f"mapq{q}" in x]
        files_q_rmdup = [x for x in files_q if "without_duplicates" in x]
        files_q_no_rmdup = [x for x in files_q if "with_duplicates" in x]
        output_rmdup = f"{dir}/stats/Summary/mapq{q}_without_duplicates.xlsx" if files_q_rmdup != [] else None
        output_no_rmdup = f"{dir}/stats/Summary/mapq{q}_with_duplicates.xlsx" if files_q_no_rmdup != [] else None
        if files_q_rmdup != []:
            save_to_excel(files_q_rmdup, output_rmdup, meta_data, rmdup=True)
        if output_no_rmdup is not None:
            save_to_excel(files_q_no_rmdup, output_no_rmdup, meta_data, rmdup=False)
            
create_stats_excel_for_different_mapq(design_matrix = snakemake.params.design_matrix, mapq = snakemake.params.mapq, dir = snakemake.params.run_dir, samples = snakemake.params.samples)
