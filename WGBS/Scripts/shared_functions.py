import os
import re
import csv
import sys
from pathlib import Path


def get_suffixes(sample, data_dir, trim=False):
    """
    Get the suffix and extension of the files for a given sample in the data directory
    
    Parameters
    ----------
    sample : str
        Sample name
    data_dir : str
        Path to the data directory
    trim : bool
        Whether function is called from the trim function
    
    Returns
    -------
    suffix : str
        Suffix of the file
    extension : str
        Extension of the file
    """
    pattern = re.compile(rf"{sample}_(R[12])(_?[^\.]*\.)?(.+)$")
    for filename in os.listdir(data_dir):
        match = pattern.search(filename)
        if match:
            suffix = match.group(2) if not trim else re.sub(r'\.$', '_', match.group(2))
            extension = match.group(3)
            return suffix, extension
    sys.exit(f"No matching files found for sample {sample} in directory {data_dir}")

def open_user_design_matrix(DESIGN_MATRIX):
    """
    Open the user design matrix and return a dictionary with the sample as key and the row as value
    
    Parameters
    ----------
    DESIGN_MATRIX : str
        Path to the user design matrix
        
    Returns
    -------
    design_matrix : dict
        Dictionary with the sample as key and the row as value
    """
    design_matrix = {}
    with open(DESIGN_MATRIX, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            sample = row["SampleID"]
            design_matrix[sample] = row
    return design_matrix

def check_if_design_matrix_is_valid(DESIGN_MATRIX, SAMPLES):
    """
    Check if the design matrix is valid by checking if all required columns are present, 
    if all samples in the design matrix are found in the data directory,
    if all samples in the data directory are found in the design matrix,
    and if all samples have values in the required columns.
    
    Parameters
    ----------
    DESIGN_MATRIX : str
        Path to the user design matrix
    SAMPLES : list
        List of samples found in the data directory
    
    Returns
    -------
    None
    """
    with open(DESIGN_MATRIX, 'r') as f:
        reader = list(csv.DictReader(f))
        samples_in_design_matrix = [samples['SampleID'] for samples in reader]
        # check if all required columns are present in the design matrix
        required_columns = ['SampleID', 'Condition', 'Replicate']
        for column in required_columns:
            if column not in reader[0].keys():
                sys.exit('\033[91m' + f"ERROR: Column '{column}' not found in design matrix." + '\033[0m')
        # check if all samples in sampleID are found in data directory
        for s in samples_in_design_matrix:
            if s not in SAMPLES:
                sys.exit('\033[91m' + f"ERROR: Sample '{s}' in design matrix but not found in data directory." + '\033[0m')
        # check if all samples in data directory are in design matrix
        for sample in SAMPLES:
            if sample not in samples_in_design_matrix:
                sys.exit('\033[91m' + f"ERROR: Sample '{sample}' in data directory but not found in design matrix." + '\033[0m')
        # check if all samples have values in the required columns
        for samples in reader:
            for column in required_columns:
                if samples[column] == '':
                    sys.exit('\033[91m' + f"ERROR: Sample '{samples['SampleID']}' has no value in column '{column}'" + '\033[0m')



def get_extension(file_list, sample):
    """
    Get the extension of the file for a given sample in the data directory
    
    Parameters
    ----------
    file_list : list
        List of files in the data directory
    sample : str
        Sample name
    
    Returns
    -------
    extension : str
        Extension of the file
    """
    path = [x for x in file_list if sample in x]
    extension = ''.join(Path(path[0]).suffixes)
    return extension