import os
import re
import pandas as pd
import matplotlib.pyplot as plt


def get_suffixes(sample, data_dir, trim=False):
    pattern = re.compile(rf"{sample}_(R[12])(_?[^\.]*\.)?(.+)$")
    for filename in os.listdir(data_dir):
        match = pattern.search(filename)
        if match:
            suffix = match.group(2) if not trim else re.sub(r'\.$', '_', match.group(2))
            extension = match.group(3)
            return suffix, extension
    raise FileNotFoundError(f"No matching files found for sample {sample} in directory {data_dir}")
