import os
import re
import pandas as pd
import matplotlib.pyplot as plt
import argparse



def mbias_plot(file_path):
    with open(file_path, 'r') as file:
        text = file.read()
    for R in ['R1', 'R2']:
        print(f'Plotting {R} M-Bias')
        pattern = rf"CpG context \({R}\)\n=+\n(position\tcount methylated\tcount unmethylated\t% methylation\tcoverage\n(?:\d+\t[^\n]+\n?)*)"    
        match = re.search(pattern, text, re.DOTALL)
        if not match:
            return None
        data_block = match.group(1)
        lines = data_block.strip().split('\n')
        data = [line.split('\t') for line in lines]
        df = pd.DataFrame(data[1:], columns=data[0])    
        df = df.astype({'position': int, 'count methylated': int, 'count unmethylated': int, '% methylation': float, 'coverage': int})
        fig, ax1 = plt.subplots()
        fig.set_size_inches(20, 6)
        fig.suptitle(f'{os.path.basename(file_path).split(".")[0]} {R} M-Bias', fontsize=16)
        ax2 = ax1.twinx()
        ax1.plot(df['position'], df['% methylation'], 'g-')
        ax2.plot(df['position'], df['coverage'], 'b-')
        ax1.set_ylim(0, 100)
        ax1.set_xlim(left=0, right=len(df)-1)
        ax1.tick_params(axis='x', labelrotation=90)
        ax1.set_xticks(df['position'])
        ax1.set_xlabel('Position')
        ax1.set_ylabel('CpG Methylation', color='g')
        ax2.set_ylabel('Total Calls', color='b')
        ouput_file = file_path.replace('.txt', f'_{R}.png')
        print(f'Saving plot to {ouput_file}')
        plt.savefig(ouput_file)
        plt.close(fig) 

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plot M-Bias')
    parser.add_argument('-p', '--path', type=str, help='Path to M-Bias file')
    args = parser.parse_args()
    mbias_plot(args.path)