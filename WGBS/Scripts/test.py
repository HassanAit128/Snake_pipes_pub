import os
import subprocess
import re
import glob
import pandas as pd

def extract_bsmap_stats(file):
    stats = {}
    with open(file, "r") as f :
        data = f.read()
    stats['SampleID'] = re.search(r'\/([^\/]+)\.bam', data).group(1)
    stats['Total read pairs'] = int(re.search(r'total read pairs: (\d+)', data).group(1))
    stats['Aligned pairs'] = int(re.search(r'aligned pairs: (\d+)', data).group(1))
    stats['Aligned pairs %'] = float(re.search(r'aligned pairs: \d+ \((\d+\.\d+)%\)', data).group(1)) 
    stats['Unique pairs'] = int(re.search(r'unique pairs: (\d+)', data).group(1))
    stats['Non unique pairs'] = int(re.search(r'non-unique pairs: (\d+)', data).group(1))
    stats['Total time consumed (h)'] = int(re.search(r'total time consumed:  (\d+)', data).group(1))/3600
    return stats

def extract_bismark_stats(file):
    stats = {}
    with open(file, "r") as f :
        data = f.read()
    stats['SampleID'] = re.search(r'\/(SRR\d+)_', data).group(1)
    stats['Total read pairs'] = int(re.search(r'Sequence pairs analysed in total:\t(\d+)', data).group(1))
    stats['Unique pairs'] = int(re.search(r'Number of paired-end alignments with a unique best hit:\t(\d+)', data).group(1))
    stats['Unique pairs %'] = float(re.search(r'Mapping efficiency:\t(\d+\.\d+)%', data).group(1))
    stats['Non unique pairs'] = int(re.search(r'Sequence pairs did not map uniquely:\t(\d+)', data).group(1))
    stats['Total time consumed (h)'] = int(re.search(r'Bismark completed in (\d+)d (\d+)h (\d+)m (\d+)s', data).group(1)) * 24 + int(re.search(r'Bismark completed in (\d+)d (\d+)h (\d+)m (\d+)s', data).group(2)) + int(re.search(r'Bismark completed in (\d+)d (\d+)h (\d+)m (\d+)s', data).group(3))/60 + int(re.search(r'Bismark completed in (\d+)d (\d+)h (\d+)m (\d+)s', data).group(4))/3600
    return stats

def create_align_reports(run_dir, run_name, tool):
    stats = []
    if tool.lower() == "bsmap":
        align_reports = glob.glob(f"{run_dir}/{run_name}/LOGS/ALIGN_LOGS/*.log")
        for file in align_reports:
            stats.append(extract_bsmap_stats(file))
    elif tool.lower() == "bismark":
        align_reports = glob.glob(f"{run_dir}/{run_name}/*.txt")
        for file in align_reports:
            stats.append(extract_bismark_stats(file))
    return stats
        
# def samtools_count(bam_file, cores):
#     cmd = f"samtools view -@ {cores} -c {bam_file}"
#     count = subprocess.check_output(cmd, shell=True).decode().strip()
#     return count 

# def get_bam_stats(samples, run_dir, run_name, cores, downsampling):
#     stats = []
#     for sample in samples:
#         dedplicated_count = samtools_count(f"{run_dir}/{run_name}/BAMs/{sample}.deduplicated.bam", cores)
#         regular_count = samtools_count(f"{run_dir}/{run_name}/BAMs/{sample}.deduplicated.regular.bam", cores)
#         final_count = samtools_count(f"{run_dir}/{run_name}/BAMs/{sample}.deduplicated.regular.sorted.bam", cores)
#         if downsampling:
#             ds_count = samtools_count(f"{run_dir}/{run_name}/BAMs/{sample}.deduplicated.regular.ds.sorted.bam", cores)
#             stats.append({"SampleID": sample, "Deduplicated": dedplicated_count, "Regular Chr": regular_count, "Final_Count": final_count, "Final Count (million)": int(final_count)/1000000, "Downsampled": ds_count})
#         else:
#             stats.append({"SampleID": sample, "Deduplicated": dedplicated_count, "Regular Chr": regular_count, "Final_Count": final_count, "Final Count (million)": int(final_count)/1000000})
#     return stats

def get_merge_and_save_stats(run_dir, run_name, cores, tool, downsampling):
    align_stats = create_align_reports(run_dir, run_name, tool)
    samples = [sample["SampleID"] for sample in align_stats]
    # bam_stats = get_bam_stats(samples, run_dir, run_name, cores, tool, downsampling)
    align_stats_df = pd.DataFrame(align_stats)
    # bam_stats_df = pd.DataFrame(bam_stats)
    # merged_df = pd.merge(align_stats_df, bam_stats_df, on="SampleID")
    align_stats_df.to_excel(f"{run_dir}/{run_name}/Stats.xlsx", index=False)
    return align_stats_df

import pandas as pd 
import os
from datetime import datetime
import base64
from io import BytesIO
import argparse

def img_to_base64(img_path):
    with open(img_path, "rb") as img_file:
        return base64.b64encode(img_file.read()).decode('utf-8')
    
def styling():
    css = f"""
        body {{
            font-family: 'Helvetica Neue', Helvetica, Arial, sans-serif;
            margin: 0;
            padding: 0;
            background-color: #f9f9f9;
            color: #333;
            line-height: 1.6;
        }}

        header {{
            background-color: #2a3f54;
            color: white;
            padding: 20px;
            text-align: center;
            border-bottom: 3px solid #1a2732;
            box-shadow: 0 2px 8px rgba(0, 0, 0, 0.1);
        }}
        header h1 {{
            margin: 0;
            font-size: 2.5rem;
        }}
        header h5 {{
            margin: 5px 0;
            font-size: 1.2rem;
            font-weight: normal;
        }}

        main {{
            padding: 20px;
            max-width: 1200px;
            margin: auto;
            margin-top: 20px;
        }}

        .Main {{
            font-size: 1.2em; 
            padding: 20px;
            background-color: #fff;
            border: 1px solid #ccc; 
            border-radius: 10px;
            box-shadow: 0 2px 8px rgba(0, 0, 0, 0.1);
            margin-top: 20px;
        }}

        h2, h3 {{
            color: #2a3f54;
            margin-top: 2rem;
        }}

        h2 {{
            font-size: 2em;
        }}

        h3 {{
            font-size: 1.5em;
        }}

        p {{
            font-size: 1.1rem;
            margin: 10px 0;
        }}

        table {{
            border-collapse: collapse;
            width: 100%;
            margin: 20px 0;
            background-color: #fff;
            border-radius: 10px;
            box-shadow: 0 2px 8px rgba(0, 0, 0, 0.1);
        }}

        th, td {{
            border: 1px solid #ccc;
            padding: 12px 15px;
            text-align: left;
        }}

        th {{
            background-color: #f0f0f5;
            color: #333;
            font-weight: bold;
            border-top: none;
        }}

        tr:nth-child(even) {{
            background-color: #f9f9f9;
        }}

        tr:hover {{
            background-color: #f1f1f1;
        }}

        .scrollable-table {{
            overflow-x: auto;
            white-space: nowrap;
            position: relative;
        }}

        th:first-child, td:first-child {{
            position: sticky;
            left: 0;
            background-color: white;
            z-index: 1;
        }}

        .image-container {{
            display: flex;
            justify-content: space-around;
            align-items: center;
            margin-top: 2rem;
        }}

        .image-container div {{
            text-align: center;
            width: 48%;
        }}

        .image-container img {{
            max-width: 100%;
            height: auto;
            border: 1px solid #ddd;
            background-color: #fff;
        }}

        .image-title {{
            font-weight: bold;
            margin-bottom: 0.5rem;
            font-size: 1rem;
            color: #2a3f54;
        }}

        .condition-content {{
            display: none;
        }}

        .condition-content.active {{
            display: block;
        }}

        footer {{
            background-color: #2a3f54;
            color: white;
            text-align: center;
            padding: 20px;
            margin-top: 20px;
        }}
    """
    return css
    
def generate_html_code_main(dir, run, cores, tool, tool_methyl, downsampling, total_samples, multiqc_path, multiqc_trimmed_pathmed_link):
    css = styling()
    if not os.path.exists(multiqc_trimmed_pathmed_link):
        multiqc_trimmed_section = ""
    else:
        multiqc_trimmed_section = f"<p>Reports for trimmed reads can be found here: <span style='color: blue;'>{multiqc_trimmed_pathmed_link}</span></p>"
    stats = get_merge_and_save_stats(dir, run, cores, tool, downsampling)
    
    volcano_img_base64 = ""
    annotation_summary_1 = ""
    if os.path.exists(f"{dir}/{run}/gr/methylkit_VolcanoPlot.png"):
        volcano_img_base64 = img_to_base64(f"{dir}/{run}/gr/methylkit_VolcanoPlot.png")
    if os.path.exists(f"{dir}/{run}/gr/methylkit_summary.png"):
        annotation_summary_1 = img_to_base64(f"{dir}/{run}/gr/methylkit_summary.png")
    
    template = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <title>{run} Report</title>
        <style>
            {css}
        </style>
    </head>
    <body>
        <header>
            <h1>{run} Report</h1>
            <h5>Pipeline by: H. Aitou</h5>
            <h5>{datetime.now().strftime("%Y-%m-%d")}</h5>
        </header>
        <main>
            <section class="Main">
                <h2>Quality Control (QC):</h2>
                <p>Total number of samples: {total_samples}</p>
                <p>MultiQC results for non-trimmed reads can be found here: <span style="color: blue;">{multiqc_path}</span></p>
                {multiqc_trimmed_section}
            </section>
            <section class="Main">
                <h2>Analyis:</h2>
                <h3>Summary of mapping statistics:</h3>
                <p>Full summary of mapping statistics can be found here: <span style="color: blue;">{dir}/{run}.xlsx</span>.</p>
                <div class="scrollable-table">
                    {stats.to_html(index=False)}
                </div>
                <h3>Differential Analysis (<span style="color: blue;">{tool_methyl}</span>):</h3>
                <p>RData file can be found here: <span style="color: blue;">{dir}/{run}/METHYLATION/DMR/methylkit_analysis.RData</span>.</p>
                <p style="font-size: 15px;"> Note: Differential analysis results are also available in BED format for all regions significant or not. Another file in CSV format is also outputted for significant DMRs only. Both found in the same <a href="{dir}/{run}/METHYLATION/DMR">directory</a> as the RData file.</p>  
                <div class="image-container">
                    <div>
                        <div class="image-title">Number of DMRs:</div>
                        <img src="data:image/png;base64,{annotation_summary_1}" class="plot" alt="There was an error (Check the Logs)">
                    </div>
                    <div>
                        <div class="image-title">Volcano Plot of DMRs:</div>
                        <img src="data:image/jpeg;base64,{volcano_img_base64}" class="plot" alt="There was an error (Check the Logs)">
                    </div>
                </div>
            </section>
            <section class="Main">
                <h2>Annotation to TEs:</h2>
                
                
        </main>
        <footer>
            <p>Generated by the WGBS Analysis Pipeline - H. Aitou</p>
        </footer>
    </body>
    </html>
    """
    with open(f"{dir}/{run}/{run}_report.html", "w") as f:
        f.write(template)
    print(f"Report generated: {run}_report.html")



dir = r"C:\Users\hassa\Desktop\Internship\Report\wgbs_stats"
run = "bismark"
cores = 1
tool = "bismark"
tool_methyl = "methylkit"
downsampling = True
total_samples = 4
multiqc_path = r"C:\Users\hassa\Desktop\Internship\Report\wgbs_stats\multiqc_report.html"
multiqc_trimmed_path = r"C:\Users\hassa\Desktop\Internship\Report\wgbs_stats"
generate_html_code_main(dir, run, cores, tool, tool_methyl, downsampling, total_samples, multiqc_path, multiqc_trimmed_path)