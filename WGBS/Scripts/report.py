import pandas as pd 
import os
from datetime import datetime
import base64
from io import BytesIO
import argparse
from final_report import *


def img_to_base64(img_path):
    """
    Convert an image to base64 for easy embedding in HTML.
    
    Parameters
    ----------
    img_path : str
        Path to the image.
    
    Returns
    -------
    str
        Base64 encoded image.
    """
    with open(img_path, "rb") as img_file:
        return base64.b64encode(img_file.read()).decode('utf-8')
    
def styling():
    """Template to return the CSS styling for the HTML report."""
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
    
def generate_html_code_main(design_matrix, dir, run, cores, tool, tool_methyl, downsampling, total_samples, multiqc_path, multiqc_trimmed_pathmed_link, genome):
    """
    Main function to generate the HTML report.
    
    Parameters
    ----------
    design_matrix : str
        Path to the design matrix.
    dir : str
        Directory of the run.
    run : str
        Name of the run.
    cores : int
        Number of cores to use for samtools count.
    tool : str
        Alignment tool used.
    tool_methyl : str
        Methylation tool used.
    downsampling : bool
        If downsampling was used.
    total_samples : int
        Total number of samples.
    multiqc_path : str
        Path to the MultiQC report.
    multiqc_trimmed_pathmed_link : str
        Path to the MultiQC report for trimmed reads.
    genome : str
        Genome used for the analysis.
    
    Returns
    -------
    None
    """
    print("--- Starting to generate HTML report ---")
    css = styling()
    print("--- Figuring out paths and HTML content ---")
    if not os.path.exists(multiqc_trimmed_pathmed_link):
        multiqc_trimmed_section = ""
    else:
        multiqc_trimmed_section = f"<p>Reports for trimmed reads can be found here: <span style='color: blue;'>{multiqc_trimmed_pathmed_link}</span></p>"
    stats = get_merge_and_save_stats(design_matrix, dir, run, cores, tool, downsampling)
    volcano_img_base64 = ""
    annotation_summary_1 = ""
    annotation_summary_2 = ""
    annotation_summary_3 = ""
    if os.path.exists(f"{dir}/{run}/METHYLATION/DMR/VolcanoPlot.png"):
        volcano_img_base64 = img_to_base64(f"{dir}/{run}/METHYLATION/DMR/VolcanoPlot.png")
    if os.path.exists(f"{dir}/{run}/METHYLATION/DMR/DMR_summary.png"):
        annotation_summary_1 = img_to_base64(f"{dir}/{run}/METHYLATION/DMR/DMR_summary.png")
    if os.path.exists(f"{dir}/{run}/METHYLATION/DMR/TE_summary.png"):
        annotation_summary_2 = img_to_base64(f"{dir}/{run}/METHYLATION/DMR/TE_summary.png")
    if os.path.exists(f"{dir}/{run}/METHYLATION/DMR/TE_summary_2.png"):
        annotation_summary_3 = img_to_base64(f"{dir}/{run}/METHYLATION/DMR/TE_summary_2.png")
    print("--- OK ---")
    print("--- Generating HTML ---")
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
                    <h2>Quality Control (QC)</h2>
                    <p>Total number of samples: <span style="color: blue;">{total_samples}</span>.</p>
                    <p>MultiQC results for non-trimmed reads can be found here: <span style="color: blue;">"{multiqc_path}"</span>.</p>
                    {multiqc_trimmed_section}
                </section>
                <section class="Main">
                    <h2>Analyis:</h2>
                    <h3>Summary of mapping statistics:</h3>
                    <p>Full summary of mapping statistics can be found here: <span style="color: blue;">"{dir}/{run}.xlsx"</span>.</p>
                    <div class="scrollable-table">
                        {stats.to_html(index=False)}
                    </div>
                    <h3>Differential Analysis (<span style="color: blue;">{tool_methyl}</span>):</h3>
                    <p>RData file can be found here: <span style="color: blue;">"{dir}/{run}/METHYLATION/DMR/methylkit_analysis.RData".</span></p>
                    <p style="font-size: 15px;"> Note: Differential analysis results are also available in BED format for all regions significant or not. Another file in CSV format is also outputted for significant DMRs only. Both found in the same directory as the RData file.</p>  
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
                    <h3>Annotation of DMRs:</h3>
                    <p>Annotation of DMRs with set filters can be found here: <span style="color: blue;">"{dir}/{run}/METHYLATION/DMR/annotated_DMRs.txt"</span>.</p>
                    <div class="image-container">
                        <div>
                            <div class="image-title">Distribution of each family of TE among the TEs that are found to be differentially enriched (right).:</div>
                            <img src="data:image/png;base64,{annotation_summary_2}" class="plot" alt="There was an error (Check the Logs)">
                        </div>
                        <div>
                            <div class="image-title">Disturbution of DMRs within TEs:</div>
                            <img src="data:image/png;base64,{annotation_summary_3}" class="plot" alt="There was an error (Check the Logs)">
                        </div>
                        
                            
                </section>
            </main>
            <footer>
                <p>Generated by the WGBS Analysis Pipeline - H. Aitou</p>
            </footer>
        </body>
        </html>
        """
    print("--- OK ---")
    print("--- Writing HTML content to file ---")
    with open(f"{dir}/{run}/FINAL_REPORT/{run}_report.html", "w") as f:
        f.write(template)
    print("--- OK ---")
    print(f"Report generated: {run}_report.html")
    
# def generate_graph(df, limit):
#     plt.figure(figsize=(10, 6))
#     plt.bar(df.index, df['_Final__Reads'])
#     plt.xlabel('Samples')
#     plt.ylabel('Final Reads')
#     plt.xticks(ticks=df.index, labels=df['SampleID'], rotation='vertical')
#     if limit != 0:
#         plt.axhline(y=limit, color='r', linestyle='--', label='ENCODE minimum requirement')
#         plt.legend(loc='upper right')
#     plt.tight_layout()
#     buffer = BytesIO()
#     plt.savefig(buffer, format='png')
#     plt.close()
#     buffer.seek(0)
#     img_base64 = base64.b64encode(buffer.read()).decode('utf-8')
#     return img_base64


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate a report for the pipeline')
    parser.add_argument('--design_matrix', type=str, help='Path to the design matrix')
    parser.add_argument('--dir', type=str, help='Directory of the run')
    parser.add_argument('--run', type=str, help='Name of the run')
    parser.add_argument('--cores', type=int, help='Number of cores to use for samtools count')
    parser.add_argument('--tool', type=str, help='Alignment tool used')
    parser.add_argument('--tool_methyl', type=str, help='Methylation tool used')
    parser.add_argument('--downsampling', type=bool, help='If downsampling was used')
    parser.add_argument('--total_samples', type=int, help='Total number of samples')
    parser.add_argument('--multiqc_path', type=str, help='Path to the MultiQC report')
    parser.add_argument('--multiqc_trimmed_path', type=str, help='Path to the MultiQC report for trimmed reads')
    parser.add_argument('--genome', type=str, help='Genome used for the analysis')
    args = parser.parse_args()
    generate_html_code_main(args.design_matrix, args.dir, args.run, args.cores, args.tool, args.tool_methyl, args.downsampling, args.total_samples, args.multiqc_path, args.multiqc_trimmed_path, args.genome)
