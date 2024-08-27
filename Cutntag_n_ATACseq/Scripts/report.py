import pandas as pd 
import os
from datetime import datetime
import base64
from io import BytesIO
import argparse


def img_to_base64(img_path):
    with open(img_path, "rb") as img_file:
        return base64.b64encode(img_file.read()).decode('utf-8')
    
def generate_html_code_main(run, total_samples, multiqc_path, multiqc_trimmed_pathmed_link, conditions):
    if not os.path.exists(multiqc_trimmed_pathmed_link):
        multiqc_trimmed_section = ""
    else:
        multiqc_trimmed_section = f"<p>MultiQC results for trimmed reads can be found here : <a href='{multiqc_trimmed_pathmed_link}'>MultiQC</a></p>"
    dropdown_options = "\n".join([f"<option value='{condition}'>{condition}</option>" for condition in conditions])    
    template = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <title>{run} Report</title>
        <style>
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
                
            .dropdown {{
                font-size: 1.2em; 
                padding: 20px;
                background-color: #fff;
                border: 1px solid #ccc; 
                border-radius: 10px;
                box-shadow: 0 2px 8px rgba(0, 0, 0, 0.1);
                margin-top: 20px;
            }}

            .dropdown h2 {{
                font-size: 1.5em;
                margin-bottom: 10px; 
                color: #2a3f54;
            }}

            .dropdown label {{
                display: block;
                margin-bottom: 10px; 
                color: #2a3f54;
            }}

            .dropdown select {{
                width: 100%; 
                padding: 12px 15px; 
                font-size: 1em; 
                border: 1px solid #ccc; 
                border-radius: 5px; 
                background-color: #f9f9f9;
                box-shadow: inset 0 1px 3px rgba(0, 0, 0, 0.1);
                transition: border-color 0.3s ease-in-out;
            }}

            .dropdown select:focus {{
                border-color: #2a3f54;
                outline: none;
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
        </style>
        <script>
            function showCondition(conditionId) {{
                var conditions = document.getElementsByClassName('condition-content');
                for (var i = 0; i < conditions.length; i++) {{
                    conditions[i].classList.remove('active');
                }}
                document.getElementById(conditionId).classList.add('active');
            }}
            window.onload = function() {{
                document.getElementById('condition-dropdown').value = '{conditions[0]}';
                showCondition('{conditions[0]}');
            }};
        </script>
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
                <p>Total number of samples: {total_samples}</p>
                <p>MultiQC results for non-trimmed reads can be found here: <a href="{multiqc_path}">MultiQC</a></p>
                {multiqc_trimmed_section}
            </section>
            <section class="Main">
                <h2>Analyis:</h2>
                <section class="dropdown">
                    <select id="condition-dropdown" onchange="showCondition(this.value)">
                        {dropdown_options}
                    </select>
                </section>
    """
    return template

def generate_html_code_per_condition(dir, summary, condition, Mode, df, limit):
    pca_img_base64 = ""
    volcano_img_base64 = ""
    annotation_summary_1 = ""
    annotation_summary_2 = ""
    deeptools_profile = ""
    deeptools_heatmap = ""
    header = []
    if os.path.exists(f"{dir}/DiffBind/diffbind_{condition}_PCA.png"):
        pca_img_base64 = img_to_base64(f"{dir}/DiffBind/diffbind_{condition}_PCA.png")
    if os.path.exists(f"{dir}/DiffBind/diffbind_{condition}_Volcano.png"):
        volcano_img_base64 = img_to_base64(f"{dir}/DiffBind/diffbind_{condition}_Volcano.png")
    if os.path.exists(f"{dir}/DiffBind/diffbind_{condition}_status_summary.png"):
        annotation_summary_1 = img_to_base64(f"{dir}/DiffBind/diffbind_{condition}_status_summary.png")
    if os.path.exists(f"{dir}/DiffBind/diffbind_{condition}_diffbind_annotated_TEs_summary_2.png"):
        annotation_summary_2 = img_to_base64(f"{dir}/DiffBind/diffbind_{condition}_diffbind_annotated_TEs_summary_2.png")
    if os.path.exists(f"{dir}/DiffBind/diffbind_{condition}_diffbind_header.txt"):
        with open(f"{dir}/DiffBind/diffbind_{condition}_diffbind_header.txt", "r") as f:
            lines = f.readlines()
            header = [line.strip() for line in lines]
    if os.path.exists(f"{dir}/Matrix/{condition}_profile.png"):
        deeptools_profile = img_to_base64(f"{dir}/Matrix/{condition}_profile.png")
    if os.path.exists(f"{dir}/Matrix/{condition}_heatmap.png"):
        deeptools_heatmap = img_to_base64(f"{dir}/Matrix/{condition}_heatmap.png")
    template = f"""
        <div id="{condition}" class="condition-content">
            <main>
                <h3>Summary of mapping statistics:</h3>
                <p>Full summary of mapping statistics for <span style="font-weight: bold; color: blue;">{condition}</span> can be found <a href="{dir}/{condition}.xlsx">here</a>.</p>
                <p style="font-size: 10px;">*Final Reads are color-coded based on the ENCODE minimum requirement of {limit} reads for the {Mode} analysis.</p>
                <div class="scrollable-table">
                    {summary}
                </div>
                <h3>Differential Analysis (Diffbind)</h3>
                <p>Diffbind results for <span style="font-weight: bold; color: blue;">{condition}</span> can be found <a href="{dir}/diffbind_{condition}_diffbind.bed">here.</a></p>
                <p style="font-size: 10px;">(Note: The Diffbind results are in BED format, the header should look like this: {header}, also found <a href="{dir}/DiffBind/diffbind_{condition}_diffbind_header.txt">here</a>.)</p>
                <div class="image-container">
                    <div>
                        <div class="image-title">PCA Plot:</div>
                        <img src="data:image/png;base64,{pca_img_base64}" class="plot" alt="There was an error (Check the Logs)">
                    </div>
                    <div>
                        <div class="image-title">Volcano Plot:</div>
                        <img src="data:image/jpeg;base64,{volcano_img_base64}" class="plot" alt="There was an error (Check the Logs)">
                    </div>
                </div>
                <h3>Annotation of Peaks</h3>
                <p>Annotation of peaks for <span style="font-weight: bold; color: blue;">{condition}</span> can be found <a href="{dir}/diffbind_{condition}_diffbind_annotated.bed">here</a>.</p>
                <div class="image-container">
                    <div>
                        <div class="image-title">Summary of peak annotation 1:</div>
                        <img src="data:image/png;base64,{annotation_summary_1}" class="plot" alt="There was an error (Check the Logs)">
                    </div>
                    <div>
                        <div class="image-title">Summary of peak annotation 2:</div>
                        <img src="data:image/jpeg;base64,{annotation_summary_2}" class="plot" alt="There was an error (Check the Logs)">
                    </div>
                </div>
                <h3>Profile Plot and Heatmap</h3>
                <p>Profile plot and heatmap for <span style="font-weight: bold; color: blue;">{condition}</span> can be found <a href="{dir}/Matrix/{condition}_profile.png">here</a> and <a href="{dir}/Matrix/{condition}_heatmap.png">here</a>.</p>
                <div class="image-container">
                    <div>
                        <div class="image-title">Profile Plot:</div>
                        <img src="data:image/png;base64,{deeptools_profile}" class="plot" alt="There was an error (Check the Logs)">
                    </div>
                    <div>
                        <div class="image-title">Heatmap:</div>
                        <img src="data:image/jpeg;base64,{deeptools_heatmap}" class="plot" alt="There was an error (Check the Logs)">
                    </div>
                </div>
            </main>
        </div>
    """
    return template

def create_html_file(dir, summaries, run, multiqc_path, multiqc_trimmed_path, total_samples, Mode):
    limit = 50000000 if Mode == "ATAC-Seq" else 2000000 if Mode == "CUT&Tag" else 0
    final_html_page = ""
    columns = ["SampleID", "Condition", "Raw_reads", "Uniq_and_MutiMapped__Reads", "Uniq_and_MutiMapped__percent", "PCR_duplicates__percent", "_Final__Reads", "_Final__Million"]
    conditions = [summary.split("/")[-1].split(".")[0] for summary in summaries]
    final_html_page += generate_html_code_main(run, total_samples, multiqc_path, multiqc_trimmed_path, conditions)
    for summary in summaries:
        condition = summary.split("/")[-1].split(".")[0]
        df = pd.read_excel(summary)
        df = df[columns]
        styled_df = apply_coloring(df, '_Final__Reads', limit)
        html_table = styled_df.to_html(index=False)
        html_code = generate_html_code_per_condition(dir, html_table, condition, Mode, df, limit)
        final_html_page += html_code
    final_html_page += f"""
        </section>
        </main>
        <footer>
            <p>Generated by the {Mode} Analysis Pipeline - H. Aitou</p>
        </footer>
    </body>
    </html>
    """
    with open(f"{dir}/{run}_report.html", "w") as f:
        f.write(final_html_page)
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

def color_scale(value, min_value, limit):
    if value >= limit:
        return 'background-color: green'
    else:
        norm_value = (value - min_value) / (limit - min_value)
        red = int(255 * (1 - norm_value))
        green = int(255 * norm_value)
        return f'background-color: rgb({red},{green},0)'

def apply_coloring(df, column, limit):
    min_value = df[column].min()
    return df.style.applymap(lambda x: color_scale(x, min_value, limit) if column in df.columns else '', subset=[column])

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate a report for the pipeline')
    parser.add_argument('--run', type=str, help='Name of the run')
    parser.add_argument('--summaries', nargs='+', type=str, help='List of summaries')
    parser.add_argument('--dir', type=str, help='Directory of the run')
    parser.add_argument('--total_samples', type=int, help='Total number of samples')
    parser.add_argument('--multiqc_path', type=str, help='Path to the MultiQC report')
    parser.add_argument('--multiqc_trimmed_path', type=str, help='Path to the MultiQC report for trimmed reads')
    parser.add_argument('--Mode', type=str, help='Type of analysis')
    args = parser.parse_args()
    create_html_file(args.dir, args.summaries, args.run, args.multiqc_path, args.multiqc_trimmed_path, args.total_samples, args.Mode)
    