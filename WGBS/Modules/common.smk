import os
import glob
from Scripts.shared_functions import get_suffixes
from snakemake.utils import min_version
import sys

min_version("8.0.1")

ansi = {'underline': '\033[4m', 'bold': '\033[1m', 'red': '\033[91m', 'end': '\033[0m'}
############################################# CONFIGURATION ################################################
RUN = config['name']
WORKING_DIR = config['working_dir'] 
if not os.path.exists(WORKING_DIR):
    sys.exit('\033[91m' + 'ERROR: You need to specify a valid working directory' + '\033[0m')
PROJECT_DIR = WORKING_DIR + '/' + RUN
RAW_DATA_DIR = config['raw_data_dir'] 
if not os.path.exists(RAW_DATA_DIR):
    sys.exit('\033[91m' + 'ERROR: You need to specify a valid path to the raw data directory' + '\033[0m')
DESIGN_MATRIX = config['design_matrix']
if not os.path.exists(DESIGN_MATRIX):
    sys.exit('\033[91m' + 'ERROR: You need to specify a valid path to the design matrix file' + '\033[0m')
TRIM = config['trim_reads']['run_trimming'] 
TRIM_ARGS = config['trim_reads']['trim_galore_args'] 
TRIM_CORES = config['trim_reads']['trim_galore_cores'] 
if TRIM_CORES is None or TRIM_CORES == 0:
    TRIM_CORES = 1

GENOME = config['alignment']['genome']
if GENOME not in ['hg38', 'mm10']:
    sys.exit('\033[91m' + 'ERROR:You need to specify a valid genome (i.e. hg38 or mm10)' + '\033[0m')
GENOME_PATH = config['alignment']['path_to_genome_fasta']
if not os.path.exists(GENOME_PATH):
    sys.exit('\033[91m' + 'ERROR: You need to specify a valid path to the genome .fasta/.fa file' + '\033[0m')

ALIGN = True
idx_output = []
if ALIGN:
    ALIGN_TOOL = config['alignment']['alignment_tool'].lower() 
    if ALIGN_TOOL == 'bismark':
        IDX_GENOME = config['bismark']['run_indexing'] 
        if not IDX_GENOME:
            GENOME_IDX_DIR = config['bismark']['path_to_bismark_idx']
            if not os.path.exists(GENOME_IDX_DIR):
                sys.exit('\033[91m' + 'ERROR: You need to specify a valid path to the genome index' + '\033[0m')
        else:
            sys.exit('\033[91m' + 'ERROR: You need to specify a path to the genome fasta file or set run_indexing to True' + '\033[0m')
        ALIGN_CORES = config['bismark']['bismark_cores']
        if ALIGN_CORES is None or ALIGN_CORES == 0:
            ALIGN_CORES = 1
        ALIGN_ARGS = config['bismark']['bismark_args']
    elif ALIGN_TOOL == 'bsmap':
        ALIGN_CORES = config['bsmap']['bsmap_cores']
        if ALIGN_CORES is None or ALIGN_CORES == 0:
            ALIGN_CORES = 1
        ALIGN_ARGS = config['bsmap']['bsmap_args']
        ALIGN_MODE = config['bsmap']['bsmap_mode'].lower()
        if ALIGN_MODE not in ['unique', 'multi']:
            sys.exit('\033[91m' + 'ERROR: Alignment mode must be either unique or multi' + '\033[0m')
    else:
        sys.exit('\033[91m' + 'ERROR: Alignment tool must be either bismark or bsmap' + '\033[0m')

DEDUP = True
CALL_METHYLATION_TOOL = config['call_methylation']['methylation_caller'].lower() 

if DEDUP:
    DEDUP_TOOL = config['deduplication']['deduplication_tool'].lower()
    if DEDUP_TOOL == 'picard':
        DEDUP_ARGS = config['deduplication']['picard_args'] 
    elif DEDUP_TOOL == 'bismark':
        if ALIGN_TOOL == 'bsmap':
            sys.exit('\033[91m' + 'ERROR: Bismark deduplication is not compatible with BSMap alignment' + '\033[0m')
        if CALL_METHYLATION_TOOL == 'methyldackel':
            sys.exit('\033[91m' + 'ERROR: Bismark deduplication seems to cause errors with Methyldackel methylation calling, please use Picard instead.' + '\033[0m')
        DEDUP_ARGS = config['deduplication']['bismark_deduplication_args'] 
    else:
        sys.exit('\033[91m' + 'ERROR: Deduplication tool must be either picard or bismark' + '\033[0m')
    
SAMTOOLS_CORES = config['samtools']['samtools_cores']
if SAMTOOLS_CORES is None or SAMTOOLS_CORES == 0:
    SAMTOOLS_CORES = 1
KEEP_REGULAR = True
DOWNSAMPLE = config['downsampling']['run_downsampling'] 

if CALL_METHYLATION_TOOL == 'bismark':
    if ALIGN_TOOL == 'bsmap':
        sys.exit('\033[91m' + 'ERROR: Bismark methylation calling is not compatible with BSMap alignment' + '\033[0m')
    CALL_METH_CORES = config['bismark_call']['bismark_methylation_cores']
    CALL_METH_ARGS = config['bismark_call']['bismark_methylation_args']
elif CALL_METHYLATION_TOOL == 'methyldackel':
    CALL_METH_CORES = config['Methyldackel']['methyldackel_cores']
    CALL_MBIAS_ARGS = config['Methyldackel']['methyldackel_mbias_args']
    CALL_METH_ARGS = config['Methyldackel']['methyldackel_args']
    AUTO_IGNORE = config['Methyldackel']['auto_ignore']
    ALSO_METHYLKIT = config['Methyldackel']['also_methylkit']
else:
    sys.exit('\033[91m' + "ERROR: Methylation caller must be either bismark or methyldackel" + "\033[0m")

CALL_METHYLATION_BIAS = config['call_methylation']['run_methylation_bias_check']
if CALL_METHYLATION_TOOL == 'methyldackel' and AUTO_IGNORE == True:
    CALL_METHYLATION_BIAS = True
STOP_ALL = config['call_methylation']['stop_after_m_bias_check'] 
CALL_METHYLATION = config['call_methylation']['run_methylation_calling'] 

DMR_analysis = config['DMR_analysis']['run_DMR_analysis']
DMR_analysis_tool = config['DMR_analysis']['DMR_analysis_tool'].lower()
if DMR_analysis :
    if not CALL_METHYLATION:
        sys.exit('\033[91m' + 'ERROR: You need to run methylation calling before DMR analysis' + '\033[0m')
    MIN_COV = config['DMR_analysis']["minimum_coverage"]
    if DMR_analysis_tool == 'metilene':
        METILENE_PATH = config['metilene']['path_to_metilene']
        METILENE_CORES = config['metilene']['metilene_cores']
        METILENE_ARGS = config['metilene']['metilene_args']
    elif DMR_analysis_tool == 'methylkit':
        pass
    else:
        sys.exit('\033[91m' + 'ERROR: DMR analysis tool must be either metilene or methylkit' + '\033[0m')

GENERATE_REPORT = True
CLEANUP = config['cleanup'] 

FILE_LIST = []
patterns = ['*.fq', '*.fq.gz', '*.fastq', '*.fastq.gz']
for pattern in patterns:
    FILE_LIST.extend(glob.glob(os.path.join(RAW_DATA_DIR, '*_R1*' + pattern)))
pattern = re.compile(r'(.+)_R1.*\.(fq|fq\.gz|fastq|fastq\.gz)$')
SAMPLES = [pattern.match(os.path.basename(x)).group(1) for x in FILE_LIST if pattern.match(os.path.basename(x))]
if FILE_LIST == [] or SAMPLES == []:
    sys.exit('\033[91m' + 'ERROR: NO FASTQ files found in the project directory. Check the path to the FASTQ files' + '\033[0m')
############################################################################################################

