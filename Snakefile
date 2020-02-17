import pandas as pd
from pathlib import Path

configfile: 'config.yaml'

include: 'rules/star.smk'
include: 'rules/rsem.smk'
include: 'rules/trim-galore.smk'
include: 'rules/fastqc.smk'

ruleorder: trim_galore_pe > trim_galore_se

manifest = pd.read_csv(config['manifest'])
RESULT_DIR = Path(config['result_dir'])

SAMPLES = manifest.name.values
SAMPLE2LIB = {r.name:r.library_layout for r in manifest.to_records()}
SE_SAMPLES = manifest[manifest.library_layout == 'single'].name.values
PE_SAMPLES = manifest[manifest.library_layout == 'paired'].name.values

RAW_QC_SE = expand(str(DATA_DIR / '{sample}_fastqc.html'), sample=SE_SAMPLES)
RAW_QC_PE = expand(str(DATA_DIR / '{sample}.read1_fastqc.html'), sample=PE_SAMPLES)
TRIMMED_QC_SE = expand(str(RESULT_DIR / '01_trim_galore' / '{sample}.trimmed_fastqc.html'), sample=SE_SAMPLES)
TRIMMED_QC_PE = expand(str(RESULT_DIR / '01_trim_galore' / '{sample}.read1.trimmed_fastqc.html'), sample=PE_SAMPLES)
ALIGNED_BAM = expand(str(RESULT_DIR / '02_star' / '{sample}.sorted.bam'), sample=SAMPLES)
EXPRESSIONS = expand(str(RESULT_DIR / '03_rsem' / '{sample}.genes.results'), sample=SAMPLES)

RESULT_FILES = []
RESULT_FILES.append(RAW_QC_SE)
RESULT_FILES.append(RAW_QC_PE)
RESULT_FILES.append(TRIMMED_QC_SE)
RESULT_FILES.append(TRIMMED_QC_PE)
RESULT_FILES.append(ALIGNED_BAM)
RESULT_FILES.append(EXPRESSIONS)

rule all:
    input: RESULT_FILES
