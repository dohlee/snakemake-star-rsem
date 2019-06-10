configfile: 'config.yaml'
REFERENCE_FASTA = config['reference']['fasta']
REFERENCE_GTF = config['reference']['gtf']
RNA_SEQ_SINGLE_TMP = 'data/test_single_full.fastq.gz'
RNA_SEQ_SINGLE = 'data/test_single.fastq.gz'
RNA_SEQ_PAIRED_TMP = ['data/test_paired_full.read1.fastq.gz', 'data/test_paired_full.read2.fastq.gz']
RNA_SEQ_PAIRED = ['data/test_paired.read1.fastq.gz', 'data/test_paired.read2.fastq.gz']

ALL = []
ALL.append(REFERENCE_FASTA)
ALL.append(REFERENCE_GTF)
ALL.append(RNA_SEQ_SINGLE)
ALL.append(RNA_SEQ_PAIRED)

rule all:
    input: ALL

rule reference:
    output: config['reference']['fasta']
    wrapper: 'http://dohlee-bio.info:9193/test/reference'

rule reference_transcriptome:
    output: config['rsem_reference']
    wrapper: 'http://dohlee-bio.info:9193/test/reference/transcriptome'

rule annotation:
    output: config['reference']['gtf']
    wrapper: 'http://dohlee-bio.info:9193/test/annotation'

rule rna_seq_single:
    output: temp(RNA_SEQ_SINGLE_TMP)
    wrapper: 'http://dohlee-bio.info:9193/test/rna-seq/se'

rule rna_seq_paired:
    output: temp(RNA_SEQ_PAIRED_TMP)
    wrapper: 'http://dohlee-bio.info:9193/test/rna-seq/pe'

rule subsample_fastq_single:
    input:
        # Required input.
        reads = RNA_SEQ_SINGLE_TMP
    output:
        RNA_SEQ_SINGLE
    threads: 1  # No more than 1 threads will be used.
    params:
        # Required parameters.
        k = 5000  # Number of sampled reads.
    wrapper:
        'http://dohlee-bio.info:9193/subsample-fastq'

rule subsample_fastq_paired:
    input:
        # Required input.
        reads = RNA_SEQ_PAIRED_TMP
    output:
        RNA_SEQ_PAIRED
    threads: 1  # No more than 1 threads will be used.
    params:
        # Required parameters.
        k = 5000  # Number of sampled reads.
    wrapper:
        'http://dohlee-bio.info:9193/subsample-fastq'

