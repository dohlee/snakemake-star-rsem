from pathlib import Path

DATA_DIR = Path(config['data_dir'])
RESULT_DIR = Path(config['result_dir']) 

rule star_genome_generate:
    input:
        # Required input.
        # NOTE: Reference genome should be uncompressed.
        reference = config['reference']['fasta']
    output:
        index_directory = directory(config['star_index_dir'])
    threads: config['threads']['star_genome_generate']
    params:
        # Optional parameters. Read through the comments carefully.
        # Provide appropriate option for your data,
        # or comment out the option if it is not needed.
        extra = config['star_genome_generate']['extra'],
        # It is recommended to give gene annotation file in at least one of
        # genomeGeneration or alignReads step.
        # Also it is recommened to use GENCODE annotation file, either of GTF of GFF3 file.
        # Don't forget to specify sjdb_gtf_tag_exon_parent_transcript = 'Parent'
        # in case of GFF3 annotation file.
        sjdb_gtf_file = config['star_genome_generate']['sjdb_gtf_file'],
        # Length of the donor/acceptor sequence on each side of the junctions,
        # ideally = (maxReadLength - 1).
        # In most cases, the default value of 100 will work well.
        sjdb_overhang = config['star_genome_generate']['sjdb_overhang'],
        # If genome is from UCSC, and annotation is from ENSEMBL,
        # use sjdb_gtf_chr_prefix = 'chr'
        sjdb_gtf_chr_prefix = config['star_genome_generate']['sjdb_gtf_chr_prefix'],
        # If you use GFF3 annotation file,
        # use sjdb_gtf_tag_exon_parent_transcript = 'Parent'
        sjdb_gtf_tag_exon_parent_transcript = config['star_genome_generate']['sjdb_gtf_tag_exon_parent_transcript'],
    log: 'logs/star/genome_generate/%s.log' % config['reference']['name']
    wrapper:
        'http://dohlee-bio.info:9193/star/genome-generate'

rule star_2_pass_se:
    input:
        # Required input.
        reads = DATA_DIR / '{sample}.fastq.gz',
        star_index = directory(config['star_index_dir']),
    output:
        # There is no need to output sam or unsorted bam file!
        # So this wrapper includes '--outSAMtype BAM SortedByCoordinate' option by default.
        RESULT_DIR / '01_star' / 'se' / '{sample}.sorted.bam'
    threads: config['threads']['star_2_pass_se']
    params:
        # Optional parameters. Read through the comments carefully.
        # Provide appropriate option for your data,
        # or comment out the option if it is not needed.
        extra = config['star_2_pass_single']['extra'],
        # It is recommended to give gene annotation file in at least one of
        # genomeGeneration or alignReads step.
        # Also it is recommened to use GENCODE annotation file, either of GTF of GFF3 file.
        # Don't forget to specify sjdb_gtf_tag_exon_parent_transcript = 'Parent'
        # in case of GFF3 annotation file.
        # NOTE: Make sure the annotation file is unzipped!
        sjdb_gtf_file = config['star_2_pass_single']['sjdb_gtf_file'],
        # Length of the donor/acceptor sequence on each side of the junctions,
        # ideally = (maxReadLength - 1).
        # In most cases, the default value of 100 will work well.
        sjdb_overhang = config['star_2_pass_single']['sjdb_overhang'],
        # If genome is from UCSC, and annotation is from ENSEMBL,
        # use sjdb_gtf_chr_prefix = 'chr'
        sjdb_gtf_chr_prefix = config['star_2_pass_single']['sjdb_gtf_chr_prefix'],
        # If you use GFF3 annotation file,
        # use sjdb_gtf_tag_exon_parent_transcript = 'Parent'
        sjdb_gtf_tag_exon_parent_transcript = config['star_2_pass_single']['sjdb_gtf_tag_exon_parent_transcript'],
    log: 'logs/star/2-pass/{sample}.log'
    wrapper:
        'http://dohlee-bio.info:9193/star/2-pass'

rule star_2_pass_pe:
    input:
        # Required input.
        reads = [
            DATA_DIR / '{sample}.read1.fastq.gz',
            DATA_DIR / '{sample}.read2.fastq.gz',
        ],
        star_index = directory(config['star_index_dir']),
    output:
        # There is no need to output sam or unsorted bam file!
        # So this wrapper includes '--outSAMtype BAM SortedByCoordinate' option by default.
        RESULT_DIR / '01_star' / 'pe' / '{sample}.sorted.bam'
    threads: config['threads']['star_2_pass_pe']
    params:
        # Optional parameters. Read through the comments carefully.
        # Provide appropriate option for your data,
        # or comment out the option if it is not needed.
        extra = config['star_2_pass_single']['extra'],
        # It is recommended to give gene annotation file in at least one of
        # genomeGeneration or alignReads step.
        # Also it is recommened to use GENCODE annotation file, either of GTF of GFF3 file.
        # Don't forget to specify sjdb_gtf_tag_exon_parent_transcript = 'Parent'
        # in case of GFF3 annotation file.
        # NOTE: Make sure the annotation file is unzipped!
        sjdb_gtf_file = config['star_2_pass_single']['sjdb_gtf_file'],
        # Length of the donor/acceptor sequence on each side of the junctions,
        # ideally = (maxReadLength - 1).
        # In most cases, the default value of 100 will work well.
        sjdb_overhang = config['star_2_pass_single']['sjdb_overhang'],
        # If genome is from UCSC, and annotation is from ENSEMBL,
        # use sjdb_gtf_chr_prefix = 'chr'
        sjdb_gtf_chr_prefix = config['star_2_pass_single']['sjdb_gtf_chr_prefix'],
        # If you use GFF3 annotation file,
        # use sjdb_gtf_tag_exon_parent_transcript = 'Parent'
        sjdb_gtf_tag_exon_parent_transcript = config['star_2_pass_single']['sjdb_gtf_tag_exon_parent_transcript'],
    log: 'logs/star/2-pass/{sample}.log'
    wrapper:
        'http://dohlee-bio.info:9193/star/2-pass'

