rule rsem_prepare_reference:
    input:
        fasta = config['reference']['fasta'],
        # You cannot provide both GTF and GFF3 files.
        gtf = config['reference']['gtf'],
    output:
        config['rsem_reference']
    params:
        extra = config['rsem_prepare_reference']['extra'],
        # Give a comma-separated list of transcript categories,
        # e.g. "mRNA,rRNA". Only transcripts that match the pattern
        # will be extracted.
        # Default: "mRNA"
        gff3_rna_patterns = config['rsem_prepare_reference']['gff3_rna_patterns'],
        # Give a comma-separated list of trusted sources,
        # e.g. "ENSEMBL,HAVANA". Only transcripts coming from
        # these sources will be extracted. If this option is off,
        # all sources are accepted.
        # Default: False
        trusted_sources = config['rsem_prepare_reference']['trusted_sources'],
        # Use information from given file to map from transcript
        # (isoform) ids to gene ids. Each line of the file should be
        # of the form:
        #
        # gene_id transcript_id
        #
        # with the two fields separated by a tab character.
        # If you are using a GTF file for the "UCSC Genes" gene set from
        # the UCSC Genome Browser, then the "knownIsoforms.txt" file
        # (obtained from the "Downloads" section of the UCSC Genome Browser
        # site) is of this format.
        # If this option is off, then the mapping of isoforms to genes depends
        # on whether the '--gtf' option is specified. If '--gtf' is specified,
        # then RSEM uses the "gene_id" and "transcript_id" attributes in the GTF
        # file. Otherwise, RSEM assumes that each sequence in the reference
        # sequence files is a separate gene.
        # Default: False
        transcript_to_gene_map = config['rsem_prepare_reference']['transcript_to_gene_map'],
        # Use information from given file to provide gene_id and transcript_id
        # information for each allele-specific transcript. Each line of the file
        # should be of the form:
        #
        # gene_id transcript_id allele_id
        #
        # with the fields separated by a tab character.
        # This option is designed for quantifying allele-specific expression.
        # It is only valid if '--gtf' option is not specified. allele_id should be
        # the sequence names presented in the Multi-FASTA-formatted files.
        # Default: False
        allele_to_gene_map = config['rsem_prepare_reference']['allele_to_gene_map'],
        # Supress the output of logging information.
        # Default: False
        quiet = config['rsem_prepare_reference']['quiet'],
    threads: config['threads']['rsem_prepare_reference']
    log: 'logs/rsem_prepare_reference/_.log'
    wrapper:
        'http://dohlee-bio.info:9193/rsem/prepare-reference'

def rsem_calculate_expression_extra(wildcards):
    lib = SAMPLE2LIB[wildcards.sample]
    if lib.upper().startswith('SINGLE'):
        return config['rsem_calculate_expression']['extra']
    else:
        return ' '.join([config['rsem_calculate_expression']['extra'], '--paired-end'])

c = config['rsem_calculate_expression']
rule rsem_calculate_expression:
    input:
        reads = [
            RESULT_DIR / '02_star' / '{sample}.transcriptome.bam',
        ],
        # NOTE: REFERENCE_PREFIX.transcipts.fa should be generated via
        # rsem-prepare-reference
        reference = config['rsem_reference'],
    output:
        genes = RESULT_DIR / '03_rsem' / '{sample}.genes.results',
        isoforms = RESULT_DIR / '03_rsem' / '{sample}.isoforms.results'
    params:
        extra = rsem_calculate_expression_extra,
        # Input reads do not contain quality scores.
        # If this option is on, fasta input is assumed.
        no_qualities = c['no_qualities'],
        # This option defines the strandedness of the RNA-Seq reads.
        # It recognizes three values: 'none', 'forward', and 'reverse'.
        # 'none': Non-strand-specific protocols.
        # 'forward': All (upstream) reads are derived from the forward strand.
        # 'reverse': All (upstream) reads are derived from the reverse strand.
        # If forward/reverse is set, the '--norc'/'--nofw' Bowtie/Bowtie2 option
        # will also be enabled to avoid aligning reads to the opposite strand.
        # For Illumina truSeq Stranded protocols, please use 'reverse'.
        # Default: 'none'
        strandedness = c['strandedness'],
        # Input file contains alignments in SAM/BAM/CRAM format.
        # The exact file format will be determined automatically.
        # Default: False
        alignments = c['alignments'],
        # Use Bowtie2 instead of Bowtie to align reads.
        # Default: False
        bowtie2 = c['bowtie2'],
        # Use STAR to align reads. Alignment parameters are from ENCODE3's
        # STAR-RSEM pipeline.
        # Default: False
        star = c['star'],
        # If gene_name/transcript_name is available, append it to the end of
        # gene_id/transcript_id (separated by '_') in files
        # 'sample_name.isoforms.results' and 'sample_name.genes.results'
        # Default: False
        append_names = c['append_names'],
        # Set the seed for the random number generators used in calculating
        # posterior mean estimates and credibility intervals.
        # The seed must be a non-negative 32 bit integer.
        # Default: 0
        seed = c['seed'],
        # By default, RSEM uses Dirichlet(1) as the prior to calculate
        # posterior mean estimates and credibility intervals. However, much less
        # genes are expressed in single cell RNA-Seq data. Thus, if you want
        # to compute posterior mean estimates and/or credibility intervals
        # and you have single-cell RNA-Seq data, you are recommended to turn
        # this option. Then RSEM will use Dirichlet(0.1) as the prior which
        # encourage the sparsity of the expression levels.
        # Default: False
        single_cell_prior = c['single_cell_prior'],
        # Run RSEM's collapsed Gibbs sampler to calculate posterior mean estimates.
        # Default: False
        calc_pme = c['calc_pme'],
        # Calculate 95% credibility intervals and posterior mean esimates.
        # The credibility level can be changed by setting '--ci-credibility-level'.
        # Default: False
        calc_ci = c['calc_ci'],
        # Supress the output of logging information.
        # Default: False
        quiet = c['quiet'],
        # Sort BAM file aligned under transcript coordinate by read name.
        # Setting this option on will produce deterministic maximum likelihood
        # estimations from independent runs. Note that osorting will take long
        # time and lots of memory.
        # Default: False
        sort_bam_by_read_name = c['sort_bam_by_read_name'],
        # Do not output any BAM file.
        # Default: False
        no_bam_output = c['no_bam_output'],
        # When RSEM generates a BAM file, instead of outputting all alignments a
        # read has with their posterior probabilities, one alignment is sampled
        # according to the posterior probabilities. The sampling procedure
        # includes the alignment to the "noise" transcript, which does not
        # appear in the BAM file. Only the sampled alignment has a weight of 1.
        # All other alignments have weight 0. If the "noise" transcript is
        # sampled, all alignments appeared in the BAM file should have weight 0.
        # Default: False
        sampling_for_bam = c['sampling_for_bam'],
        # Generate a BAM file, 'sample_name.genome.bam', with alignments ampped
        # to genoic coordinates and annotated with their posterior probabilities.
        # In addition, RSEM will call samtools (included in RSEM package) to sort
        # and index the bam file.
        # 'sample_name.genome.sorted.bam' and 'sample_name.genome.sorted.bam.bai'
        # will be generated.
        # Default: False
        output_genome_bam = c['output_genome_bam'],
        # Sort RSEM generated transcript and genome BAM files by coordinates and
        # build associated indices.
        # Default: False
        sort_bam_by_coordinate = c['sort_bam_by_coordinate'],
        # Set the maximum memory per thread that can be used by 'samtools sort',
        # Maximum memory accepts sufficies 'K/M/G'. RSEM will pass it to the '-m'
        # option of 'samtools sort'. Note that the default used here is different
        # from the default used by samtools.
        # Default '1G'.
        sort_bam_memory_per_thread = c['sort_bam_memory_per_thread'],
    threads: config['threads']['rsem_calculate_expression']
    log: 'logs/rsem_calculate_expression_pe/{sample}.log'
    benchmark: 'benchmarks/rsem_calculate_expression_pe/{sample}.benchmark'
    wrapper:
        'http://dohlee-bio.info:9193/rsem/calculate-expression'

