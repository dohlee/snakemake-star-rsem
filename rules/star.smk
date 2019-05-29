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
        # genomeGeneration or 2_pass step.
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

rule star_2_pass_single:
    input:
        # Required input.
        reads = [
            DATA_DIR / '{sample}.fastq.gz',
        ],
        star_index = config['star_index_dir'],
    output:
        # There is no need to output sam or unsorted bam file!
        # So this wrapper includes '--outSAMtype BAM SortedByCoordinate' option by default.
        genome_alignment = RESULT_DIR / '01_star' / 'se' / '{sample}.sorted.bam',
        # Aligning to transcriptome is optional.
        # NOTE: You should set quantMode to 'TranscriptomeSAM' to get 
        # this output.
        transcriptome_alignment = RESULT_DIR / '01_star' / 'se' / '{sample}.transcriptome.bam',
    threads: config['threads']['star_2_pass_single']
    params:
        # Optional parameters. Read through the comments carefully.
        # Provide appropriate option for your data,
        # or comment out the option if it is not needed.
        extra = config['star_2_pass_single']['extra'],
        # Random number generator seed.
        # Default: 777
        run_rng_seed = config['star_2_pass_single']['run_rng_seed'],
        ### Quantification and Annotations
        # Types of quantification requested
        # False : none
        # TranscriptomeSAM: output SAM/BAM alignments to transcriptome into a separate file.
        # GeneCounts: count reads per gene
        # Default: False
        quant_mode = config['star_2_pass_single']['quant_mode'],
        # Mode of shared memory usage for the genome files.
        # Only used with --runMode 2_pass.
        # Default: NoSharedMemory
        genome_load = config['star_2_pass_single']['genome_load'],
        # Path to the files with genomic coordinates
        # (chr <tab> start <tab> end <tab> strand)
        # for the splice junction introns. Multiple files can be supplied
        # and will be concatenated.
        # Default: False
        sjdb_file_chr_start_end = config['star_2_pass_single']['sjdb_file_chr_start_end'],
        # It is recommended to give gene annotation file in at least one of
        # genomeGeneration or 2_pass step.
        # Also it is recommened to use GENCODE annotation file, either of GTF of GFF3 file.
        # Don't forget to specify sjdb_gtf_tag_exon_parent_transcript = 'Parent'
        # in case of GFF3 annotation file.
        # NOTE: Make sure the annotation file is unzipped!
        sjdb_gtf_file = config['star_2_pass_single']['sjdb_gtf_file'],
        # If genome is from UCSC, and annotation is from ENSEMBL,
        # use sjdb_gtf_chr_prefix = 'chr'
        sjdb_gtf_chr_prefix = config['star_2_pass_single']['sjdb_gtf_chr_prefix'],
        # Length of the donor/acceptor sequence on each side of the junctions,
        # ideally, (maxReadLength - 1).
        # In most cases, the default value of 100 will work well.
        sjdb_overhang = config['star_2_pass_single']['sjdb_overhang'],
        # If you use GFF3 annotation file,
        # use sjdb_gtf_tag_exon_parent_transcript = 'Parent'
        sjdb_gtf_tag_exon_parent_transcript = config['star_2_pass_single']['sjdb_gtf_tag_exon_parent_transcript'],
        ### Output filtering
        # Type of filtering.
        # Normal: standard filtering using only current alignment.
        # BySJout: keep only those reads that contain junctions that passed 
        # filtering into SJ.out.tab
        # Default: Normal
        out_filter_type = config['star_2_pass_single']['out_filter_type'],
        # The score range below the maximum score for multimapping alignments.
        # Default: 1
        out_filter_multimap_score_range = config['star_2_pass_single']['out_filter_multimap_score_range'],
        # Maximum number of loci the read is allowed to map to. Alignments (all of them)
        # will be output only if the read maps to no more loci than this value.
        # Otherwise no alignments will be output, and the read will be counted as
        # "mapped to too many loci" in the Log.final.out.
        # Default: 10
        out_filter_multimap_n_max = config['star_2_pass_single']['out_filter_multimap_n_max'],
        # Alignment will be output only if it has no more mismatches than this value.
        # Default: 10
        out_filter_mismatch_n_max = config['star_2_pass_single']['out_filter_mismatch_n_max'],
        # Alignment will be output only if its ratio of mismatches to *mapped* length
        # is less than or equal to this value.
        # Default: 0.3
        out_filter_mismatch_n_over_l_max = config['star_2_pass_single']['out_filter_mismatch_n_over_l_max'],
        # Alignment will be output only if its ratio of mismatches to *read* length
        # is less than or equal to this value.
        # Default: 1.0
        out_filter_mismatch_n_over_read_l_max = config['star_2_pass_single']['out_filter_mismatch_n_over_read_l_max'],
        # Alignment will be output only if its score is higher than or equal
        # to this value.
        # Default: 0
        out_filter_score_min = config['star_2_pass_single']['out_filter_score_min'],
        # Same as outFilterScoreMin, but normalized to read length.
        # (sum of mates' lengths for paired-end reads)
        # Default: 0.66
        out_filter_score_min_over_l_read = config['star_2_pass_single']['out_filter_score_min_over_l_read'],
        # Alignment will be output only if the number of matched bases is higher than or equal to this value.
        # Default: 0
        out_filter_match_n_min = config['star_2_pass_single']['out_filter_match_n_min'],
        # Same as outFilterMatchNmin, but normalized to the read length.
        # (sum of mates' lengths for pared-end reads)
        # Default: 0.66
        out_filter_match_n_min_over_l_read = config['star_2_pass_single']['out_filter_match_n_min_over_l_read'],
        # Filter alignment using their motifs
        # None: no filtering
        # RemoveNoncanonical: filter out alignments that contain non-canonical junctions
        # RemoveNoncanonicalUnannotated: filter out alignments that contain non-canonical
        #   unannotated junctions when using annotated splice junctions database.
        #   The annotated non-canonical junctions will be kept.
        # Default: None
        out_filter_intron_motifs = config['star_2_pass_single']['out_filter_intron_motifs'],
        # Filter alignments
        # RemoveInconsistentStrands: remove alignments that have junctions with
        # inconsistent strands. 
        # None: no filtering
        # Default: RemoveInconsistentStrands
        out_filter_intron_strands = config['star_2_pass_single']['out_filter_intron_strands'],
        ### Output filtering: Splice junctions
        # Which reads to consider for collapsed splice junctions output.
        # All: all reads, unique- and multi-mappers
        # Unique: uniquely mapping reads only
        # Default: All
        out_sj_filter_reads = config['star_2_pass_single']['out_sj_filter_reads'],
        # Minimum overhang length for splice junctions on both sides for:
        # (1) non-canonical motifs,
        # (2) GT/AG and CT/AC motif,
        # (3) GC/AG and CT/GC motif,
        # (4) AT/AC and /GT/AT motif.
        # -1 means no output for that motif.
        # This does not apply to annotated junctions
        # Default: '30 12 12 12'
        out_sj_filter_overhang_min = config['star_2_pass_single']['out_sj_filter_overhang_min'],
        # Minimum uniquely mapping read count per junction for:
        # (1) non-canonical motifs,
        # (2) GT/AG and CT/AC motif,
        # (3) GC/AG and CT/GC motif,
        # (4) AT/AC and GT/AT motif.
        # -1 means no output for that motif.
        # Junctions are output if one of outSJfilterCountUniqueMin OR outSJfilterCountTotalMin
        # conditions are satisfied.
        # This does not apply to annotated junctions.
        # Default: '3 1 1 1'
        out_sj_filter_count_unique_min = config['star_2_pass_single']['out_sj_filter_count_unique_min'],
        # Minimum total (multi-mapping+unique) read count per junction for:
        # (1) non-canonical motifs,
        # (2) GT/AG and CT/AC motif,
        # (3) GC/AG and CT/GC motif,
        # (4) AT/AC and GT/AT motif.
        # -1 means no output for that motif.
        # Junctions are output if one of outSJfilterCountUniqueMin OR outSJfilterCountTotalMin
        # conditions are satisfied.
        # This does not apply to annotated junctions.
        # Default: '3 1 1 1'
        out_sj_filter_count_total_min = config['star_2_pass_single']['out_sj_filter_count_total_min'],
        # Minimum allowed distance to other junctions' donor/acceptor.
        # This does not apply to annotated junctions.
        # Default: '10 0 5 10'
        out_sj_filter_dist_to_other_sj_min = config['star_2_pass_single']['out_sj_filter_dist_to_other_sj_min'],
        # Maximum gap allowed for junctions supported by 1,2,3,,,N reads.
        # i.e. by default junctions supported
        # by 1 read can have gaps <= 50000,
        # by 2 reads: <= 100000,
        # by 3 reads: <= 200000,
        # by >=4 reads any gap <= alignIntronMax
        # This does not apply to annotated junctions.
        # Default: '50000 100000 200000'
        out_sj_filter_intron_max_vs_read_n = config['star_2_pass_single']['out_sj_filter_intron_max_vs_read_n'],
        ### Alignments and seeding
        # Minimum intron size: genomic gap is considered intron if
        # its length>=alignIntronMin, otherwise it is considered deletion.
        # Default: 21
        align_intron_min = config['star_2_pass_single']['align_intron_min'],
        # maximum intron size: if 0, max intron size will be determined by
        # (2^winBinNbits)*winAnchorDistNbins
        # Default: 0
        align_intron_max = config['star_2_pass_single']['align_intron_max'],
        # Minimum overhang (i.e. block size) for spliced alignments.
        # Default: 5
        align_sj_overhang_min = config['star_2_pass_single']['align_sj_overhang_min'],
    log: 'logs/star_2_pass_single/{sample}.log'
    benchmark: 'benchmarks/star_2_pass_single/{sample}.benchmark'
    wrapper:
        'http://dohlee-bio.info:9193/star/2-pass'

rule star_2_pass_paired:
    input:
        # Required input.
        reads = [
            DATA_DIR / '{sample}.read1.fastq.gz',
            DATA_DIR / '{sample}.read2.fastq.gz',
        ],
        star_index = config['star_index_dir'],
    output:
        # There is no need to output sam or unsorted bam file!
        # So this wrapper includes '--outSAMtype BAM SortedByCoordinate' option by default.
        genome_alignment = RESULT_DIR / '01_star' / 'pe' / '{sample}.sorted.bam',
        # Aligning to transcriptome is optional.
        # NOTE: You should set quantMode to 'TranscriptomeSAM' to get 
        # this output.
        transcriptome_alignment = RESULT_DIR / '01_star' / 'pe' / '{sample}.transcriptome.bam',
    threads: config['threads']['star_2_pass_paired']
    params:
        # Optional parameters. Read through the comments carefully.
        # Provide appropriate option for your data,
        # or comment out the option if it is not needed.
        extra = config['star_2_pass_paired']['extra'],
        # Random number generator seed.
        # Default: 777
        run_rng_seed = config['star_2_pass_paired']['run_rng_seed'],
        ### Quantification and Annotations
        # Types of quantification requested
        # False : none
        # TranscriptomeSAM: output SAM/BAM alignments to transcriptome into a separate file.
        # GeneCounts: count reads per gene
        # Default: False
        quant_mode = config['star_2_pass_paired']['quant_mode'],
        # Mode of shared memory usage for the genome files.
        # Only used with --runMode 2_pass.
        # Default: NoSharedMemory
        genome_load = config['star_2_pass_paired']['genome_load'],
        # Path to the files with genomic coordinates
        # (chr <tab> start <tab> end <tab> strand)
        # for the splice junction introns. Multiple files can be supplied
        # and will be concatenated.
        # Default: False
        sjdb_file_chr_start_end = config['star_2_pass_paired']['sjdb_file_chr_start_end'],
        # It is recommended to give gene annotation file in at least one of
        # genomeGeneration or 2_pass step.
        # Also it is recommened to use GENCODE annotation file, either of GTF of GFF3 file.
        # Don't forget to specify sjdb_gtf_tag_exon_parent_transcript = 'Parent'
        # in case of GFF3 annotation file.
        # NOTE: Make sure the annotation file is unzipped!
        sjdb_gtf_file = config['star_2_pass_paired']['sjdb_gtf_file'],
        # If genome is from UCSC, and annotation is from ENSEMBL,
        # use sjdb_gtf_chr_prefix = 'chr'
        sjdb_gtf_chr_prefix = config['star_2_pass_paired']['sjdb_gtf_chr_prefix'],
        # Length of the donor/acceptor sequence on each side of the junctions,
        # ideally, (maxReadLength - 1).
        # In most cases, the default value of 100 will work well.
        sjdb_overhang = config['star_2_pass_paired']['sjdb_overhang'],
        # If you use GFF3 annotation file,
        # use sjdb_gtf_tag_exon_parent_transcript = 'Parent'
        sjdb_gtf_tag_exon_parent_transcript = config['star_2_pass_paired']['sjdb_gtf_tag_exon_parent_transcript'],
        ### Output filtering
        # Type of filtering.
        # Normal: standard filtering using only current alignment.
        # BySJout: keep only those reads that contain junctions that passed 
        # filtering into SJ.out.tab
        # Default: Normal
        out_filter_type = config['star_2_pass_paired']['out_filter_type'],
        # The score range below the maximum score for multimapping alignments.
        # Default: 1
        out_filter_multimap_score_range = config['star_2_pass_paired']['out_filter_multimap_score_range'],
        # Maximum number of loci the read is allowed to map to. Alignments (all of them)
        # will be output only if the read maps to no more loci than this value.
        # Otherwise no alignments will be output, and the read will be counted as
        # "mapped to too many loci" in the Log.final.out.
        # Default: 10
        out_filter_multimap_n_max = config['star_2_pass_paired']['out_filter_multimap_n_max'],
        # Alignment will be output only if it has no more mismatches than this value.
        # Default: 10
        out_filter_mismatch_n_max = config['star_2_pass_paired']['out_filter_mismatch_n_max'],
        # Alignment will be output only if its ratio of mismatches to *mapped* length
        # is less than or equal to this value.
        # Default: 0.3
        out_filter_mismatch_n_over_l_max = config['star_2_pass_paired']['out_filter_mismatch_n_over_l_max'],
        # Alignment will be output only if its ratio of mismatches to *read* length
        # is less than or equal to this value.
        # Default: 1.0
        out_filter_mismatch_n_over_read_l_max = config['star_2_pass_paired']['out_filter_mismatch_n_over_read_l_max'],
        # Alignment will be output only if its score is higher than or equal
        # to this value.
        # Default: 0
        out_filter_score_min = config['star_2_pass_paired']['out_filter_score_min'],
        # Same as outFilterScoreMin, but normalized to read length.
        # (sum of mates' lengths for paired-end reads)
        # Default: 0.66
        out_filter_score_min_over_l_read = config['star_2_pass_paired']['out_filter_score_min_over_l_read'],
        # Alignment will be output only if the number of matched bases is higher than or equal to this value.
        # Default: 0
        out_filter_match_n_min = config['star_2_pass_paired']['out_filter_match_n_min'],
        # Same as outFilterMatchNmin, but normalized to the read length.
        # (sum of mates' lengths for pared-end reads)
        # Default: 0.66
        out_filter_match_n_min_over_l_read = config['star_2_pass_paired']['out_filter_match_n_min_over_l_read'],
        # Filter alignment using their motifs
        # None: no filtering
        # RemoveNoncanonical: filter out alignments that contain non-canonical junctions
        # RemoveNoncanonicalUnannotated: filter out alignments that contain non-canonical
        #   unannotated junctions when using annotated splice junctions database.
        #   The annotated non-canonical junctions will be kept.
        # Default: None
        out_filter_intron_motifs = config['star_2_pass_paired']['out_filter_intron_motifs'],
        # Filter alignments
        # RemoveInconsistentStrands: remove alignments that have junctions with
        # inconsistent strands. 
        # None: no filtering
        # Default: RemoveInconsistentStrands
        out_filter_intron_strands = config['star_2_pass_paired']['out_filter_intron_strands'],
        ### Output filtering: Splice junctions
        # Which reads to consider for collapsed splice junctions output.
        # All: all reads, unique- and multi-mappers
        # Unique: uniquely mapping reads only
        # Default: All
        out_sj_filter_reads = config['star_2_pass_paired']['out_sj_filter_reads'],
        # Minimum overhang length for splice junctions on both sides for:
        # (1) non-canonical motifs,
        # (2) GT/AG and CT/AC motif,
        # (3) GC/AG and CT/GC motif,
        # (4) AT/AC and /GT/AT motif.
        # -1 means no output for that motif.
        # This does not apply to annotated junctions
        # Default: '30 12 12 12'
        out_sj_filter_overhang_min = config['star_2_pass_paired']['out_sj_filter_overhang_min'],
        # Minimum uniquely mapping read count per junction for:
        # (1) non-canonical motifs,
        # (2) GT/AG and CT/AC motif,
        # (3) GC/AG and CT/GC motif,
        # (4) AT/AC and GT/AT motif.
        # -1 means no output for that motif.
        # Junctions are output if one of outSJfilterCountUniqueMin OR outSJfilterCountTotalMin
        # conditions are satisfied.
        # This does not apply to annotated junctions.
        # Default: '3 1 1 1'
        out_sj_filter_count_unique_min = config['star_2_pass_paired']['out_sj_filter_count_unique_min'],
        # Minimum total (multi-mapping+unique) read count per junction for:
        # (1) non-canonical motifs,
        # (2) GT/AG and CT/AC motif,
        # (3) GC/AG and CT/GC motif,
        # (4) AT/AC and GT/AT motif.
        # -1 means no output for that motif.
        # Junctions are output if one of outSJfilterCountUniqueMin OR outSJfilterCountTotalMin
        # conditions are satisfied.
        # This does not apply to annotated junctions.
        # Default: '3 1 1 1'
        out_sj_filter_count_total_min = config['star_2_pass_paired']['out_sj_filter_count_total_min'],
        # Minimum allowed distance to other junctions' donor/acceptor.
        # This does not apply to annotated junctions.
        # Default: '10 0 5 10'
        out_sj_filter_dist_to_other_sj_min = config['star_2_pass_paired']['out_sj_filter_dist_to_other_sj_min'],
        # Maximum gap allowed for junctions supported by 1,2,3,,,N reads.
        # i.e. by default junctions supported
        # by 1 read can have gaps <= 50000,
        # by 2 reads: <= 100000,
        # by 3 reads: <= 200000,
        # by >=4 reads any gap <= alignIntronMax
        # This does not apply to annotated junctions.
        # Default: '50000 100000 200000'
        out_sj_filter_intron_max_vs_read_n = config['star_2_pass_paired']['out_sj_filter_intron_max_vs_read_n'],
        ### Alignments and seeding
        # Minimum intron size: genomic gap is considered intron if
        # its length>=alignIntronMin, otherwise it is considered deletion.
        # Default: 21
        align_intron_min = config['star_2_pass_paired']['align_intron_min'],
        # maximum intron size: if 0, max intron size will be determined by
        # (2^winBinNbits)*winAnchorDistNbins
        # Default: 0
        align_intron_max = config['star_2_pass_paired']['align_intron_max'],
        # Minimum overhang (i.e. block size) for spliced alignments.
        # Default: 5
        align_sj_overhang_min = config['star_2_pass_paired']['align_sj_overhang_min'],
    log: 'logs/star_2_pass_paired/{sample}.log'
    benchmark: 'benchmarks/star_2_pass_paired/{sample}.benchmark'
    wrapper:
        'http://dohlee-bio.info:9193/star/2-pass'

