/*
 * Detect paired-end FASTQ files in a directory and emit a TSV manifest.
 */

process DETECT_FASTQ_PAIRS {
    label 'process_low'
    publishDir "${params.outdir}/00_fastq_pairs", mode: params.publish_mode

    input:
    path  input_dir
    val   pattern

    output:
    path  'fastq_pairs.tsv',  emit: pairs_tsv
    path  'fastq_pairs.log',  emit: log

    script:
    """
    detect_fastq_pairs.py \\
        --input_dir ${input_dir} \\
        --pattern   '${pattern}' \\
        --out_tsv   fastq_pairs.tsv \\
        --log       fastq_pairs.log
    """
}
