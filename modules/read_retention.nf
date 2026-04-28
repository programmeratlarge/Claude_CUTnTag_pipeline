/*
 * Combine per-sample read counts from cutadapt, bowtie2, filter, and dedup
 * into one TSV + a MultiQC custom-content YAML.
 */

process READ_RETENTION {
    tag "${sample_id}"
    label 'process_low'
    publishDir "${params.outdir}/05_filtering", mode: params.publish_mode

    input:
    tuple val(sample_id),
          path(bowtie2_log),
          path(filter_counts),
          path(dedup_log),
          path(cutadapt_log)

    output:
    tuple val(sample_id), path("${sample_id}.retention.tsv"), emit: tsv
    path  "${sample_id}.retention_mqc.yaml", emit: mqc

    script:
    """
    read_retention_summary.py \\
        --sample      ${sample_id} \\
        --cutadapt    ${cutadapt_log} \\
        --bowtie2     ${bowtie2_log} \\
        --filter      ${filter_counts} \\
        --dedup       ${dedup_log} \\
        --out_tsv     ${sample_id}.retention.tsv \\
        --out_mqc     ${sample_id}.retention_mqc.yaml
    """
}
