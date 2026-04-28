/*
 * FRiP (Fraction of Reads in Peaks). Two flavors:
 *   FRIP_SAMPLE: per-sample BAM vs (own | merged-group) peaks
 *   FRIP_GROUP : merged-group BAM vs merged-group peaks
 * Each emits a TSV row + a MultiQC custom-content YAML so they show up in the final report.
 */

process FRIP_SAMPLE {
    tag "${sample_id} (${peak_type})"
    label 'process_low'
    publishDir "${params.outdir}/11_frip/sample_${peak_type}", mode: params.publish_mode

    input:
    tuple val(sample_id), val(meta), path(bam), path(bai), path(peaks), val(peak_type)

    output:
    path "${sample_id}.${peak_type}.frip.tsv", emit: tsv
    path "${sample_id}.${peak_type}.frip_mqc.yaml", emit: mqc

    script:
    """
    compute_frip.py \\
        --sample ${sample_id} \\
        --bam    ${bam} \\
        --peaks  ${peaks} \\
        --label  ${peak_type} \\
        --out_tsv ${sample_id}.${peak_type}.frip.tsv \\
        --out_mqc ${sample_id}.${peak_type}.frip_mqc.yaml
    """
}

process FRIP_GROUP {
    tag "${merge_group_id}"
    label 'process_low'
    publishDir "${params.outdir}/11_frip/group", mode: params.publish_mode

    input:
    tuple val(merge_group_id), val(meta), path(bam), path(bai), path(peaks)

    output:
    path "${merge_group_id}.group.frip.tsv", emit: tsv
    path "${merge_group_id}.group.frip_mqc.yaml", emit: mqc

    script:
    """
    compute_frip.py \\
        --sample ${merge_group_id} \\
        --bam    ${bam} \\
        --peaks  ${peaks} \\
        --label  group \\
        --out_tsv ${merge_group_id}.group.frip.tsv \\
        --out_mqc ${merge_group_id}.group.frip_mqc.yaml
    """
}
