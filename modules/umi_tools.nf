/*
 * UMI-aware deduplication with umi_tools (used when --use_umi is set).
 * Assumes UMI has already been moved to the read name (e.g., by umi_tools extract).
 *
 * NOTE: For raw-FASTQ UMI extraction, run `umi_tools extract` upstream of cutadapt;
 * this module dedupes a coordinate-sorted, indexed BAM whose reads carry the UMI
 * appended to the read name with `_<UMI>`.
 */

process UMI_DEDUP {
    tag "${sample_id}"
    label 'process_medium'
    publishDir "${params.outdir}/05_filtering/dedup", mode: params.publish_mode

    input:
    tuple val(sample_id), val(meta), path(bam), path(bai)

    output:
    tuple val(sample_id), val(meta),
          path("${sample_id}.dedup.bam"),
          path("${sample_id}.dedup.bam.bai"), emit: bam
    tuple val(sample_id), val(meta),
          path("${sample_id}.umitools.log"), emit: log

    script:
    """
    umi_tools dedup \\
        --paired \\
        --stdin=${bam} \\
        --stdout=${sample_id}.dedup.bam \\
        --log=${sample_id}.umitools.log
    samtools index ${sample_id}.dedup.bam
    """
}
