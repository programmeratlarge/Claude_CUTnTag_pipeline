/*
 * Picard MarkDuplicates — removes PCR/optical duplicates.
 * For CUT&Tag this is the standard dedup approach unless UMIs are used.
 */

process MARK_DUPLICATES {
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
          path("${sample_id}.markdup.metrics.txt"), emit: metrics

    script:
    def avail_mem = (task.memory ? task.memory.toGiga() : 8)
    """
    picard -Xmx${avail_mem}g MarkDuplicates \\
        I=${bam} \\
        O=${sample_id}.dedup.bam \\
        M=${sample_id}.markdup.metrics.txt \\
        REMOVE_DUPLICATES=true \\
        ASSUME_SORT_ORDER=coordinate \\
        VALIDATION_STRINGENCY=LENIENT
    samtools index ${sample_id}.dedup.bam
    """
}
