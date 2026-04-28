/*
 * Samtools utility processes: sort+index, flagstat, idxstats, stats.
 * Each one tags reports with `stage` so they can be tracked through the pipeline.
 */

process SAMTOOLS_SORT_INDEX {
    tag "${sample_id} (${stage})"
    label 'process_medium'
    publishDir "${params.outdir}/04_alignment", mode: params.publish_mode

    input:
    tuple val(sample_id), val(meta), path(bam)
    val stage

    output:
    tuple val(sample_id), val(meta),
          path("${sample_id}.sorted.bam"),
          path("${sample_id}.sorted.bam.bai"), emit: bam

    script:
    """
    samtools sort -@ ${task.cpus} -o ${sample_id}.sorted.bam ${bam}
    samtools index -@ ${task.cpus} ${sample_id}.sorted.bam
    """
}

process SAMTOOLS_FLAGSTAT {
    tag "${sample_id} (${stage})"
    label 'process_low'
    publishDir "${params.outdir}/04_alignment/flagstat_${stage}", mode: params.publish_mode

    input:
    tuple val(sample_id), val(meta), path(bam), path(bai)
    val stage

    output:
    tuple val(sample_id), path("${sample_id}.${stage}.flagstat.txt"), emit: report

    script:
    """
    samtools flagstat ${bam} > ${sample_id}.${stage}.flagstat.txt
    """
}

process SAMTOOLS_IDXSTATS {
    tag "${sample_id} (${stage})"
    label 'process_low'
    publishDir "${params.outdir}/04_alignment/idxstats_${stage}", mode: params.publish_mode

    input:
    tuple val(sample_id), val(meta), path(bam), path(bai)
    val stage

    output:
    tuple val(sample_id), path("${sample_id}.${stage}.idxstats.txt"), emit: report

    script:
    """
    samtools idxstats ${bam} > ${sample_id}.${stage}.idxstats.txt
    """
}

process SAMTOOLS_STATS {
    tag "${sample_id} (${stage})"
    label 'process_low'
    publishDir "${params.outdir}/04_alignment/stats_${stage}", mode: params.publish_mode

    input:
    tuple val(sample_id), val(meta), path(bam), path(bai)
    val stage

    output:
    tuple val(sample_id), path("${sample_id}.${stage}.stats.txt"), emit: report

    script:
    """
    samtools stats ${bam} > ${sample_id}.${stage}.stats.txt
    """
}
