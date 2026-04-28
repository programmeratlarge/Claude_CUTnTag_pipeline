/*
 * Merge multiple coordinate-sorted BAM files into a single BAM (samtools merge),
 * then index. Used for both treatment merge (by merge_group_id) and control
 * merge (by control's group_id).
 *
 * `kind` is just a tag for publishing ('treatment' or 'control').
 */

process MERGE_BAMS {
    tag "${group_id} (${kind})"
    label 'process_medium'
    publishDir "${params.outdir}/05_filtering/merged_${kind}", mode: params.publish_mode

    input:
    tuple val(group_id), val(meta), path(bams)
    val   kind

    output:
    tuple val(group_id), val(meta),
          path("${group_id}.merged.bam"),
          path("${group_id}.merged.bam.bai"), emit: bam

    script:
    if (bams instanceof List && bams.size() == 1) {
        // Single-sample group: just rename + index
        """
        cp ${bams[0]} ${group_id}.merged.bam
        samtools index -@ ${task.cpus} ${group_id}.merged.bam
        """
    } else {
        """
        samtools merge -@ ${task.cpus} -f ${group_id}.merged.bam ${bams}
        samtools index -@ ${task.cpus} ${group_id}.merged.bam
        """
    }
}
