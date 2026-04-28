/*
 * Paired-end alignment with Bowtie2 (CUT&Tag-tuned defaults).
 * Recommended Bowtie2 parameters for CUT&Tag:
 *   --end-to-end --very-sensitive --no-mixed --no-discordant
 *   -I 10 -X 700  (CUT&Tag fragments are typically 50-500 bp)
 */

process BOWTIE2_ALIGN {
    tag "${sample_id}"
    label 'process_high'
    publishDir "${params.outdir}/04_alignment", mode: params.publish_mode,
        saveAs: { fn -> fn.endsWith('.log') ? "logs/${fn}" : fn }

    input:
    tuple val(sample_id), val(meta), path(r1), path(r2)
    path  bt2_index_dir
    val   bt2_index_name

    output:
    tuple val(sample_id), val(meta), path("${sample_id}.bam"),         emit: bam
    tuple val(sample_id),            path("${sample_id}.bowtie2.log"), emit: log

    script:
    def bt2_args = params.bowtie2_args ?: '--end-to-end --very-sensitive --no-mixed --no-discordant -I 10 -X 700'
    """
    bowtie2 \\
        -x ${bt2_index_dir}/${bt2_index_name} \\
        -1 ${r1} -2 ${r2} \\
        -p ${task.cpus} \\
        ${bt2_args} \\
        2> ${sample_id}.bowtie2.log \\
      | samtools view -@ ${task.cpus} -b -o ${sample_id}.bam -
    """
}
