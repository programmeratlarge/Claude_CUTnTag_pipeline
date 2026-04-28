/*
 * Run FastQC on a list of FASTQ files for a sample.
 * `stage` is a string such as 'raw' or 'trimmed' used for the output subdirectory.
 */

process FASTQC {
    tag "${sample_id} (${stage})"
    label 'process_low'
    publishDir "${params.outdir}/${stage == 'raw' ? '01_fastqc_raw' : '03_fastqc_trimmed'}",
        mode: params.publish_mode

    input:
    tuple val(sample_id), path(fastqs)
    val stage

    output:
    tuple val(sample_id), path('*_fastqc.zip'),  emit: zips
    tuple val(sample_id), path('*_fastqc.html'), emit: html

    script:
    """
    fastqc --quiet --threads ${task.cpus} ${fastqs}
    """
}
