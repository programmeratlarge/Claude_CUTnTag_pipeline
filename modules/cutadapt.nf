/*
 * Paired-end 3' adapter trimming with cutadapt.
 * Default adapters are Nextera/Tn5 (CTGTCTCTTATACACATCT) — appropriate for
 * CUT&Tag because Tn5 introduces these sequences during tagmentation.
 */

process CUTADAPT {
    tag "${sample_id}"
    label 'process_medium'
    publishDir "${params.outdir}/02_trimmed", mode: params.publish_mode,
        saveAs: { fn -> fn.endsWith('.log') ? "logs/${fn}" : fn }

    input:
    tuple val(sample_id), val(meta), path(r1), path(r2)

    output:
    tuple val(sample_id), val(meta), path("${sample_id}_R1.trim.fastq.gz"),
                                     path("${sample_id}_R2.trim.fastq.gz"), emit: trimmed
    tuple val(sample_id), val(meta), path("${sample_id}.cutadapt.log"),     emit: log

    script:
    def adapter_fwd     = params.adapter_fwd      ?: 'CTGTCTCTTATACACATCT'
    def adapter_rev     = params.adapter_rev      ?: 'CTGTCTCTTATACACATCT'
    def min_read_length = params.min_read_length  ?: 20
    """
    cutadapt \\
        -j ${task.cpus} \\
        -a ${adapter_fwd} \\
        -A ${adapter_rev} \\
        --minimum-length ${min_read_length} \\
        -q 20 \\
        -o ${sample_id}_R1.trim.fastq.gz \\
        -p ${sample_id}_R2.trim.fastq.gz \\
        ${r1} ${r2} \\
        > ${sample_id}.cutadapt.log
    """
}
