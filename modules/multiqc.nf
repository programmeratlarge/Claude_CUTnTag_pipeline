/*
 * Aggregate QC metrics with MultiQC.
 * `stage` is a string ('raw','trimmed','alignment','final') used for naming.
 * `custom_content` is a list of *_mqc.{tsv,yaml,json} files (may be empty).
 */

process MULTIQC {
    label 'process_low'
    publishDir "${params.outdir}/13_multiqc/${stage}", mode: params.publish_mode

    input:
    path qc_files
    path mqc_config
    val  stage
    path custom_content

    output:
    path "multiqc_${stage}_report.html", emit: report
    path "multiqc_${stage}_data",        emit: data

    script:
    def cc_dir_arg = custom_content ? '.' : '.'
    """
    multiqc \\
        --config ${mqc_config} \\
        --filename multiqc_${stage}_report.html \\
        --force \\
        .
    """
}
