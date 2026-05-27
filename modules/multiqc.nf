/*
 * Aggregate QC metrics with MultiQC.
 * `stage` is a string ('raw','trimmed','alignment','final') used for naming.
 * `custom_content` is a list of *_mqc.{tsv,yaml,json} files (may be empty).
 *
 * MultiQC's output naming convention when invoked with --filename FOO.html:
 *   FOO.html             (report)
 *   FOO_data/            (data tables — always written)
 *   FOO_plots/           (static plots — written when export_plots: true)
 * So with --filename multiqc_<stage>_report.html the directories are
 * multiqc_<stage>_report_data and multiqc_<stage>_report_plots.
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
    path "multiqc_${stage}_report.html",       emit: report
    path "multiqc_${stage}_report_data",       emit: data
    path "multiqc_${stage}_report_plots", optional: true, emit: plots

    script:
    """
    multiqc \\
        --config ${mqc_config} \\
        --filename multiqc_${stage}_report.html \\
        --force \\
        .
    """
}
