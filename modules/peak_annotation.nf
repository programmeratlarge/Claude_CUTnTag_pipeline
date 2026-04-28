/*
 * Annotate peaks with ChIPseeker (R) — produces:
 *   - annotation TSV
 *   - per-set pie chart of categories
 *   - feature distribution barplot
 *   - MultiQC custom-content TSV
 */

process ANNOTATE_PEAKS {
    tag "${peakset_id}"
    label 'process_medium'
    publishDir "${params.outdir}/10_peak_annotation", mode: params.publish_mode

    input:
    tuple val(peakset_id), path(peaks)
    path  annotation_gtf
    val   tss_window

    output:
    tuple val(peakset_id), path("${peakset_id}_annotation.tsv"), emit: annotation
    path  "${peakset_id}_anno_pie.pdf",       optional: true
    path  "${peakset_id}_anno_features.pdf",  optional: true
    path  "${peakset_id}_anno_categories_mqc.tsv", emit: mqc

    script:
    """
    annotate_peaks.R \\
        --peaks       ${peaks} \\
        --gtf         ${annotation_gtf} \\
        --tss_window  ${tss_window} \\
        --label       ${peakset_id} \\
        --out_anno    ${peakset_id}_annotation.tsv \\
        --out_pie     ${peakset_id}_anno_pie.pdf \\
        --out_feat    ${peakset_id}_anno_features.pdf \\
        --out_mqc     ${peakset_id}_anno_categories_mqc.tsv
    """
}
