/*
 * deepTools processes: bamCoverage, computeMatrix (TSS + peak-centered),
 * plotProfile, plotHeatmap, plus bigWigToBedGraph helper.
 */

process BAMCOVERAGE {
    tag "${sample_id}"
    label 'process_medium'
    publishDir "${params.outdir}/06_bigwig", mode: params.publish_mode

    input:
    tuple val(sample_id), val(meta), path(bam), path(bai)
    val   binsize
    val   norm_method
    path  chrom_sizes

    output:
    tuple val(sample_id), val(meta), path("${sample_id}.bw"), emit: bigwig

    script:
    def norm_arg = (norm_method == 'None' || norm_method == 'none') ?
        '' : "--normalizeUsing ${norm_method}"
    def egs = norm_method == 'RPGC' ? "--effectiveGenomeSize ${params.effective_genome_size ?: 2913022398}" : ''
    """
    bamCoverage \\
        -b ${bam} \\
        -o ${sample_id}.bw \\
        --binSize ${binsize} \\
        ${norm_arg} \\
        ${egs} \\
        --extendReads \\
        --centerReads \\
        --numberOfProcessors ${task.cpus}
    """
}

process BIGWIG_TO_BEDGRAPH {
    tag "${sample_id}"
    label 'process_low'
    publishDir "${params.outdir}/07_bedgraph", mode: params.publish_mode

    input:
    tuple val(sample_id), val(meta), path(bigwig)
    path  chrom_sizes

    output:
    tuple val(sample_id), val(meta), path("${sample_id}.bedgraph"), emit: bedgraph

    script:
    """
    bigWigToBedGraph ${bigwig} ${sample_id}.bedgraph
    """
}

process COMPUTE_MATRIX_TSS {
    label 'process_high'
    publishDir "${params.outdir}/12_deeptools/matrices", mode: params.publish_mode

    input:
    tuple val(sample_ids), path(bigwigs)
    path  tss_bed
    val   tss_window

    output:
    path 'tss_matrix.gz', emit: matrix

    script:
    def labels = sample_ids.join(' ')
    """
    computeMatrix reference-point \\
        --referencePoint TSS \\
        -S ${bigwigs} \\
        -R ${tss_bed} \\
        -b ${tss_window} -a ${tss_window} \\
        --skipZeros \\
        --samplesLabel ${labels} \\
        -o tss_matrix.gz \\
        -p ${task.cpus}
    """
}

process COMPUTE_MATRIX_PEAKS {
    tag "${antibody}"
    label 'process_high'
    publishDir "${params.outdir}/12_deeptools/matrices", mode: params.publish_mode

    input:
    tuple val(sample_ids), path(bigwigs), val(antibody), path(peaks_bed)
    val peak_window

    output:
    tuple val(antibody), path("${antibody}_peaks_matrix.gz"), emit: matrix

    script:
    def labels = sample_ids.join(' ')
    """
    computeMatrix reference-point \\
        --referencePoint center \\
        -S ${bigwigs} \\
        -R ${peaks_bed} \\
        -b ${peak_window} -a ${peak_window} \\
        --skipZeros \\
        --samplesLabel ${labels} \\
        -o ${antibody}_peaks_matrix.gz \\
        -p ${task.cpus}
    """
}

process PLOT_PROFILE {
    label 'process_low'
    publishDir "${params.outdir}/12_deeptools/profiles", mode: params.publish_mode

    input:
    path matrix
    val  tag_

    output:
    path "${tag_}_profile.pdf", emit: pdf
    path "${tag_}_profile.png", emit: png

    script:
    """
    plotProfile -m ${matrix} -o ${tag_}_profile.pdf --plotType lines
    plotProfile -m ${matrix} -o ${tag_}_profile.png --plotType lines
    """
}

process PLOT_HEATMAP {
    label 'process_low'
    publishDir "${params.outdir}/12_deeptools/heatmaps", mode: params.publish_mode

    input:
    path matrix
    val  tag_

    output:
    path "${tag_}_heatmap.pdf", emit: pdf
    path "${tag_}_heatmap.png", emit: png

    script:
    """
    plotHeatmap -m ${matrix} -o ${tag_}_heatmap.pdf --colorMap viridis
    plotHeatmap -m ${matrix} -o ${tag_}_heatmap.png --colorMap viridis
    """
}
