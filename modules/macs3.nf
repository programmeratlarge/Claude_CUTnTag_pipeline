/*
 * MACS3 peak-calling processes for both per-sample and per-merged-group calls.
 * Uses BAMPE format because all inputs here are paired-end CUT&Tag.
 *
 * peak_calling_mode (from association_csv):
 *   'narrow' -> MACS3 default (no --broad)
 *   'broad'  -> MACS3 --broad
 *   'auto'   -> infer from antibody:
 *                 H3K27ac, H3K4me3, CTCF, TFs        => narrow
 *                 H3K27me3, H3K9me3, H3K36me3, H3K4me1, H3K9me2 => broad
 *                 default => narrow
 */

def infer_peak_mode(meta) {
    def m = meta.peak_calling_mode?.toLowerCase()
    if (m == 'broad' || m == 'narrow') return m
    def ab = meta.antibody?.toUpperCase() ?: ''
    if (ab in ['H3K27ME3','H3K9ME3','H3K36ME3','H3K4ME1','H3K9ME2']) return 'broad'
    return 'narrow'
}

process MACS3_PER_SAMPLE {
    tag "${sample_id}"
    label 'process_medium'
    publishDir "${params.outdir}/08_peaks_per_sample", mode: params.publish_mode,
        saveAs: { fn -> fn.endsWith('.log') ? "logs/${fn}" : fn }

    input:
    tuple val(sample_id), val(meta), path(bam), path(bai), path(ctrl_bams)

    output:
    tuple val(sample_id), val(meta), path("${sample_id}_peaks.{narrowPeak,broadPeak}"), emit: peaks
    tuple val(sample_id), path("${sample_id}.macs3.log"),                                emit: log
    path "${sample_id}_*.{xls,bed,bdg}",                                                 optional: true, emit: extras

    script:
    def gsize            = params.genome_size      ?: 'hs'
    def qvalue           = params.macs_qvalue      ?: 0.05
    def allow_no_control = (params.allow_no_control == true || params.allow_no_control == 'true')
    def mode  = infer_peak_mode(meta)
    def broad = mode == 'broad' ? '--broad --broad-cutoff 0.1' : ''
    def ctrl  = (ctrl_bams && ctrl_bams.size() > 0) ? "-c ${ctrl_bams.join(' ')}" : ''
    if (!ctrl && !meta.is_control && !allow_no_control && meta.control_group_id) {
        log.warn "[MACS3_PER_SAMPLE] ${sample_id} has control_group_id=${meta.control_group_id} but no control BAMs were resolved — falling back to no-control call (set --allow_no_control true to silence)"
    }
    """
    macs3 callpeak \\
        -t ${bam} ${ctrl} \\
        -f BAMPE \\
        -g ${gsize} \\
        -n ${sample_id} \\
        -q ${qvalue} \\
        ${broad} \\
        --keep-dup all \\
        --bdg \\
        --SPMR \\
        2> ${sample_id}.macs3.log
    """
}

process MACS3_GROUP {
    tag "${merge_group_id}"
    label 'process_medium'
    publishDir "${params.outdir}/09_peaks_merged_groups", mode: params.publish_mode,
        saveAs: { fn -> fn.endsWith('.log') ? "logs/${fn}" : fn }

    input:
    tuple val(merge_group_id), val(meta), path(treat_bam), path(treat_bai), path(ctrl_bam), path(ctrl_bai)

    output:
    tuple val(merge_group_id), val(meta), path("${merge_group_id}_peaks.{narrowPeak,broadPeak}"), emit: peaks
    tuple val(merge_group_id), path("${merge_group_id}.macs3.log"),                                emit: log
    path  "${merge_group_id}_*.{xls,bed,bdg}", optional: true, emit: extras

    script:
    def gsize  = params.genome_size  ?: 'hs'
    def qvalue = params.macs_qvalue  ?: 0.05
    def mode   = infer_peak_mode(meta)
    def broad  = mode == 'broad' ? '--broad --broad-cutoff 0.1' : ''
    def ctrl   = (ctrl_bam && ctrl_bam.size() > 0) ? "-c ${ctrl_bam}" : ''
    """
    macs3 callpeak \\
        -t ${treat_bam} ${ctrl} \\
        -f BAMPE \\
        -g ${gsize} \\
        -n ${merge_group_id} \\
        -q ${qvalue} \\
        ${broad} \\
        --keep-dup all \\
        --bdg \\
        --SPMR \\
        2> ${merge_group_id}.macs3.log
    """
}

process CONSENSUS_PEAKS {
    tag "${antibody}"
    label 'process_low'
    publishDir "${params.outdir}/09_peaks_merged_groups/consensus", mode: params.publish_mode

    input:
    tuple val(antibody), path(peak_files)
    path  chrom_sizes

    output:
    tuple val(antibody), path("${antibody}_consensus.bed"), emit: consensus

    script:
    """
    cat ${peak_files} \\
      | cut -f1-3 \\
      | sort -k1,1 -k2,2n \\
      | bedtools merge -i - \\
      > ${antibody}_consensus.bed
    """
}
