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
    def bt2_args  = params.bowtie2_args ?: '--end-to-end --very-sensitive --no-mixed --no-discordant -I 10 -X 700'
    // Read group: set at alignment time so Picard MarkDuplicates (which requires @RG)
    // works downstream. Fields:
    //   ID = sample_id            (unique per file — Picard / GATK key on this)
    //   SM = sample_id            (sample-level identifier)
    //   LB = sample_id            (library — used to scope optical-duplicate detection)
    //   PL = ILLUMINA             (platform; CUT&Tag is Illumina in practice)
    //   PU = <flowcell>.<lane>    (best-effort: derive from sample_id, else fall back to sample_id)
    def pu = (sample_id =~ /([A-Z0-9]+)_L(\d+)$/)
    def pu_str = pu ? "${pu[0][1]}.${pu[0][2]}" : sample_id
    """
    bowtie2 \\
        -x ${bt2_index_dir}/${bt2_index_name} \\
        -1 ${r1} -2 ${r2} \\
        -p ${task.cpus} \\
        --rg-id ${sample_id} \\
        --rg SM:${sample_id} \\
        --rg LB:${sample_id} \\
        --rg PL:ILLUMINA \\
        --rg PU:${pu_str} \\
        ${bt2_args} \\
        2> ${sample_id}.bowtie2.log \\
      | samtools view -@ ${task.cpus} -b -o ${sample_id}.bam -
    """
}
