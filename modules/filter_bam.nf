/*
 * Filter BAM:
 *   1) keep only properly paired primary alignments with MAPQ >= threshold
 *   2) remove reads on mitochondrial chromosomes (mito_chroms is comma-separated)
 *   3) optional blacklist removal via bedtools intersect -v
 *
 * Emits per-stage read-count file used downstream by READ_RETENTION.
 */

process FILTER_MITO_BLACKLIST {
    tag "${sample_id}"
    label 'process_medium'
    publishDir "${params.outdir}/05_filtering", mode: params.publish_mode

    input:
    tuple val(sample_id), val(meta), path(bam), path(bai)
    path  blacklist

    output:
    tuple val(sample_id), val(meta),
          path("${sample_id}.filtered.bam"),
          path("${sample_id}.filtered.bam.bai"), emit: bam
    tuple val(sample_id), path("${sample_id}.filter_counts.tsv"), emit: counts

    script:
    def mito_chroms = params.mito_chroms ?: 'chrM,MT,M,Mt,mitochondrion'
    def mapq        = params.mapq         ?: 30
    def mito_pattern = mito_chroms.split(',').collect { "^${it}\$" }.join('|')
    def blacklist_step = (blacklist.name == 'empty.bed' || blacklist.size() == 0) ?
        "cp ${sample_id}.nomito.bam ${sample_id}.preDup.bam" :
        "bedtools intersect -v -abam ${sample_id}.nomito.bam -b ${blacklist} > ${sample_id}.preDup.bam"
    """
    set -euo pipefail

    # 1. Total in input BAM
    aligned=\$(samtools view -c ${bam})

    # 2. MAPQ + properly paired filtering (-f 2 properly paired, -F 0x4 not-unmapped, -F 0x100 not-secondary)
    samtools view -@ ${task.cpus} -b -f 2 -F 0x004 -F 0x100 -F 0x800 -q ${mapq} ${bam} \\
        > ${sample_id}.qfilt.bam
    qfilt=\$(samtools view -c ${sample_id}.qfilt.bam)

    # 3. Remove mitochondrial chromosomes
    samtools index ${sample_id}.qfilt.bam
    chroms=\$(samtools idxstats ${sample_id}.qfilt.bam | cut -f1 | grep -Ev '${mito_pattern}' || true)
    samtools view -@ ${task.cpus} -b ${sample_id}.qfilt.bam \$chroms > ${sample_id}.nomito.bam
    nomito=\$(samtools view -c ${sample_id}.nomito.bam)

    # 4. Blacklist filter (optional)
    ${blacklist_step}
    samtools index ${sample_id}.preDup.bam
    nobl=\$(samtools view -c ${sample_id}.preDup.bam)

    # Final output (still has duplicates — dedup is the next stage)
    mv ${sample_id}.preDup.bam ${sample_id}.filtered.bam
    samtools index ${sample_id}.filtered.bam

    # Per-sample counts TSV (used by READ_RETENTION)
    {
      echo -e "sample_id\\taligned\\tmapq_pp\\tnomito\\tnoblacklist"
      echo -e "${sample_id}\\t\${aligned}\\t\${qfilt}\\t\${nomito}\\t\${nobl}"
    } > ${sample_id}.filter_counts.tsv
    """
}
