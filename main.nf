#!/usr/bin/env nextflow
/*
========================================================================================
    CUT&Tag DSL2 Pipeline — main.nf
    -----------------------------------------------------------------------------
    Paired-end CUT&Tag analysis: QC -> trim -> align -> filter -> peaks -> annotate
    Author : Senior Bioinformatics Pipeline Architect
    License: MIT
========================================================================================
*/

nextflow.enable.dsl = 2

// ----------------------------------------------------------------------------
// Help message
// ----------------------------------------------------------------------------
def helpMessage() {
    log.info """
    ============================================================
      CUT&Tag DSL2 Pipeline
    ============================================================
    Usage:
      nextflow run main.nf \\
          --input_dir         path/to/fastq/dir \\
          --outdir            results \\
          --genome            hg38 \\
          --bowtie2_index     /path/to/bt2/index_basename \\
          --chrom_sizes       hg38.chrom.sizes \\
          --annotation_gtf    hg38.gtf \\
          --tss_bed           hg38_tss.bed \\
          --association_csv   samplesheet.csv \\
          [--blacklist hg38.blacklist.bed] \\
          [--adapter_fwd CTGTCTCTTATACACATCT] \\
          [--adapter_rev CTGTCTCTTATACACATCT] \\
          [--threads 8] \\
          -profile singularity

    Required:
      --input_dir          Directory containing paired-end FASTQ files
      --outdir             Output directory
      --genome             Genome identifier (e.g. hg38, mm10, rn6, custom)
      --bowtie2_index      Bowtie2 index basename
      --chrom_sizes        Chromosome sizes file
      --annotation_gtf     GTF/GFF annotation for ChIPseeker / HOMER
      --tss_bed            BED of TSS coordinates for deepTools enrichment
      --association_csv    Sample/antibody/control mapping (see README)

    Optional:
      --blacklist          BED of blacklisted regions
      --adapter_fwd        Forward adapter (default Nextera)
      --adapter_rev        Reverse adapter (default Nextera)
      --paired_pattern     R1/R2 detection regex (default auto)
      --threads            Default thread count (default 8)
      --use_umi            Enable UMI-aware deduplication (default false)
      --umi_pattern        UMI pattern for umi_tools (default 'NNNNNN')
      --mapq               MAPQ threshold (default 30)
      --mito_chroms        Mitochondrial chromosome names CSV (default chrM,MT,M)
      --macs_qvalue        MACS3 q-value (default 0.05)
      --bigwig_binsize     BigWig bin size (default 10)
      --bigwig_norm        BigWig normalization (CPM,RPGC,RPKM,BPM,None) (default CPM)
      --tss_window         TSS window for deepTools matrix (default 3000)
      --peak_window        Peak-center window for heatmaps (default 3000)
      --allow_no_control   Allow MACS3 without control (default false)
      -profile             singularity | docker | conda | slurm | test
    """.stripIndent()
}

if (params.help) {
    helpMessage()
    exit 0
}

// ----------------------------------------------------------------------------
// Parameter validation
// ----------------------------------------------------------------------------
def required = [
    'input_dir'       : params.input_dir,
    'outdir'          : params.outdir,
    'genome'          : params.genome,
    'bowtie2_index'   : params.bowtie2_index,
    'chrom_sizes'     : params.chrom_sizes,
    'annotation_gtf'  : params.annotation_gtf,
    'tss_bed'         : params.tss_bed,
    'association_csv' : params.association_csv,
]
def missing = required.findAll { k, v -> v == null || v == '' || v == false }
if (missing) {
    log.error "Missing required parameters: ${missing.keySet().join(', ')}"
    helpMessage()
    exit 1
}

// File existence checks (Bowtie2 index basename excluded - it's a prefix, not a single file)
[
    'input_dir'      : params.input_dir,
    'chrom_sizes'    : params.chrom_sizes,
    'annotation_gtf' : params.annotation_gtf,
    'tss_bed'        : params.tss_bed,
    'association_csv': params.association_csv,
].each { name, p ->
    if (!file(p).exists()) {
        log.error "[${name}] file/dir does not exist: ${p}"
        exit 1
    }
}
if (params.blacklist && !file(params.blacklist).exists()) {
    log.error "[blacklist] file does not exist: ${params.blacklist}"; exit 1
}

// ----------------------------------------------------------------------------
// Module imports
// ----------------------------------------------------------------------------
include { DETECT_FASTQ_PAIRS    } from './modules/fastq_pairs.nf'
include { VALIDATE_ASSOCIATIONS } from './modules/validate_associations.nf'
include { FASTQC as FASTQC_RAW  } from './modules/fastqc.nf'
include { FASTQC as FASTQC_TRIM } from './modules/fastqc.nf'
include { CUTADAPT              } from './modules/cutadapt.nf'
include { BOWTIE2_ALIGN         } from './modules/bowtie2.nf'
include { SAMTOOLS_SORT_INDEX   } from './modules/samtools.nf'
include { SAMTOOLS_FLAGSTAT     } from './modules/samtools.nf'
include { SAMTOOLS_IDXSTATS     } from './modules/samtools.nf'
include { SAMTOOLS_STATS        } from './modules/samtools.nf'
include { FILTER_MITO_BLACKLIST } from './modules/filter_bam.nf'
include { MARK_DUPLICATES       } from './modules/markduplicates.nf'
include { UMI_DEDUP             } from './modules/umi_tools.nf'
include { READ_RETENTION        } from './modules/read_retention.nf'
include { BAMCOVERAGE           } from './modules/deeptools.nf'
include { COMPUTE_MATRIX_TSS    } from './modules/deeptools.nf'
include { COMPUTE_MATRIX_PEAKS  } from './modules/deeptools.nf'
include { PLOT_PROFILE          } from './modules/deeptools.nf'
include { PLOT_HEATMAP          } from './modules/deeptools.nf'
include { BIGWIG_TO_BEDGRAPH    } from './modules/deeptools.nf'
include { MACS3_PER_SAMPLE      } from './modules/macs3.nf'
include { MERGE_BAMS as MERGE_TREAT } from './modules/merge_bams.nf'
include { MERGE_BAMS as MERGE_CTRL  } from './modules/merge_bams.nf'
include { MACS3_GROUP           } from './modules/macs3.nf'
include { CONSENSUS_PEAKS       } from './modules/macs3.nf'
include { ANNOTATE_PEAKS        } from './modules/peak_annotation.nf'
include { FRIP_SAMPLE           } from './modules/frip.nf'
include { FRIP_GROUP            } from './modules/frip.nf'
include { MULTIQC as MULTIQC_RAW   } from './modules/multiqc.nf'
include { MULTIQC as MULTIQC_TRIM  } from './modules/multiqc.nf'
include { MULTIQC as MULTIQC_ALIGN } from './modules/multiqc.nf'
include { MULTIQC as MULTIQC_FINAL } from './modules/multiqc.nf'

// ----------------------------------------------------------------------------
// Workflow
// ----------------------------------------------------------------------------
workflow {

    // -----------------------------
    // Stage 0: pair FASTQs + validate associations
    // -----------------------------
    DETECT_FASTQ_PAIRS(
        file(params.input_dir),
        params.paired_pattern ?: 'auto'
    )

    VALIDATE_ASSOCIATIONS(
        DETECT_FASTQ_PAIRS.out.pairs_tsv,
        file(params.association_csv),
        params.allow_no_control ? 'true' : 'false'
    )

    // Read validated metadata into channels
    sample_meta_ch = VALIDATE_ASSOCIATIONS.out.validated_csv
        .splitCsv(header: true)
        .map { row ->
            def meta = [
                sample_id        : row.sample_id,
                species          : row.species,
                genome           : row.genome,
                antibody         : row.antibody,
                condition        : row.condition,
                replicate        : row.replicate,
                group_id         : row.group_id,
                is_control       : (row.is_control?.toLowerCase() == 'true'),
                control_group_id : row.control_group_id,
                merge_group_id   : row.merge_group_id,
                peak_calling_mode: row.peak_calling_mode ?: 'auto',
                notes            : row.notes ?: ''
            ]
            tuple(meta.sample_id, meta)
        }

    // FASTQ pairs channel from validated TSV
    raw_reads_ch = DETECT_FASTQ_PAIRS.out.pairs_tsv
        .splitCsv(header: true, sep: '\t')
        .map { row -> tuple(row.sample_id, file(row.r1), file(row.r2)) }

    // Join meta + reads -> (sample_id, meta, R1, R2)
    reads_with_meta_ch = raw_reads_ch
        .join(sample_meta_ch, by: 0)
        .map { sid, r1, r2, meta -> tuple(sid, meta, r1, r2) }

    // -----------------------------
    // Stage 1: raw FASTQ QC
    // -----------------------------
    FASTQC_RAW(reads_with_meta_ch.map { sid, meta, r1, r2 -> tuple(sid, [r1, r2]) }, 'raw')

    MULTIQC_RAW(
        FASTQC_RAW.out.zips.map { sid, files -> files }.flatten().collect(),
        file("${projectDir}/assets/multiqc_config.yaml"),
        'raw',
        Channel.empty().collect().ifEmpty([])
    )

    // -----------------------------
    // Stage 2: adapter trimming
    // -----------------------------
    CUTADAPT(
        reads_with_meta_ch,
        params.adapter_fwd,
        params.adapter_rev,
        params.min_read_length
    )

    // FastQC on trimmed reads
    FASTQC_TRIM(
        CUTADAPT.out.trimmed.map { sid, meta, r1, r2 -> tuple(sid, [r1, r2]) },
        'trimmed'
    )

    MULTIQC_TRIM(
        FASTQC_TRIM.out.zips.map { sid, files -> files }.flatten().mix(
            CUTADAPT.out.log.map { sid, meta, f -> f }
        ).collect(),
        file("${projectDir}/assets/multiqc_config.yaml"),
        'trimmed',
        Channel.empty().collect().ifEmpty([])
    )

    // -----------------------------
    // Stage 3: alignment
    // -----------------------------
    BOWTIE2_ALIGN(
        CUTADAPT.out.trimmed,
        file("${params.bowtie2_index}").parent,
        file("${params.bowtie2_index}").name
    )

    SAMTOOLS_SORT_INDEX(BOWTIE2_ALIGN.out.bam, 'aligned')
    SAMTOOLS_FLAGSTAT(SAMTOOLS_SORT_INDEX.out.bam, 'aligned')
    SAMTOOLS_IDXSTATS(SAMTOOLS_SORT_INDEX.out.bam, 'aligned')
    SAMTOOLS_STATS(SAMTOOLS_SORT_INDEX.out.bam, 'aligned')

    // -----------------------------
    // Stage 4: filter -> dedup -> final BAM
    // -----------------------------
    FILTER_MITO_BLACKLIST(
        SAMTOOLS_SORT_INDEX.out.bam,
        params.blacklist ? file(params.blacklist) : file("${projectDir}/assets/empty.bed"),
        params.mito_chroms,
        params.mapq
    )

    // UMI dedup or Picard MarkDuplicates
    if (params.use_umi) {
        UMI_DEDUP(FILTER_MITO_BLACKLIST.out.bam)
        final_bam_ch = UMI_DEDUP.out.bam
        dedup_log_ch = UMI_DEDUP.out.log
    } else {
        MARK_DUPLICATES(FILTER_MITO_BLACKLIST.out.bam)
        final_bam_ch = MARK_DUPLICATES.out.bam
        dedup_log_ch = MARK_DUPLICATES.out.metrics
    }

    // Read retention summary across all stages
    READ_RETENTION(
        BOWTIE2_ALIGN.out.log
            .join(FILTER_MITO_BLACKLIST.out.counts, by: 0)
            .join(dedup_log_ch.map { sid, meta, f -> tuple(sid, f) }, by: 0)
            .join(CUTADAPT.out.log.map { sid, meta, f -> tuple(sid, f) }, by: 0)
    )

    // -----------------------------
    // Stage 5: BigWig + BedGraph
    // -----------------------------
    BAMCOVERAGE(
        final_bam_ch,
        params.bigwig_binsize,
        params.bigwig_norm,
        file(params.chrom_sizes)
    )

    BIGWIG_TO_BEDGRAPH(BAMCOVERAGE.out.bigwig, file(params.chrom_sizes))

    // -----------------------------
    // Stage 6: per-sample MACS3 peak calling
    // -----------------------------
    // Build channel for per-sample peak calling. Treatment samples may have a control_group_id.
    // First group control BAMs by their group_id (control_group_id reference).
    ctrl_bams_by_group_ch = final_bam_ch
        .filter { sid, meta, bam, bai -> meta.is_control }
        .map    { sid, meta, bam, bai -> tuple(meta.group_id, bam) }
        .groupTuple()

    // Treatment samples (not control) joined with control bam list
    treat_with_ctrl_ch = final_bam_ch
        .filter { sid, meta, bam, bai -> !meta.is_control }
        .map    { sid, meta, bam, bai -> tuple(meta.control_group_id ?: '__none__', sid, meta, bam, bai) }
        .combine(ctrl_bams_by_group_ch.ifEmpty(tuple('__none__', [])), by: 0)
        .map    { ctrl_grp, sid, meta, bam, bai, ctrl_bams ->
            tuple(sid, meta, bam, bai, ctrl_bams ?: [])
        }

    // Controls also get peaks called against themselves (no control)
    ctrl_for_peakcall_ch = final_bam_ch
        .filter { sid, meta, bam, bai -> meta.is_control }
        .map    { sid, meta, bam, bai -> tuple(sid, meta, bam, bai, []) }

    macs_input_ch = treat_with_ctrl_ch.mix(ctrl_for_peakcall_ch)

    MACS3_PER_SAMPLE(
        macs_input_ch,
        params.genome_size,
        params.macs_qvalue,
        params.allow_no_control ? 'true' : 'false'
    )

    // -----------------------------
    // Stage 7: group-level merge + group peak calling
    // -----------------------------
    // Merge treatment BAMs by merge_group_id
    treat_to_merge_ch = final_bam_ch
        .filter { sid, meta, bam, bai -> !meta.is_control && meta.merge_group_id }
        .map    { sid, meta, bam, bai -> tuple(meta.merge_group_id, meta, bam) }
        .groupTuple()
        .map    { mg, metas, bams -> tuple(mg, metas[0], bams) }

    MERGE_TREAT(treat_to_merge_ch, 'treatment')

    // Merge control BAMs by control's group_id (which is referenced via control_group_id)
    ctrl_to_merge_ch = final_bam_ch
        .filter { sid, meta, bam, bai -> meta.is_control }
        .map    { sid, meta, bam, bai -> tuple(meta.group_id, meta, bam) }
        .groupTuple()
        .map    { gid, metas, bams -> tuple(gid, metas[0], bams) }

    MERGE_CTRL(ctrl_to_merge_ch, 'control')

    // Map each merge_group_id -> control_group_id (take first sample's control_group_id)
    mg_to_ctrlgrp_ch = sample_meta_ch
        .map { sid, meta -> tuple(meta.merge_group_id, meta.control_group_id ?: '__none__', meta.is_control) }
        .filter { mg, cg, isctrl -> !isctrl && mg }
        .unique()
        .map { mg, cg, isctrl -> tuple(mg, cg) }

    // Join merged treatment BAM with control BAM
    treat_merged_ch = MERGE_TREAT.out.bam
        .map { mg, meta, bam, bai -> tuple(mg, meta, bam, bai) }

    ctrl_merged_ch = MERGE_CTRL.out.bam
        .map { gid, meta, bam, bai -> tuple(gid, bam, bai) }

    // (mg, meta, treat_bam, treat_bai) JOIN (mg, cg) -> (mg, meta, treat_bam, treat_bai, cg)
    treat_with_cg_ch = treat_merged_ch.join(mg_to_ctrlgrp_ch, by: 0)
    // Now we need to attach the control merged BAM by cg key
    // Re-key by cg
    treat_keyed_by_cg = treat_with_cg_ch
        .map { mg, meta, tbam, tbai, cg -> tuple(cg, mg, meta, tbam, tbai) }

    group_macs_input_ch = treat_keyed_by_cg
        .combine(ctrl_merged_ch.mix(Channel.of(tuple('__none__', [], []))), by: 0)
        .map { cg, mg, meta, tbam, tbai, cbam, cbai ->
            tuple(mg, meta, tbam, tbai, cbam ?: [], cbai ?: [])
        }

    MACS3_GROUP(
        group_macs_input_ch,
        params.genome_size,
        params.macs_qvalue,
        params.allow_no_control ? 'true' : 'false'
    )

    // Consensus peaks across all merged groups (per antibody)
    peaks_per_antibody_ch = MACS3_GROUP.out.peaks
        .map { mg, meta, peaks -> tuple(meta.antibody, peaks) }
        .groupTuple()

    CONSENSUS_PEAKS(peaks_per_antibody_ch, file(params.chrom_sizes))

    // -----------------------------
    // Stage 8: peak annotation (per-sample + per-group)
    // -----------------------------
    sample_peaks_for_annot = MACS3_PER_SAMPLE.out.peaks
        .map { sid, meta, peaks -> tuple("sample_${sid}", peaks) }

    group_peaks_for_annot = MACS3_GROUP.out.peaks
        .map { mg, meta, peaks -> tuple("group_${mg}", peaks) }

    consensus_peaks_for_annot = CONSENSUS_PEAKS.out.consensus
        .map { ab, peaks -> tuple("consensus_${ab}", peaks) }

    all_peaks_for_annot = sample_peaks_for_annot
        .mix(group_peaks_for_annot)
        .mix(consensus_peaks_for_annot)

    ANNOTATE_PEAKS(
        all_peaks_for_annot,
        file(params.annotation_gtf),
        params.tss_window
    )

    // -----------------------------
    // Stage 9: FRiP — per sample (vs own peaks, vs group peaks, vs consensus)
    // -----------------------------
    // Per-sample FRiP vs own peaks
    sample_frip_own = final_bam_ch
        .filter { sid, meta, bam, bai -> !meta.is_control }
        .join(MACS3_PER_SAMPLE.out.peaks.map { sid, meta, peaks -> tuple(sid, peaks) }, by: 0)
        .map { sid, meta, bam, bai, peaks -> tuple(sid, meta, bam, bai, peaks, 'own') }

    // Per-sample FRiP vs merged-group peaks (key by merge_group_id)
    sample_frip_group = final_bam_ch
        .filter { sid, meta, bam, bai -> !meta.is_control && meta.merge_group_id }
        .map { sid, meta, bam, bai -> tuple(meta.merge_group_id, sid, meta, bam, bai) }
        .combine(MACS3_GROUP.out.peaks.map { mg, meta, peaks -> tuple(mg, peaks) }, by: 0)
        .map { mg, sid, meta, bam, bai, peaks -> tuple(sid, meta, bam, bai, peaks, 'group') }

    FRIP_SAMPLE(sample_frip_own.mix(sample_frip_group))

    // Group-level FRiP: merged BAM vs its merged peaks
    group_frip_in = MERGE_TREAT.out.bam
        .map { mg, meta, bam, bai -> tuple(mg, meta, bam, bai) }
        .join(MACS3_GROUP.out.peaks.map { mg, meta, peaks -> tuple(mg, peaks) }, by: 0)
        .map { mg, meta, bam, bai, peaks -> tuple(mg, meta, bam, bai, peaks) }

    FRIP_GROUP(group_frip_in)

    // -----------------------------
    // Stage 10: deepTools matrices and plots (TSS + peak-centered)
    // -----------------------------
    bw_collected_ch = BAMCOVERAGE.out.bigwig
        .map { sid, meta, bw -> tuple(sid, bw) }
        .collect(flat: false)
        .map { items ->
            def sids = items.collect { it[0] }
            def bws  = items.collect { it[1] }
            tuple(sids, bws)
        }

    COMPUTE_MATRIX_TSS(
        bw_collected_ch,
        file(params.tss_bed),
        params.tss_window
    )

    PLOT_PROFILE(COMPUTE_MATRIX_TSS.out.matrix, 'tss')
    PLOT_HEATMAP(COMPUTE_MATRIX_TSS.out.matrix, 'tss')

    // Peak-centered: use consensus peaks per antibody
    peak_center_input_ch = bw_collected_ch
        .combine(CONSENSUS_PEAKS.out.consensus)
        .map { sids, bws, ab, peaks -> tuple(sids, bws, ab, peaks) }

    COMPUTE_MATRIX_PEAKS(peak_center_input_ch, params.peak_window)

    // -----------------------------
    // Stage 11: final MultiQC integrated report
    // -----------------------------
    align_qc_ch = SAMTOOLS_FLAGSTAT.out.report.map { sid, f -> f }
        .mix(SAMTOOLS_IDXSTATS.out.report.map { sid, f -> f })
        .mix(SAMTOOLS_STATS.out.report.map    { sid, f -> f })
        .mix(BOWTIE2_ALIGN.out.log.map        { sid, f -> f })
        .collect()

    MULTIQC_ALIGN(
        align_qc_ch,
        file("${projectDir}/assets/multiqc_config.yaml"),
        'alignment',
        Channel.empty().collect().ifEmpty([])
    )

    // Aggregate all QC inputs for the FINAL report (paths only — strip sample IDs)
    final_qc_ch = FASTQC_RAW.out.zips.map        { sid, files -> files }.flatten()
        .mix(FASTQC_TRIM.out.zips.map            { sid, files -> files }.flatten())
        .mix(CUTADAPT.out.log.map                { sid, meta, f -> f })
        .mix(BOWTIE2_ALIGN.out.log.map           { sid, f -> f })
        .mix(SAMTOOLS_FLAGSTAT.out.report.map    { sid, f -> f })
        .mix(SAMTOOLS_IDXSTATS.out.report.map    { sid, f -> f })
        .mix(SAMTOOLS_STATS.out.report.map       { sid, f -> f })
        .mix(dedup_log_ch.map                    { sid, meta, f -> f })
        .mix(FILTER_MITO_BLACKLIST.out.counts.map { sid, f -> f })
        .collect()

    custom_content_ch = READ_RETENTION.out.mqc
        .mix(FRIP_SAMPLE.out.mqc)
        .mix(FRIP_GROUP.out.mqc)
        .mix(ANNOTATE_PEAKS.out.mqc)
        .collect()

    MULTIQC_FINAL(
        final_qc_ch,
        file("${projectDir}/assets/multiqc_config.yaml"),
        'final',
        custom_content_ch
    )
}

// ----------------------------------------------------------------------------
// Completion handler
// ----------------------------------------------------------------------------
workflow.onComplete {
    log.info """
    ============================================================
      Pipeline finished: ${workflow.success ? 'SUCCESS' : 'FAILED'}
      Duration : ${workflow.duration}
      Output   : ${params.outdir}
      Reports  : ${params.outdir}/13_multiqc/
    ============================================================
    """.stripIndent()
}

workflow.onError {
    log.error "Pipeline error: ${workflow.errorMessage}"
}
