# Changelog

All notable changes to `cuttag-dsl2` will be documented in this file.
The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [1.0.0] — 2026-04

### Added
- Initial production release.
- Full DSL2 workflow: FASTQ pairing → FastQC → Cutadapt → FastQC → Bowtie2 →
  filter (MAPQ + properly paired + mito + blacklist) → Picard MarkDuplicates
  (or `umi_tools dedup` with `--use_umi`) → bamCoverage/BigWig + bedGraph →
  per-sample MACS3 (with control matching) → group merge → group MACS3 →
  ChIPseeker annotation → FRiP (sample-own + sample-vs-group + group-level) →
  deepTools TSS + peak-centered matrices/profiles/heatmaps → 4-stage MultiQC.
- 16 process modules and 5 helper scripts (4 Python, 1 R) in `bin/`.
- 5 profiles: `standard` / `slurm` / `singularity` / `docker` / `apptainer` /
  `conda` / `test` / `debug`.
- 8-rule association-CSV validator with informative error messages.
- Auto narrow/broad inference from antibody name (broad: H3K27me3, H3K9me3,
  H3K36me3, H3K4me1, H3K9me2; narrow otherwise).
- `--allow_no_control` escape hatch for IgG-free designs.
- `--use_umi` toggle to swap Picard MarkDuplicates for UMI-aware dedup.
- `bin/make_tss_bed.py` helper to derive a TSS BED from any GTF/GFF.
- ENCODE blacklist support (optional, off by default).
- Custom-content MultiQC sections: read retention, FRiP per-sample/per-group,
  peak annotation distribution.
- Docs: `README.md`, `docs/DESIGN_NOTES.md`, `docs/REFERENCES.md`.
- Container images pinned to specific versions in
  `conf/{docker,singularity}.config` and matching the conda spec in
  `environment.yml`.
- BioHPC launch script `run_cuttag_pipeline.sh` with TODO-marked path
  placeholders, scratch under `/workdir/$USER/`.

### Notes
- Default adapter is Tn5/Nextera (`CTGTCTCTTATACACATCT`).
- Default Bowtie2 args: `--end-to-end --very-sensitive --no-mixed --no-discordant -I 10 -X 700`.
- Default MAPQ filter: 30.
- Default BigWig normalization: CPM (binsize 10 bp).
- Default `effective_genome_size` is hg38 (2,913,022,398); override for other
  organisms when using `--bigwig_norm RPGC`.
