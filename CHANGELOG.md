# Changelog

All notable changes to `cuttag-dsl2` will be documented in this file.
The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [Unreleased]

### Added
- **Genome filter**: `validate_associations.py` now drops CSV rows whose
  `genome` column doesn't match the `--genome` parameter. This lets a single
  association CSV cover multiple genomes; run the pipeline once per genome
  with the corresponding references. FASTQ files corresponding to filtered
  rows are skipped silently (with an INFO line in the validation log). The
  validator fails loudly if `--genome` matches zero rows, listing the genomes
  it did find — useful for catching `hg38` vs `GRCh38` typos.

### Changed
- `bin/detect_fastq_pairs.py` now searches `--input_dir` **recursively**,
  including all subdirectories (hidden directories starting with `.` are
  skipped). Behaviour is unchanged for flat layouts.
- The pairs.log now lists every subdirectory that contained FASTQ files, and
  duplicate-sample-ID errors list every offending file path (not just the
  derived ID).

### Fixed
- **`A process input channel evaluates to null -- Invalid declaration 'val ...'`** at the `CUTADAPT` call (and would have fired downstream at MACS3, BAMCOVERAGE, etc.). The root cause was Nextflow 23.x having trouble auto-broadcasting multiple trailing `val` inputs against a queue channel that's reused by other processes upstream (e.g. `reads_with_meta_ch` consumed by both `FASTQC_RAW` and `CUTADAPT`). Refactored 7 modules to read `params.*` values directly inside the process script block instead of declaring them as `val` inputs: `CUTADAPT`, `FILTER_MITO_BLACKLIST`, `MACS3_PER_SAMPLE`, `MACS3_GROUP`, `BAMCOVERAGE`, `COMPUTE_MATRIX_TSS`, `COMPUTE_MATRIX_PEAKS`, `ANNOTATE_PEAKS`. Behaviour is unchanged; the value channels are eliminated, so the broadcast issue disappears.
- `find_mate()` now anchors the R1→R2 substitution at the END of the filename.
  Previously, sample names that themselves contained `_R1_` (e.g. a replicate
  marker like `H3K27ac_WT_R1_S1_L001_R1_001.fastq.gz`) would have the WRONG
  occurrence replaced, so the R2 mate was never found.
- Orphan-R2 detection now compares full paths instead of bare filenames, so
  two subdirectories containing R2 files with the same name no longer mask
  each other.
- `conf/base.config` no longer uses `withName: 'X' { label '...' }`, which
  Nextflow rejects with "Unknown method invocation `label`". Labels are
  declared in each module's process body, which is the only valid place.
- `conf/slurm.config` no longer hard-codes `--account=default`; the partition
  and account are now parameterised via `--slurm_partition` / `--slurm_account`
  with defaults in the main params block.

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
