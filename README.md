# cuttag-dsl2 — Production Nextflow DSL2 pipeline for paired-end CUT&Tag

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A522.10-23aa62.svg)](https://www.nextflow.io/)
[![Run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg)](https://sylabs.io/singularity/)

A reproducible, container-aware Nextflow DSL2 pipeline that takes paired-end
CUT&Tag FASTQ files all the way through QC → trimming → alignment →
filtering → deduplication → BigWig signal tracks → per-sample peak calling →
group-level merged peak calling (with control matching) → peak annotation →
FRiP → deepTools profile/heatmap visualizations → integrated MultiQC report.

---

## Pipeline overview

```
   raw FASTQ ─▶ FastQC ─▶ Cutadapt ─▶ FastQC ─▶ Bowtie2 ─▶ samtools sort/index
                                                                  │
                                            ┌─────────────────────┘
                                            ▼
                       Filter (MAPQ + properly-paired + mito + blacklist)
                                            │
                                            ▼
                            Picard MarkDuplicates  /  umi_tools dedup
                                            │
                ┌───────────────────────────┼─────────────────────────┐
                ▼                           ▼                         ▼
   bamCoverage → BigWig/BedGraph   MACS3 per sample          Merge BAMs by
                                  (vs sample's control)      merge_group_id
                                                                      │
                                                                      ▼
                                                        MACS3 per merged group
                                                        (vs merged control)
                                                                      │
                                                ┌─────────────────────┴─────────┐
                                                ▼                               ▼
                                       ChIPseeker annotation              Consensus peaks
                                                │                               │
                                                └────────────► FRiP, deepTools matrices/heatmaps
                                                                                │
                                                                                ▼
                                                                       Final MultiQC report
```

---

## Required inputs

| Parameter            | Required | Description                                                      |
|----------------------|:--------:|------------------------------------------------------------------|
| `--input_dir`        | ✅       | Directory containing paired-end FASTQ files                      |
| `--outdir`           | ✅       | Output directory                                                 |
| `--genome`           | ✅       | Genome label (e.g. `hg38`, `mm10`, `rn6`, `custom`)              |
| `--bowtie2_index`    | ✅       | Bowtie2 index basename (e.g. `/path/to/idx/hg38`)                |
| `--chrom_sizes`      | ✅       | Tab-delimited chromosome sizes file                              |
| `--annotation_gtf`   | ✅       | GTF/GFF for ChIPseeker peak annotation                           |
| `--tss_bed`          | ✅       | BED of TSS coordinates (deepTools matrix)                        |
| `--association_csv`  | ✅       | Sample/antibody/control mapping table (see below)                |
| `--blacklist`        | ⬜       | BED of blacklisted regions                                       |
| `--adapter_fwd/rev`  | ⬜       | Override adapter sequences (default Nextera/Tn5)                 |
| `--use_umi`          | ⬜       | Enable UMI-aware deduplication                                   |

Run `nextflow run main.nf --help` for the full parameter list.

---

## Association table specification

The `--association_csv` is the heart of sample-to-control and replicate-to-group
matching. **Required columns**:

| Column              | Type    | Description |
|---------------------|---------|-------------|
| `sample_id`         | string  | Must match the FASTQ-derived sample ID |
| `species`           | string  | `human`, `mouse`, `rat`, `custom`, … |
| `genome`            | string  | `hg38`, `mm10`, `rn6`, … |
| `antibody`          | string  | Target antibody, e.g. `H3K27ac`, `H3K4me3`, `IgG`, `Input` |
| `condition`         | string  | Treatment / cell type / experimental group |
| `replicate`         | string  | Biological replicate ID (`1`, `2`, `R1`, …) |
| `group_id`          | string  | **Unique** sample-group identifier |
| `is_control`        | bool    | `true` if IgG/Input/control, otherwise `false` |
| `control_group_id`  | string  | For non-control rows: the `group_id` of the control to use |
| `merge_group_id`    | string  | Samples with the same value are merged for group-level peak calling |
| `peak_calling_mode` | enum    | `narrow`, `broad`, or `auto` (auto infers from antibody) |
| `notes`             | string  | Optional free text |

### Validation rules (enforced by `validate_associations.py`)

1. Every `sample_id` from the FASTQ pair manifest must appear in the CSV.
2. No extra rows in the CSV unless `--allow_extra true`.
3. `is_control` must be exactly `true` or `false` (case-insensitive).
4. `peak_calling_mode` must be one of `narrow`, `broad`, `auto`, or empty (treated as `auto`).
5. All rows sharing a `merge_group_id` must agree on `species`, `genome`, `antibody`, and `condition`.
6. Every non-control row must have a `control_group_id` matching some `is_control=true` row's `group_id`, **unless** `--allow_no_control true`.
7. `is_control=true` rows must leave `control_group_id` empty (controls don't need controls).

The pipeline fails fast with a clear error message if any rule is violated.

---

## Example `association_csv`

```csv
sample_id,species,genome,antibody,condition,replicate,group_id,is_control,control_group_id,merge_group_id,peak_calling_mode,notes
H3K27ac_WT_R1,human,hg38,H3K27ac,WT,1,grp_H3K27ac_WT_R1,false,grp_IgG_WT,merge_H3K27ac_WT,narrow,
H3K27ac_WT_R2,human,hg38,H3K27ac,WT,2,grp_H3K27ac_WT_R2,false,grp_IgG_WT,merge_H3K27ac_WT,narrow,
H3K27ac_KO_R1,human,hg38,H3K27ac,KO,1,grp_H3K27ac_KO_R1,false,grp_IgG_KO,merge_H3K27ac_KO,narrow,
H3K27ac_KO_R2,human,hg38,H3K27ac,KO,2,grp_H3K27ac_KO_R2,false,grp_IgG_KO,merge_H3K27ac_KO,narrow,
H3K27me3_WT_R1,human,hg38,H3K27me3,WT,1,grp_H3K27me3_WT_R1,false,grp_IgG_WT,merge_H3K27me3_WT,broad,
H3K27me3_WT_R2,human,hg38,H3K27me3,WT,2,grp_H3K27me3_WT_R2,false,grp_IgG_WT,merge_H3K27me3_WT,broad,
IgG_WT_R1,human,hg38,IgG,WT,1,grp_IgG_WT,true,,,auto,IgG control for WT
IgG_WT_R2,human,hg38,IgG,WT,2,grp_IgG_WT,true,,,auto,IgG control for WT
IgG_KO_R1,human,hg38,IgG,KO,1,grp_IgG_KO,true,,,auto,IgG control for KO
IgG_KO_R2,human,hg38,IgG,KO,2,grp_IgG_KO,true,,,auto,IgG control for KO
```

In this example:
- 4 H3K27ac samples are split into 2 merge groups (WT and KO), each compared against the merged IgG controls of its condition.
- 2 H3K27me3 WT samples are merged and compared against merged IgG WT controls; broad peak mode is used.
- IgG samples sharing the same `group_id` are merged into a single control BAM.

---

## Output directory structure

```
results/
├── 00_fastq_pairs/                 fastq_pairs.tsv, validation log, validated CSV
├── 01_fastqc_raw/                  FastQC reports on raw reads
├── 02_trimmed/                     trimmed FASTQs + cutadapt logs
├── 03_fastqc_trimmed/              FastQC reports on trimmed reads
├── 04_alignment/                   sorted BAMs + bowtie2 logs + flagstat/idxstats/stats
├── 05_filtering/                   filtered BAMs, dedup metrics, retention TSVs
│   ├── dedup/                      MarkDuplicates / umi_tools dedup outputs
│   ├── merged_treatment/           merged treatment BAMs
│   └── merged_control/             merged control BAMs
├── 06_bigwig/                      .bw files (normalized signal)
├── 07_bedgraph/                    .bedgraph
├── 08_peaks_per_sample/            MACS3 narrowPeak/broadPeak per sample
├── 09_peaks_merged_groups/         merged-group MACS3 peaks
│   └── consensus/                  consensus peaks per antibody
├── 10_peak_annotation/             ChIPseeker output + plots
├── 11_frip/                        FRiP TSVs (per-sample own/group/consensus, per-group)
├── 12_deeptools/
│   ├── matrices/
│   ├── profiles/
│   ├── heatmaps/
│   └── tss_enrichment/
├── 13_multiqc/                     raw / trimmed / alignment / final MultiQC reports
└── pipeline_info/                  trace, timeline, report, DAG
```

---

## Quick start

### 1. Clone the pipeline

```bash
mkdir -p /workdir/$USER && cd /workdir/$USER
git clone <this-repo> cuttag_pipeline
cd cuttag_pipeline
```

### 2. Prepare reference files

For human hg38, you need:
- Bowtie2 index (e.g. `/local/storage/data/human_hg38/Sequence/Bowtie2Index/genome`)
- `hg38.chrom.sizes`
- `hg38.gtf`
- `hg38_tss.bed` (extract from GTF)
- `hg38-blacklist.v2.bed.gz` (ENCODE blacklist, optional)

### 3. Stage your FASTQs and write your association CSV

Place your `*_R1*.fastq.gz` / `*_R2*.fastq.gz` files anywhere under `--input_dir` — the pipeline **recursively scans all subdirectories** (skipping hidden ones). A common layout is one subdirectory per sequencing run or lane:

```
fastq/
├── runA/
│   ├── lane1/SampleA_S1_L001_R1_001.fastq.gz
│   └── lane1/SampleA_S1_L001_R2_001.fastq.gz
├── runB/
│   ├── SampleB_S2_L002_R1_001.fastq.gz
│   └── SampleB_S2_L002_R2_001.fastq.gz
└── runC/Sample_input.R1.fastq.gz  # different naming convention also OK
```

R1/R2 mates must live in the **same** directory. The recognized naming
patterns are listed in the docstring of `bin/detect_fastq_pairs.py`. Sample IDs
are derived from the R1 filename and must be unique across the whole tree; if
the same sample appears in two subdirectories, `detect_fastq_pairs.py` fails
with a clear error that lists every offending path.

### 4. Launch (Cornell BioHPC example)

```bash
module load nextflow/25.4.3

mkdir -p /workdir/$USER/tmp /workdir/$USER/.singularity_cache
export NXF_SINGULARITY_CACHEDIR=/workdir/$USER/.singularity_cache
export TMPDIR=/workdir/$USER/tmp

nextflow run main.nf \
    --input_dir       /workdir/$USER/cuttag/fastq \
    --outdir          /workdir/$USER/cuttag/results \
    --genome          hg38 \
    --bowtie2_index   /workdir/$USER/refs/hg38/bowtie2/hg38 \
    --chrom_sizes     /workdir/$USER/refs/hg38/hg38.chrom.sizes \
    --annotation_gtf  /workdir/$USER/refs/hg38/hg38.gtf \
    --tss_bed         /workdir/$USER/refs/hg38/hg38_tss.bed \
    --blacklist       /workdir/$USER/refs/hg38/hg38-blacklist.v2.bed \
    --association_csv /workdir/$USER/cuttag/samplesheet.csv \
    -profile          singularity \
    -resume
```

### 5. SLURM execution

```bash
nextflow run main.nf \
    ... (same params as above) \
    -profile singularity,slurm
```

### 6. Conda execution (no containers)

```bash
nextflow run main.nf \
    ... \
    -profile conda
```

---

## Installation

### Cornell BioHPC users

Pre-installed software is documented at <https://biohpc.cornell.edu/lab/userguide.aspx?a=software>. You only need to:

1. `module load nextflow/25.4.3`
2. Set `NXF_SINGULARITY_CACHEDIR` to a path under `/workdir/$USER/`
3. Use `-profile singularity` (Singularity 1.4.0-1.el9 is already installed)

No conda environment is required if you use `-profile singularity`.

### Other systems

- **Nextflow** ≥ 22.10 (DSL2)
- **Singularity / Apptainer** ≥ 3.5 (recommended) **or** Docker **or** Conda
- Optional: SLURM/SGE for HPC

If using `-profile conda`, you'll need conda or mamba in PATH; the pipeline
sets `process.conda` to `environment.yml` and lazily creates the env.

---

## Required software (full list)

| Tool                | Version  | Use                                                  |
|---------------------|----------|------------------------------------------------------|
| FastQC              | 0.12.1   | Per-FASTQ quality                                    |
| MultiQC             | 1.25.1   | Aggregate reports                                    |
| Cutadapt            | 4.9      | Adapter trimming                                     |
| Bowtie2             | 2.5.4    | Read alignment                                       |
| Samtools            | 1.20     | BAM I/O, sorting, filtering                          |
| BEDtools            | 2.31.1   | BED operations, blacklist intersect                  |
| Picard              | 3.2.0    | MarkDuplicates                                       |
| umi_tools           | 1.1.5    | UMI-aware dedup (optional)                           |
| deepTools           | 3.5.5    | bamCoverage, computeMatrix, plotProfile/Heatmap      |
| MACS3               | 3.0.2    | Peak calling                                         |
| ucsc-bigwigtobedgraph |        | BigWig → BedGraph                                    |
| ChIPseeker          | 1.40.0   | Peak annotation (R/Bioconductor)                     |
| Python              | ≥ 3.10   | helper scripts                                       |
| R                   | ≥ 4.4    | annotate_peaks.R                                     |

All of these are pinned in `environment.yml` and the per-process container
images in `conf/{docker,singularity}.config`.

---

## Test dataset / synthetic test strategy

Because real CUT&Tag FASTQs are large, we recommend building a tiny synthetic
test set:

```bash
# 1. Take 10 000 reads from one R1/R2 pair
seqtk sample -s100 sample_R1.fastq.gz 10000 | gzip > test/fastq/H3K27ac_WT_R1_R1.fastq.gz
seqtk sample -s100 sample_R2.fastq.gz 10000 | gzip > test/fastq/H3K27ac_WT_R1_R2.fastq.gz
# repeat for all entries in test/test_associations.csv

# 2. Use a 1-chromosome Bowtie2 index (e.g. chr22)
bowtie2-build chr22.fa test/refs/chr22

# 3. Run with -profile test,singularity
nextflow run main.nf -profile test,singularity \
    --bowtie2_index $PWD/test/refs/chr22 \
    --chrom_sizes   $PWD/test/refs/chr22.chrom.sizes \
    --annotation_gtf $PWD/test/refs/chr22.gtf \
    --tss_bed        $PWD/test/refs/chr22_tss.bed
```

A complete reference test bundle isn't shipped here — it's expected to live
outside Git/registry. See `test/README.md` for instructions.

---

## Recommendations / sensible defaults

| Choice                        | Default                          | Rationale |
|-------------------------------|----------------------------------|-----------|
| Adapter (`adapter_fwd/rev`)   | `CTGTCTCTTATACACATCT`            | Tn5/Nextera adapter from CUT&Tag tagmentation |
| Bowtie2 args                  | `--end-to-end --very-sensitive --no-mixed --no-discordant -I 10 -X 700` | CUT&Tag fragments are 50–500 bp; concordant only |
| MAPQ threshold                | 30                               | Standard CUT&Tag/ChIP filter |
| Mitochondrial chroms          | `chrM,MT,M,Mt,mitochondrion`     | Covers UCSC, Ensembl, NCBI naming |
| Dedup                         | Picard MarkDuplicates (REMOVE)   | UMIs are uncommon in CUT&Tag; switch to `--use_umi` if present |
| MACS3 q-value                 | 0.05                             | Standard threshold |
| BigWig binsize                | 10 bp                            | Sharp resolution at typical CUT&Tag peak widths |
| BigWig normalization          | CPM                              | Comparable across libraries; switch to RPGC for absolute coverage |
| Peak mode (auto)              | broad for H3K27me3/H3K9me3/H3K36me3/H3K4me1/H3K9me2; narrow otherwise | Reflects expected mark biology |
| TSS / peak-center window      | ±3000 bp                         | Standard for profile/heatmap |

---

## Validation and failure-mode handling

- **Missing required parameter** → fail before launch with a help message.
- **Non-existent input file** → fail before launch.
- **Ambiguous FASTQ pairing** → `detect_fastq_pairs.py` exits non-zero and lists offenders.
- **Duplicate sample_id** → `detect_fastq_pairs.py` flags every duplicate row.
- **Missing/extra association_csv rows** → `validate_associations.py` lists every offender.
- **Bad peak_calling_mode / is_control values** → `validate_associations.py` rejects.
- **merge_group_id inconsistency** → `validate_associations.py` rejects.
- **Missing control for non-control sample** → fail unless `--allow_no_control true`.
- **Empty peak set** → annotation/FRiP modules emit empty outputs but don't crash the workflow.

All processes use `errorStrategy = retry` for transient failures (codes 130–145, 104, 125, 137, 139, 140) up to `maxRetries=2`. `-resume` works at every stage because outputs go to deterministic publish directories.

### Troubleshooting

- **`ERROR ~ Unknown method invocation `label` on _parse_closure...`** — caused by setting `label '...'` inside a config `withName` block. `label` is a process-body directive only; in config, override resources with `cpus`/`memory`/`time` directly, or use `withLabel: process_X { ... }` against the labels already declared in each module. The shipped `conf/base.config` does this correctly; if you've added a custom `withName` block and copied the wrong pattern, that's the cause.
- **`Singularity command line option(s) not recognized` on BioHPC** — make sure `module load nextflow/25.4.3` is run before `nextflow run`. BioHPC's bundled Singularity 1.4 is found automatically once Nextflow is on `PATH`.
- **`Cannot invoke method ... on null object` referencing `params.outdir`** — you're missing `--outdir` (or one of the other required flags). The pre-flight in `main.nf` should catch this with a clear message before that error appears; if it doesn't, you've probably edited `main.nf` and removed the validation block.
- **Stuck at `DETECT_FASTQ_PAIRS`** — check `results/00_fastq_pairs/fastq_pairs.log`. Common causes: filenames don't match any of the recognized R1/R2 patterns, or two samples share the same derived `sample_id` after lane-suffix stripping.

---

## Testing strategy

1. **Schema test** — edit `test/test_associations.csv` to introduce a deliberate
   error (e.g. unmatched `control_group_id`) and confirm the pipeline aborts at
   the `VALIDATE_ASSOCIATIONS` stage.
2. **Pairing test** — drop one R2 file from the test directory; confirm
   `DETECT_FASTQ_PAIRS` warns and exits non-zero.
3. **Control-free run** — set `--allow_no_control true` with no IgG samples and
   confirm peaks are still called.
4. **Resume test** — interrupt mid-run; rerun with `-resume` and confirm only
   downstream stages re-execute.
5. **Profile-switch test** — same data on `-profile singularity` and `-profile conda` should yield identical peak BED files (modulo non-deterministic
   tie-breaking).

---

## Version history

- **v1.0.0** (2026-04) — Initial production release.

---

## Citation

If you use this pipeline, please cite:
- Kaya-Okur HS et al. (2019). CUT&Tag for efficient epigenomic profiling of
  small samples and single cells. *Nat. Commun.* 10:1930.
- DiTommaso P. et al. (2017). Nextflow enables reproducible computational workflows. *Nat. Biotechnol.* 35:316–319.
- And the individual tool publications (Bowtie2, MACS3, deepTools, ChIPseeker, MultiQC, etc.).

---

## Documentation

- `README.md` — this file (overview + quick start).
- `docs/REFERENCES.md` — how to prepare reference files (Bowtie2 index, GTF, TSS BED, blacklist) on Cornell BioHPC.
- `docs/DESIGN_NOTES.md` — channel topology, control-resolution logic, and where to make common changes.
- `test/README.md` — running the smoke test.
- `CHANGELOG.md` — version history.

The pipeline ships a small CLI helper for the `--tss_bed` requirement:

```bash
bin/make_tss_bed.py --gtf hg38.gtf --out hg38_tss.bed --protein_coding --one_per_gene
```

---

## License

MIT — see `LICENSE`.
