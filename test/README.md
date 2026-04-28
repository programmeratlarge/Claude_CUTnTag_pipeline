# Test bundle

This directory contains everything needed to run a smoke test of the pipeline,
**except** the actual FASTQ data and chromosome FASTA — both of which are too
large to ship in the repo and are organism/specific.

## Files

| File | Purpose |
|------|---------|
| `test_associations.csv`     | Example association CSV: 4 H3K27ac samples (WT/KO × 2 reps), 2 H3K27me3 WT samples (broad mode), 4 IgG controls (WT/KO × 2 reps) |
| `make_synth_test_data.sh`   | One-shot script that subsamples a real FASTQ pair, builds a chr22-only Bowtie2 index, and lays everything out under `test/fastq/` and `test/refs/` |

## How to run a test

You need one real CUT&Tag paired-end FASTQ pair (any sample) and a reference
FASTA (e.g. hg38). On Cornell BioHPC:

```bash
cd /workdir/$USER/cuttag_pipeline
chmod +x test/make_synth_test_data.sh

bash test/make_synth_test_data.sh \
    --r1   /local/storage/data/example/SAMPLE_R1.fastq.gz \
    --r2   /local/storage/data/example/SAMPLE_R2.fastq.gz \
    --hg38 /local/storage/data/short_read_index/hg38/genome.fa \
    --n    10000

module load nextflow/25.4.3
mkdir -p /workdir/$USER/tmp /workdir/$USER/.singularity_cache
export NXF_SINGULARITY_CACHEDIR=/workdir/$USER/.singularity_cache
export TMPDIR=/workdir/$USER/tmp

nextflow run main.nf -profile test,singularity \
    --bowtie2_index   $PWD/test/refs/chr22 \
    --chrom_sizes     $PWD/test/refs/chr22.chrom.sizes \
    --annotation_gtf  $PWD/test/refs/chr22.gtf \
    --tss_bed         $PWD/test/refs/chr22_tss.bed
```

A successful test run takes ~10–20 minutes on a 2-CPU node and exercises every
process at least once. Inspect:

- `test/results/13_multiqc/final/multiqc_final_report.html`
- `test/results/00_fastq_pairs/fastq_pairs.tsv`
- `test/results/00_fastq_pairs/associations_validated.csv`
- `test/results/08_peaks_per_sample/`
- `test/results/09_peaks_merged_groups/`

## What the test exercises

- ✅ FASTQ pairing across 10 samples
- ✅ Association CSV validation (all 8 rules)
- ✅ FastQC on raw + trimmed reads
- ✅ Cutadapt with Tn5 adapters
- ✅ Bowtie2 alignment
- ✅ Filtering: MAPQ + properly paired + mito + (no blacklist)
- ✅ Picard MarkDuplicates
- ✅ bamCoverage → BigWig + bedGraph
- ✅ Per-sample MACS3 with control matching (IgG WT vs WT samples; IgG KO vs KO samples)
- ✅ Per-merge-group MACS3 (4 merge groups × 2 antibodies)
- ✅ Consensus peaks (per antibody)
- ✅ ChIPseeker annotation (per-sample, per-group, per-consensus)
- ✅ FRiP scoring at all three levels
- ✅ deepTools matrices + plotProfile + plotHeatmap
- ✅ All four MultiQC stages (raw, trimmed, alignment, final)

## Failure-mode tests

To verify the pipeline fails gracefully:

1. **Drop one R2 file** from `test/fastq/`:
   ```bash
   rm test/fastq/H3K27ac_WT_R1_R2.fastq.gz
   nextflow run main.nf -profile test,singularity ...
   # → DETECT_FASTQ_PAIRS exits with "no R2 mate found" + non-zero
   ```

2. **Break the CSV**: edit `test_associations.csv` and change one
   `control_group_id` to a value that doesn't exist:
   ```bash
   sed -i 's/grp_IgG_WT/grp_NONEXISTENT/' test/test_associations.csv
   nextflow run main.nf -profile test,singularity ...
   # → VALIDATE_ASSOCIATIONS exits with the rule-6 error
   ```

3. **No-control mode**: keep only the H3K27ac rows in the CSV (delete IgG rows)
   and add `--allow_no_control true`. The pipeline should run, MACS3 should
   call peaks without `-c`, and the FRiP/annotation modules should still
   complete.
