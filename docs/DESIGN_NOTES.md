# Design notes — `cuttag-dsl2`

This document explains the channel topology, key design decisions, and where to
make changes when extending the pipeline.

---

## 1. Why DSL2 and why this layout

- **DSL2** lets each tool live in its own `modules/*.nf` file with explicit
  `input:` / `output:` blocks, which makes the workflow graph readable.
- The orchestration logic lives in **one place** (`main.nf`); modules are
  agnostic to upstream/downstream context, so swapping (e.g.) `bowtie2` for
  `bwa-mem2` is a single-file change.
- Helper scripts go in `bin/` because Nextflow auto-prepends the pipeline's
  `bin/` to `PATH` inside every task — no need for absolute paths.
- Container assignments are factored into per-engine config files
  (`conf/{docker,singularity}.config`). The same image versions appear in both;
  the conda spec in `environment.yml` is the single source of truth for tool
  versions.

---

## 2. Channel topology

The most subtle piece is how sample → control matching is resolved. Three
patterns are reused throughout `main.nf`:

### 2a. `(sample_id, meta)` join key

Almost every channel that carries per-sample data is keyed on `sample_id`. The
canonical "rich" record is `(sample_id, meta, bam, bai)` where `meta` is a
groovy `Map` parsed once from `associations_validated.csv`. Whenever a process
emits something that needs to be re-joined later, it emits `(sample_id, ...)`
so a `.join(..., by: 0)` works.

### 2b. Control resolution at per-sample peak calling

Controls are pooled by their own `group_id`, then joined to treatment samples
via `control_group_id`. Conceptually:

```
controls   : (group_id, [bam1, bam2, ...])    # via groupTuple
treatments : (control_group_id, sid, meta, bam, bai)
combine by control_group_id == group_id
=>           (sid, meta, bam, bai, [ctrl_bam1, ...])
```

A sentinel `'__none__'` is used so treatments without controls can still join
(via an `ifEmpty(tuple('__none__', []))`) without dropping rows.

### 2c. Group merging — two-level keying

Treatment BAMs are grouped on `merge_group_id`; control BAMs are grouped on the
controls' own `group_id`. The two are joined through a small mapping channel
`mg_to_ctrlgrp_ch` derived from the validated CSV that translates one to the
other.

```
mg_to_ctrlgrp_ch : (merge_group_id, control_group_id)
treat_merged_ch  : (merge_group_id, meta, treat_bam, treat_bai)
ctrl_merged_ch   : (control_group_id,        ctrl_bam, ctrl_bai)

step 1 : treat_merged_ch.join(mg_to_ctrlgrp_ch, by: 0)
         => (mg, meta, tbam, tbai, cg)

step 2 : re-key by cg, then combine with ctrl_merged_ch (also keyed by cg)
         => (cg, mg, meta, tbam, tbai, cbam, cbai)
```

This is why `MACS3_GROUP` consumes a 7-tuple: it carries both keys plus the
metadata needed for narrow/broad inference.

---

## 3. Why FRiP runs three times per sample

The pipeline reports FRiP against:

1. The sample's **own** narrow/broad peak set — measures self-consistency.
2. The **merged-group** peak set the sample contributed to — measures how
   "typical" the sample is for its replicate group.
3. The **merged-group BAM** vs. its own peak set (`FRIP_GROUP`) — gives a
   single number per replicate group, useful for cross-condition comparisons.

A consensus-peak FRiP isn't computed by default because that would weight
ubiquitous peaks heavily and isn't directly meaningful for sample-level QC.
You can add it by extending `sample_frip_*` in `main.nf` to also `.combine`
with `CONSENSUS_PEAKS.out.consensus`.

---

## 4. Why narrow/broad is inferred from antibody name

The CSV's `peak_calling_mode` column is the source of truth, but most users
leave it as `auto`. The default rule (in `modules/macs3.nf::infer_peak_mode`):

| Antibody                                                     | Mode    |
|--------------------------------------------------------------|---------|
| `H3K27me3`, `H3K9me3`, `H3K36me3`, `H3K4me1`, `H3K9me2`      | broad   |
| Anything else                                                | narrow  |

This reflects the standard CUT&Tag/ChIP literature (broad marks vs. point-source
marks). Override per-sample by setting `peak_calling_mode=broad|narrow` in the
CSV.

---

## 5. Why two MACS3 calls per merged group

Merged-group MACS3 is called once with the merged control BAM as `-c`. There is
intentionally **no** "merged treatment vs. own merged treatment" call; that
would just be calling peaks against itself. If you want a no-control reference
call, set `--allow_no_control true` and ensure the control_group_id column is
empty (or use `__none__` semantics).

---

## 6. Resource scaling and retries

- All processes use labels (`process_low/medium/high/long`) with resources that
  scale on `task.attempt` (4 GB → 8 GB → 16 GB on retry, etc.).
- Transient cluster signals (130–145, 104, 125, 137, 139, 140) trigger an
  automatic retry up to `maxRetries=2`. Persistent failures (exit 1, 2)
  terminate the run.
- `check_max()` in `nextflow.config` caps any computed resource at the
  user-supplied `--max_cpus`, `--max_memory`, `--max_time`. This makes the same
  pipeline portable from a laptop to a 48-core HPC node without editing
  module-level defaults.

---

## 7. Where to make common changes

| You want to…                                  | Edit                                       |
|-----------------------------------------------|--------------------------------------------|
| Change Bowtie2 args                           | `params.bowtie2_args` in `nextflow.config` |
| Change MAPQ threshold                         | `params.mapq`                              |
| Add a new mitochondrial chrom name            | `params.mito_chroms`                       |
| Switch from CPM to RPGC normalization         | `params.bigwig_norm` + `params.effective_genome_size` |
| Add a new aligner                             | new module + swap `BOWTIE2_ALIGN` in `main.nf` |
| Add a new peak caller (e.g. SEACR)            | new module + branch in stage 6/7           |
| Pin a different container image               | `conf/{docker,singularity}.config`         |
| Change how narrow/broad is inferred           | `modules/macs3.nf::infer_peak_mode`        |
| Add a new validation rule                     | `bin/validate_associations.py`             |
| Add a new MultiQC custom-content section      | emit `*_mqc.{yaml,tsv}` from a process; pass it into `MULTIQC_FINAL`'s `custom_content` |
