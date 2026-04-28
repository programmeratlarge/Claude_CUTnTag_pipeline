# Preparing reference files

The pipeline requires a Bowtie2 index, chromosome sizes, GTF, TSS BED, and
optionally an ENCODE blacklist. This document shows how to assemble those on
**Cornell BioHPC** — paths under `/local/storage/data/` are the curated mirrors
maintained by BioHPC; you are not required to re-download anything.

> **Storage rule on BioHPC**: writable scratch and outputs must live under
> `/workdir/$USER/`. Do **not** stage references under `/home/$USER/` — that
> filesystem has tight quotas and is not meant for large data.

---

## hg38 (human)

Most of the human references are pre-built on BioHPC under
`/local/storage/data/short_read_index/hg38/`. Verify the exact path on your
compute node, then copy or symlink:

```bash
mkdir -p /workdir/$USER/refs/hg38
cd       /workdir/$USER/refs/hg38

# Bowtie2 index — basename is whatever the existing files share (often 'hg38' or 'genome')
ln -s /local/storage/data/short_read_index/hg38/bowtie2/genome.* .
BT2_BASENAME=$PWD/genome   # adjust to match the symlinked basename

# Chromosome sizes
ln -s /local/storage/data/short_read_index/hg38/hg38.chrom.sizes .

# GTF (Gencode is the usual choice)
ln -s /local/storage/data/annotation/hg38/gencode.v45.annotation.gtf hg38.gtf
```

If the BioHPC mirror name differs on your node, search with:

```bash
find /local/storage/data -maxdepth 4 -name 'hg38*' 2>/dev/null
```

### TSS BED — derive from the GTF

The pipeline ships a helper:

```bash
$PIPELINE_DIR/bin/make_tss_bed.py \
    --gtf hg38.gtf \
    --out hg38_tss.bed \
    --protein_coding \
    --one_per_gene
```

`--protein_coding` keeps roughly 20 k transcripts; `--one_per_gene` further
collapses to one TSS per gene (recommended for cleaner deepTools profiles).

### ENCODE blacklist (optional but recommended)

```bash
wget -O hg38-blacklist.v2.bed.gz \
    https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg38-blacklist.v2.bed.gz
gunzip hg38-blacklist.v2.bed.gz
```

### Verify

```bash
ls -la /workdir/$USER/refs/hg38/
# Expected: bowtie2 index files, hg38.chrom.sizes, hg38.gtf, hg38_tss.bed, hg38-blacklist.v2.bed
```

---

## mm10 (mouse)

Same layout, different paths:

```bash
mkdir -p /workdir/$USER/refs/mm10 && cd /workdir/$USER/refs/mm10

ln -s /local/storage/data/short_read_index/mm10/bowtie2/genome.* .
ln -s /local/storage/data/short_read_index/mm10/mm10.chrom.sizes .
ln -s /local/storage/data/annotation/mm10/gencode.vM25.annotation.gtf mm10.gtf

$PIPELINE_DIR/bin/make_tss_bed.py --gtf mm10.gtf --out mm10_tss.bed --protein_coding --one_per_gene

wget -O mm10-blacklist.v2.bed.gz \
    https://github.com/Boyle-Lab/Blacklist/raw/master/lists/mm10-blacklist.v2.bed.gz
gunzip mm10-blacklist.v2.bed.gz
```

When running for mm10, also pass:
```
--genome_size           mm
--effective_genome_size 2652783500
```

---

## Building a Bowtie2 index from scratch

Only needed if BioHPC doesn't already have an index for your genome:

```bash
cd /workdir/$USER/refs/myorganism
bowtie2-build --threads 8 genome.fa genome
# now BT2_BASENAME=/workdir/$USER/refs/myorganism/genome
```

Use 8–16 threads; for hg38 this takes ~1.5 hours and ~5 GB RAM.

---

## Getting `chrom.sizes` from a FASTA

```bash
samtools faidx genome.fa
cut -f1,2 genome.fa.fai > genome.chrom.sizes
```

---

## Validating your reference set

Before launching, sanity-check that everything resolves:

```bash
# Bowtie2 index files exist?
ls $BT2_BASENAME.*.bt2 || echo "missing index"

# chrom.sizes is non-empty and tab-delimited?
head -3 hg38.chrom.sizes

# GTF parses?
head -3 hg38.gtf

# TSS BED has reasonable content?
wc -l hg38_tss.bed         # ~20-60k typical
head -3 hg38_tss.bed
```

If any of those look wrong, the pipeline will exit at the parameter-validation
stage with a pointed error message rather than crashing mid-run.

---

## Effective genome size cheat sheet

For `--bigwig_norm RPGC`, set `--effective_genome_size` to the unique-mapping
genome size. From the deepTools documentation table:

| Genome | `--effective_genome_size` | `--genome_size` (MACS3) |
|--------|--------------------------:|:------------------------|
| hg19   |              2,864,785,220 | hs                      |
| hg38   |              2,913,022,398 | hs                      |
| mm10   |              2,652,783,500 | mm                      |
| mm39   |              2,654,621,837 | mm                      |
| dm6    |                142,573,017 | dm                      |
| ce11   |                 95,159,452 | ce                      |
