#!/usr/bin/env bash
# ----------------------------------------------------------------------------
# run_cuttag_pipeline.sh
#
# Launch script for the cuttag-dsl2 Nextflow pipeline on Cornell BioHPC.
#
# REVIEW EVERY PATH BELOW BEFORE EXECUTING — placeholders are marked TODO.
# Run this script in a separate terminal session (not via Claude/automation),
# inside e.g. `screen` or `tmux`, since the pipeline can take many hours.
# ----------------------------------------------------------------------------

set -euo pipefail

# 1. Load Nextflow (BioHPC module system)
module load nextflow/25.4.3 || {
    echo "ERROR: failed to module load nextflow/25.4.3" >&2
    exit 1
}

# 2. Workdir / scratch / cache — keep everything under /workdir/$USER
export NXF_WORK="/workdir/${USER}/cuttag_run/work"
export NXF_SINGULARITY_CACHEDIR="/workdir/${USER}/.singularity_cache"
export TMPDIR="/workdir/${USER}/tmp"
mkdir -p "${NXF_WORK}" "${NXF_SINGULARITY_CACHEDIR}" "${TMPDIR}"

# 3. Pipeline location and outputs
PIPELINE_DIR="/workdir/${USER}/cuttag_pipeline"   # TODO: confirm
OUTDIR="/workdir/${USER}/cuttag_run/results"
mkdir -p "${OUTDIR}"

# 4. Inputs — EDIT THESE
INPUT_DIR="/workdir/${USER}/cuttag_run/fastq"               # TODO: directory of raw FASTQs
ASSOC_CSV="/workdir/${USER}/cuttag_run/samplesheet.csv"     # TODO: your association CSV
GENOME="hg38"                                                # TODO: hg38 | mm10 | rn6 | custom
BT2_INDEX="/workdir/${USER}/refs/hg38/bowtie2/hg38"          # TODO: bowtie2 index basename
CHROM_SIZES="/workdir/${USER}/refs/hg38/hg38.chrom.sizes"    # TODO
GTF="/workdir/${USER}/refs/hg38/hg38.gtf"                    # TODO
TSS_BED="/workdir/${USER}/refs/hg38/hg38_tss.bed"            # TODO
BLACKLIST="/workdir/${USER}/refs/hg38/hg38-blacklist.v2.bed" # TODO (omit --blacklist if not used)

# 5. Optional knobs
GENOME_SIZE="hs"   # hs | mm | rn | dm | ce | <int>
THREADS=8

cd "${PIPELINE_DIR}"

# 6. Launch
nextflow run main.nf \
    --input_dir       "${INPUT_DIR}" \
    --outdir          "${OUTDIR}" \
    --genome          "${GENOME}" \
    --bowtie2_index   "${BT2_INDEX}" \
    --chrom_sizes     "${CHROM_SIZES}" \
    --annotation_gtf  "${GTF}" \
    --tss_bed         "${TSS_BED}" \
    --blacklist       "${BLACKLIST}" \
    --association_csv "${ASSOC_CSV}" \
    --genome_size     "${GENOME_SIZE}" \
    --threads         "${THREADS}" \
    -profile          singularity \
    -work-dir         "${NXF_WORK}" \
    -resume

echo "Done. Reports in ${OUTDIR}/13_multiqc/"
