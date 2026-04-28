#!/usr/bin/env bash
# ----------------------------------------------------------------------------
# make_synth_test_data.sh
# Generate a minimal synthetic test dataset for the cuttag-dsl2 pipeline.
#
# Strategy:
#   1. Subsample 10 000 paired reads from a real CUT&Tag FASTQ pair using seqtk.
#   2. Optionally restrict to chr22 by aligning + extracting reads.
#   3. Build a chr22-only Bowtie2 index from a real reference FASTA.
#   4. Drop the resulting test files into test/fastq/ and test/refs/ so that
#      `nextflow run main.nf -profile test,singularity ...` works.
#
# Run from the pipeline root, e.g.:
#   bash test/make_synth_test_data.sh \
#       --r1   /workdir/$USER/seed/SAMPLE_R1.fastq.gz \
#       --r2   /workdir/$USER/seed/SAMPLE_R2.fastq.gz \
#       --hg38 /local/storage/data/short_read_index/hg38/genome.fa
# ----------------------------------------------------------------------------
set -euo pipefail

R1=""
R2=""
HG38_FA=""
N_READS=10000
SEED=42

usage() {
    cat <<EOF
Usage: $0 --r1 <R1.fastq.gz> --r2 <R2.fastq.gz> --hg38 <genome.fa> [--n 10000]
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --r1)   R1="$2"; shift 2 ;;
        --r2)   R2="$2"; shift 2 ;;
        --hg38) HG38_FA="$2"; shift 2 ;;
        --n)    N_READS="$2"; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "Unknown option: $1" >&2; usage; exit 1 ;;
    esac
done

if [[ -z "$R1" || -z "$R2" || -z "$HG38_FA" ]]; then
    usage
    exit 1
fi

for tool in seqtk samtools bowtie2-build; do
    if ! command -v $tool >/dev/null 2>&1; then
        echo "ERROR: $tool not in PATH" >&2
        exit 1
    fi
done

PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT_FASTQ="${PIPELINE_DIR}/test/fastq"
OUT_REFS="${PIPELINE_DIR}/test/refs"
mkdir -p "${OUT_FASTQ}" "${OUT_REFS}"

# Read the test_associations.csv to get sample names
SAMPLES=$(awk -F, 'NR > 1 { print $1 }' "${PIPELINE_DIR}/test/test_associations.csv")

echo "## 1/3  Subsampling $N_READS read pairs per sample (seed=$SEED)"
for sid in $SAMPLES; do
    seqtk sample -s${SEED} "$R1" "$N_READS" | gzip > "${OUT_FASTQ}/${sid}_R1.fastq.gz"
    seqtk sample -s${SEED} "$R2" "$N_READS" | gzip > "${OUT_FASTQ}/${sid}_R2.fastq.gz"
    echo "    $sid"
    SEED=$((SEED + 1))     # different seed per sample → not identical files
done

echo "## 2/3  Building chr22-only Bowtie2 index"
samtools faidx "$HG38_FA" chr22 > "${OUT_REFS}/chr22.fa"
samtools faidx "${OUT_REFS}/chr22.fa"
cut -f1,2 "${OUT_REFS}/chr22.fa.fai" > "${OUT_REFS}/chr22.chrom.sizes"
bowtie2-build --threads 4 "${OUT_REFS}/chr22.fa" "${OUT_REFS}/chr22"

echo "## 3/3  Mini GTF + TSS BED (chr22 only)"
cat > "${OUT_REFS}/chr22.gtf" <<'EOF'
chr22	test	gene	16000000	16001000	.	+	.	gene_id "g1"; transcript_id "t1";
chr22	test	transcript	16000000	16001000	.	+	.	gene_id "g1"; transcript_id "t1"; transcript_biotype "protein_coding";
chr22	test	exon	16000000	16001000	.	+	.	gene_id "g1"; transcript_id "t1";
chr22	test	gene	17500000	17501500	.	-	.	gene_id "g2"; transcript_id "t2";
chr22	test	transcript	17500000	17501500	.	-	.	gene_id "g2"; transcript_id "t2"; transcript_biotype "protein_coding";
chr22	test	exon	17500000	17501500	.	-	.	gene_id "g2"; transcript_id "t2";
EOF
"${PIPELINE_DIR}/bin/make_tss_bed.py" \
    --gtf "${OUT_REFS}/chr22.gtf" \
    --out "${OUT_REFS}/chr22_tss.bed" \
    --protein_coding

echo
echo "Done. Run the test profile with:"
echo "  cd ${PIPELINE_DIR}"
echo "  module load nextflow/25.4.3"
echo "  nextflow run main.nf -profile test,singularity \\"
echo "      --bowtie2_index   \$PWD/test/refs/chr22 \\"
echo "      --chrom_sizes     \$PWD/test/refs/chr22.chrom.sizes \\"
echo "      --annotation_gtf  \$PWD/test/refs/chr22.gtf \\"
echo "      --tss_bed         \$PWD/test/refs/chr22_tss.bed"
