#!/usr/bin/env python3
"""
compute_frip.py
Compute the Fraction of Reads in Peaks (FRiP) for a BAM against a peaks BED/narrowPeak/broadPeak file.

Strategy: count total mapped reads in BAM and reads-in-peaks via samtools view + bedtools intersect.
Emits both a TSV row and a MultiQC custom-content YAML.
"""
from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path


def run(cmd: list[str]) -> str:
    res = subprocess.run(cmd, check=True, capture_output=True, text=True)
    return res.stdout


def count_total_reads(bam: str) -> int:
    out = run(["samtools", "view", "-c", "-F", "0x4", bam])
    return int(out.strip())


def count_reads_in_peaks(bam: str, peaks: str) -> int:
    """Use samtools + bedtools to count reads overlapping peak regions."""
    # samtools view -L peaks.bed -c excludes secondary/supplementary by default? No - we add -F 0x100 -F 0x800
    out = run(["samtools", "view", "-c", "-F", "0x4", "-F", "0x100", "-F", "0x800", "-L", peaks, bam])
    return int(out.strip())


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--sample", required=True)
    p.add_argument("--bam", required=True)
    p.add_argument("--peaks", required=True)
    p.add_argument("--label", required=True, help="own | group | consensus")
    p.add_argument("--out_tsv", required=True)
    p.add_argument("--out_mqc", required=True)
    args = p.parse_args()

    if not shutil.which("samtools"):
        print("ERROR: samtools not found", file=sys.stderr)
        return 1

    peaks_path = Path(args.peaks)
    if peaks_path.stat().st_size == 0:
        total = count_total_reads(args.bam)
        in_peaks = 0
        frip = 0.0
    else:
        total = count_total_reads(args.bam)
        in_peaks = count_reads_in_peaks(args.bam, args.peaks)
        frip = (in_peaks / total) if total > 0 else 0.0

    with open(args.out_tsv, "w") as out:
        out.write("sample\tlabel\ttotal_reads\treads_in_peaks\tFRiP\n")
        out.write(f"{args.sample}\t{args.label}\t{total}\t{in_peaks}\t{frip:.6f}\n")

    mqc = f"""id: 'frip_{args.sample}_{args.label}'
section_name: 'FRiP scores'
description: 'Fraction of Reads in Peaks for {args.sample} vs {args.label} peaks.'
plot_type: 'bargraph'
pconfig:
    id: 'frip_{args.sample}_{args.label}_bg'
    title: 'FRiP: {args.sample} ({args.label})'
    ylab: 'Reads'
data:
    '{args.sample} ({args.label})':
        reads_in_peaks: {in_peaks}
        reads_outside_peaks: {max(total - in_peaks, 0)}
"""
    Path(args.out_mqc).write_text(mqc)
    return 0


if __name__ == "__main__":
    sys.exit(main())
