#!/usr/bin/env python3
"""
read_retention_summary.py
Combine per-stage read counts for one sample into a TSV and a MultiQC custom-content
YAML so the final report can show how many reads survived each filtering step.
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path


def parse_cutadapt_log(p: Path) -> dict:
    """Parse a cutadapt log for paired-end stats. Returns counts in read PAIRS."""
    txt = p.read_text(errors="ignore")
    out = {"input_pairs": 0, "passing_pairs": 0}
    m = re.search(r"Total read pairs processed:\s+([\d,]+)", txt)
    if m:
        out["input_pairs"] = int(m.group(1).replace(",", ""))
    m = re.search(r"Pairs written \(passing filters\):\s+([\d,]+)", txt)
    if m:
        out["passing_pairs"] = int(m.group(1).replace(",", ""))
    return out


def parse_bowtie2_log(p: Path) -> dict:
    """Parse Bowtie2's stderr summary for paired-end alignment stats."""
    txt = p.read_text(errors="ignore")
    out = {"reads_processed_pairs": 0, "aligned_at_least_once_pairs": 0}
    m = re.search(r"^(\d+)\s+reads;\s+of these:", txt, flags=re.M)
    if m:
        # Bowtie2 reports total reads where each pair == 1 read entry, so this is pairs already.
        out["reads_processed_pairs"] = int(m.group(1))
    aligned = 0
    for label in [r"aligned concordantly exactly 1 time",
                  r"aligned concordantly >1 times",
                  r"aligned discordantly 1 time"]:
        m = re.search(rf"\s+(\d+)\s+\(.+?\) {label}", txt)
        if m:
            aligned += int(m.group(1))
    out["aligned_at_least_once_pairs"] = aligned
    return out


def parse_filter_counts(p: Path) -> dict:
    """Filter counts TSV produced by FILTER_MITO_BLACKLIST. Counts are in read records."""
    out = {"aligned_reads": 0, "mapq_pp_reads": 0, "nomito_reads": 0, "noblacklist_reads": 0}
    with open(p) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            row = dict(zip(header, line.rstrip("\n").split("\t")))
            out["aligned_reads"]    = int(row.get("aligned", 0))
            out["mapq_pp_reads"]    = int(row.get("mapq_pp", 0))
            out["nomito_reads"]     = int(row.get("nomito", 0))
            out["noblacklist_reads"] = int(row.get("noblacklist", 0))
    return out


def parse_dedup(p: Path) -> dict:
    """
    Parse Picard MarkDuplicates metrics OR umi_tools dedup log for the
    'after dedup' read count.
    """
    txt = p.read_text(errors="ignore")
    # Picard
    m = re.search(r"^LIBRARY\s+UNPAIRED_READS_EXAMINED\s+READ_PAIRS_EXAMINED.*$",
                  txt, flags=re.M)
    if m:
        # Find the data line: header line then next non-empty line
        lines = txt.splitlines()
        for i, ln in enumerate(lines):
            if ln.startswith("LIBRARY") and "READ_PAIRS_EXAMINED" in ln:
                hdr = ln.split("\t")
                if i + 1 < len(lines):
                    data = lines[i + 1].split("\t")
                    rec = dict(zip(hdr, data))
                    pairs_examined = int(rec.get("READ_PAIRS_EXAMINED", 0))
                    pair_dups      = int(rec.get("READ_PAIR_DUPLICATES", 0))
                    return {
                        "after_dedup_reads": (pairs_examined - pair_dups) * 2,
                        "duplicate_reads":   pair_dups * 2,
                    }
    # umi_tools
    m = re.search(r"Number of reads out:\s+(\d+)", txt)
    if m:
        return {"after_dedup_reads": int(m.group(1)), "duplicate_reads": 0}
    return {"after_dedup_reads": 0, "duplicate_reads": 0}


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--sample", required=True)
    p.add_argument("--cutadapt", required=True)
    p.add_argument("--bowtie2",  required=True)
    p.add_argument("--filter",   required=True)
    p.add_argument("--dedup",    required=True)
    p.add_argument("--out_tsv",  required=True)
    p.add_argument("--out_mqc",  required=True)
    args = p.parse_args()

    cut    = parse_cutadapt_log(Path(args.cutadapt))
    bt2    = parse_bowtie2_log(Path(args.bowtie2))
    filt   = parse_filter_counts(Path(args.filter))
    dedup  = parse_dedup(Path(args.dedup))

    # Convert pair counts to read counts where downstream is in reads
    raw_reads          = cut["input_pairs"] * 2
    trimmed_reads      = cut["passing_pairs"] * 2
    aligned_reads      = filt["aligned_reads"]
    mapq_pp_reads      = filt["mapq_pp_reads"]
    nomito_reads       = filt["nomito_reads"]
    noblacklist_reads  = filt["noblacklist_reads"]
    after_dedup_reads  = dedup["after_dedup_reads"]

    rows = [
        ("raw_reads",         raw_reads),
        ("trimmed_reads",     trimmed_reads),
        ("aligned_reads",     aligned_reads),
        ("properly_paired",   mapq_pp_reads),
        ("non_mito_reads",    nomito_reads),
        ("non_blacklist",     noblacklist_reads),
        ("final_after_dedup", after_dedup_reads),
    ]

    with open(args.out_tsv, "w") as out:
        out.write("sample\t" + "\t".join(k for k, _ in rows) + "\n")
        out.write(args.sample + "\t" + "\t".join(str(v) for _, v in rows) + "\n")

    # MultiQC custom content (YAML format, plot type bargraph)
    mqc = f"""id: 'read_retention_{args.sample}'
section_name: 'Read retention through pipeline'
description: 'Number of reads remaining after each filtering step for sample {args.sample}.'
plot_type: 'bargraph'
pconfig:
    id: 'read_retention_bargraph_{args.sample}'
    title: 'Read retention: {args.sample}'
    ylab: 'Reads'
    cpswitch_counts_label: 'Reads'
data:
    {args.sample}:
"""
    for k, v in rows:
        mqc += f"        {k}: {v}\n"
    Path(args.out_mqc).write_text(mqc)
    return 0


if __name__ == "__main__":
    sys.exit(main())
