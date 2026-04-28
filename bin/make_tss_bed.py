#!/usr/bin/env python3
"""
make_tss_bed.py
Extract a TSS BED file from a GTF/GFF annotation. The pipeline requires --tss_bed
for the deepTools TSS-anchored matrix; this helper builds one in the format
deepTools expects (BED6 with strand).

For each transcript, emits a single 1-bp BED record at the TSS:
    +-strand transcript: chrom  start    start+1  transcript_id  .  +
    --strand transcript: chrom  end-1    end      transcript_id  .  -

Optionally collapses to the protein-coding subset, or to one TSS per gene.

Usage:
    make_tss_bed.py --gtf hg38.gtf --out hg38_tss.bed
    make_tss_bed.py --gtf hg38.gtf --out hg38_tss_protein_coding.bed --protein_coding
    make_tss_bed.py --gtf hg38.gtf --out hg38_tss_one_per_gene.bed --one_per_gene
"""
from __future__ import annotations

import argparse
import gzip
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Iterator, Tuple


def open_text(p: Path):
    if str(p).endswith(".gz"):
        return gzip.open(p, "rt")
    return open(p, "rt")


ATTR_RE = re.compile(r'(\w+)\s+"([^"]*)"')
ATTR_RE_GFF3 = re.compile(r"(\w+)=([^;]+)")


def parse_attrs(field: str) -> dict:
    """Parse the 9th GTF/GFF column into a dict."""
    out = dict(ATTR_RE.findall(field))
    if not out:
        out = dict(ATTR_RE_GFF3.findall(field))
    return out


def iter_transcripts(gtf: Path) -> Iterator[Tuple[str, int, int, str, str, str, str]]:
    """Yield (chrom, start, end, transcript_id, gene_id, strand, biotype) for each transcript line."""
    with open_text(gtf) as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            if f[2] != "transcript":
                continue
            chrom  = f[0]
            start  = int(f[3]) - 1   # GTF is 1-based, BED is 0-based
            end    = int(f[4])
            strand = f[6]
            attrs  = parse_attrs(f[8])
            tid = attrs.get("transcript_id") or attrs.get("ID") or ""
            gid = attrs.get("gene_id")        or attrs.get("Parent") or ""
            biotype = (attrs.get("transcript_biotype")
                       or attrs.get("transcript_type")
                       or attrs.get("gene_biotype")
                       or attrs.get("gene_type")
                       or "")
            if not tid:
                continue
            yield chrom, start, end, tid, gid, strand, biotype


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--gtf", required=True, help="Input GTF/GFF (optionally gzipped)")
    p.add_argument("--out", required=True, help="Output BED file (tab-delimited)")
    p.add_argument("--protein_coding", action="store_true",
                   help="Restrict to transcripts whose biotype is 'protein_coding'")
    p.add_argument("--one_per_gene", action="store_true",
                   help="Keep only one TSS per gene (the most upstream one)")
    args = p.parse_args()

    gtf = Path(args.gtf)
    if not gtf.exists():
        print(f"ERROR: GTF not found: {gtf}", file=sys.stderr)
        return 1

    records = []
    n_total = 0
    for chrom, start, end, tid, gid, strand, biotype in iter_transcripts(gtf):
        n_total += 1
        if args.protein_coding and biotype != "protein_coding":
            continue
        if strand == "+":
            tss_start = start
            tss_end   = start + 1
        elif strand == "-":
            tss_start = end - 1
            tss_end   = end
        else:
            continue  # skip features without strand
        records.append((chrom, tss_start, tss_end, tid, gid, strand))

    if args.one_per_gene:
        # Keep the most-upstream TSS per gene_id (smallest coordinate on +, largest on -)
        by_gene: dict = defaultdict(list)
        for r in records:
            chrom, tss_start, tss_end, tid, gid, strand = r
            by_gene[gid or tid].append(r)
        kept = []
        for gid, rows in by_gene.items():
            plus  = [r for r in rows if r[5] == "+"]
            minus = [r for r in rows if r[5] == "-"]
            if plus:
                kept.append(min(plus, key=lambda r: r[1]))
            if minus:
                kept.append(max(minus, key=lambda r: r[1]))
        records = kept

    # Sort by chrom, start
    records.sort(key=lambda r: (r[0], r[1]))

    with open(args.out, "w") as out:
        for chrom, s, e, tid, gid, strand in records:
            # BED6: chrom start end name score strand
            out.write(f"{chrom}\t{s}\t{e}\t{tid}\t.\t{strand}\n")

    print(f"OK: parsed {n_total} transcripts from GTF; wrote {len(records)} TSS records to {args.out}",
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
