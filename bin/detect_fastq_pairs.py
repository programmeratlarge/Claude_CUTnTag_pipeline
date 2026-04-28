#!/usr/bin/env python3
"""
detect_fastq_pairs.py
Scan a directory for paired-end FASTQ files, validate every R1 has a matching R2,
and emit a TSV manifest plus a log file.

Recognized R1/R2 patterns (all case-insensitive):
    *_R1*.fastq.gz / *_R2*.fastq.gz
    *_R1*.fq.gz    / *_R2*.fq.gz
    *_1.fastq.gz   / *_2.fastq.gz       (also .fq.gz)
    *_1.fq.gz      / *_2.fq.gz
    *.R1.fastq.gz  / *.R2.fastq.gz      (dot-separated)
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path
from typing import Dict, List, Tuple

EXTS = (".fastq.gz", ".fq.gz", ".fastq", ".fq")

R1_RE = [
    re.compile(r"^(?P<base>.+?)_R1(?P<suffix>(_\d+)?)(?P<ext>\.f(ast)?q(\.gz)?)$", re.I),
    re.compile(r"^(?P<base>.+?)_1(?P<ext>\.f(ast)?q(\.gz)?)$", re.I),
    re.compile(r"^(?P<base>.+?)\.R1(?P<ext>\.f(ast)?q(\.gz)?)$", re.I),
]
R2_TEMPLATES = [
    ("R1", "R2"),  # the first capture-group naming
    ("_1.", "_2."),
    (".R1.", ".R2."),
]


def find_mate(r1_path: Path) -> Path | None:
    """Try to locate the matching R2 file for a given R1 file."""
    name = r1_path.name
    candidates = []
    # Substitution-based: replace the first occurrence of R1/_1./.R1. with R2/_2./.R2.
    for old, new in R2_TEMPLATES:
        if old in name:
            cand = r1_path.with_name(name.replace(old, new, 1))
            if cand.exists():
                candidates.append(cand)
    # De-duplicate while preserving order
    seen = set()
    deduped = []
    for c in candidates:
        if c not in seen:
            seen.add(c)
            deduped.append(c)
    if len(deduped) == 1:
        return deduped[0]
    return None


def derive_sample_id(r1_name: str) -> str:
    """Strip lane/_R1/_001/ext suffix to derive a sample ID."""
    s = r1_name
    for ext in EXTS:
        if s.lower().endswith(ext):
            s = s[: -len(ext)]
            break
    s = re.sub(r"_R1(_\d+)?$", "", s, flags=re.I)
    s = re.sub(r"_1$", "", s)
    s = re.sub(r"\.R1$", "", s, flags=re.I)
    s = re.sub(r"_S\d+_L\d{3}$", "", s)  # _S1_L001
    s = re.sub(r"_L\d{3}$", "", s)
    return s


def is_r1(name: str) -> bool:
    for rx in R1_RE:
        if rx.match(name):
            return True
    return False


def is_r2(name: str) -> bool:
    """Heuristic: skip files that are clearly the R2 mate of a pair."""
    return (
        re.match(r"^.+?_R2(_\d+)?\.f(ast)?q(\.gz)?$", name, re.I) is not None
        or re.match(r"^.+?_2\.f(ast)?q(\.gz)?$", name, re.I) is not None
        or re.match(r"^.+?\.R2\.f(ast)?q(\.gz)?$", name, re.I) is not None
    )


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input_dir", required=True)
    p.add_argument("--pattern", default="auto",
                   help="auto, or a regex with named groups 'base' and 'mate'")
    p.add_argument("--out_tsv", required=True)
    p.add_argument("--log", required=True)
    args = p.parse_args()

    in_dir = Path(args.input_dir).resolve()
    if not in_dir.is_dir():
        print(f"ERROR: input_dir not a directory: {in_dir}", file=sys.stderr)
        return 1

    log_lines: List[str] = []
    log_lines.append(f"# detect_fastq_pairs.py")
    log_lines.append(f"# input_dir = {in_dir}")
    log_lines.append(f"# pattern   = {args.pattern}")

    all_files = sorted(p for p in in_dir.iterdir() if p.is_file()
                       and any(p.name.lower().endswith(e) for e in EXTS))
    log_lines.append(f"# total FASTQ candidates: {len(all_files)}")

    pairs: List[Tuple[str, Path, Path, str]] = []
    seen_r2: set = set()
    skipped: List[str] = []
    warnings: List[str] = []

    for f in all_files:
        if is_r2(f.name):
            continue
        if not is_r1(f.name):
            skipped.append(f"SKIP (no R1 pattern match): {f.name}")
            continue
        mate = find_mate(f)
        if mate is None:
            warnings.append(f"WARN: no R2 mate found for {f.name}")
            continue
        sample_id = derive_sample_id(f.name)
        pairs.append((sample_id, f, mate, "auto"))
        seen_r2.add(mate.name)

    # Detect ambiguous: R2 files not paired
    for f in all_files:
        if is_r2(f.name) and f.name not in seen_r2:
            warnings.append(f"WARN: orphan R2 file (no R1 partner): {f.name}")

    # Detect duplicate sample_ids
    sample_counts: Dict[str, int] = {}
    for s, _, _, _ in pairs:
        sample_counts[s] = sample_counts.get(s, 0) + 1
    duplicates = [k for k, v in sample_counts.items() if v > 1]
    for d in duplicates:
        warnings.append(f"ERROR: duplicate sample_id derived: {d}")

    # Write manifest
    with open(args.out_tsv, "w") as out:
        out.write("sample_id\tr1\tr2\tinput_dir\tpattern\tstatus\n")
        for sid, r1, r2, pat in pairs:
            status = "ERROR_DUPLICATE" if sid in duplicates else "OK"
            out.write(f"{sid}\t{r1}\t{r2}\t{in_dir}\t{pat}\t{status}\n")

    # Write log
    with open(args.log, "w") as logf:
        for ln in log_lines:
            logf.write(ln + "\n")
        logf.write(f"# pairs_detected: {len(pairs)}\n")
        for sid, r1, r2, pat in pairs:
            logf.write(f"PAIR\t{sid}\t{r1.name}\t{r2.name}\n")
        for s in skipped:
            logf.write(s + "\n")
        for w in warnings:
            logf.write(w + "\n")

    if duplicates:
        print(f"ERROR: duplicate sample_ids: {duplicates}", file=sys.stderr)
        return 2
    if not pairs:
        print(f"ERROR: no valid paired-end FASTQ files found in {in_dir}", file=sys.stderr)
        return 3

    print(f"OK: {len(pairs)} pairs detected. See {args.out_tsv} and {args.log}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
