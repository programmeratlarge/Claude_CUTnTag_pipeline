#!/usr/bin/env python3
"""
detect_fastq_pairs.py
Recursively scan a directory tree for paired-end FASTQ files, validate every R1
has a matching R2, and emit a TSV manifest plus a log file.

The search descends into ALL subdirectories under --input_dir. Hidden
directories (those whose name starts with '.') are skipped. R2 mates are
expected in the same directory as their R1 (standard sequencer layout); the
script does NOT cross-link an R1 in one subdir to an R2 in another.

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
from typing import Dict, Iterator, List, Tuple

EXTS = (".fastq.gz", ".fq.gz", ".fastq", ".fq")

R1_RE = [
    re.compile(r"^(?P<base>.+?)_R1(?P<suffix>(_\d+)?)(?P<ext>\.f(ast)?q(\.gz)?)$", re.I),
    re.compile(r"^(?P<base>.+?)_1(?P<ext>\.f(ast)?q(\.gz)?)$", re.I),
    re.compile(r"^(?P<base>.+?)\.R1(?P<ext>\.f(ast)?q(\.gz)?)$", re.I),
]


R2_PATTERNS: List[Tuple[re.Pattern, str]] = [
    # Each pattern matches the read-direction marker at the END of the filename
    # (possibly followed by Illumina's _NNN block, e.g. _R1_001.fastq.gz).
    # The first capture group is the part to replace; the second is the
    # suffix that must be preserved.
    (re.compile(r"_R1(_\d+)?(\.f(?:ast)?q(?:\.gz)?)$", re.I), "_R2"),
    (re.compile(r"_1(\.f(?:ast)?q(?:\.gz)?)$",         re.I), "_2"),
    (re.compile(r"\.R1(\.f(?:ast)?q(?:\.gz)?)$",       re.I), ".R2"),
]


def find_mate(r1_path: Path) -> Path | None:
    """Locate the R2 mate of an R1 file by anchoring the R1→R2 substitution at
    the END of the filename.

    This avoids miscounting sample names that themselves contain '_R1_' as a
    replicate marker (e.g. ``H3K27ac_WT_R1_S1_L001_R1_001.fastq.gz`` — only
    the LAST ``_R1`` is the read-direction marker).
    """
    name = r1_path.name
    candidates: List[Path] = []
    for rx, replacement in R2_PATTERNS:
        m = rx.search(name)
        if not m:
            continue
        # Splice: name[:m.start()] + replacement + (groups after the marker)
        # group(1) is the optional _NNN block (or None); group(-1) is always the extension
        ext = m.group(m.lastindex)
        mid = m.group(1) if m.lastindex >= 2 else ""
        mid = mid or ""
        cand_name = name[:m.start()] + replacement + mid + ext
        cand = r1_path.with_name(cand_name)
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


def iter_fastq_files(root: Path) -> Iterator[Path]:
    """Yield FASTQ files anywhere under ``root``, skipping hidden directories.

    A directory is "hidden" if its name (any component of the path relative to
    root) starts with '.'. The root itself is not subject to this rule.
    """
    for p in root.rglob("*"):
        rel_parts = p.relative_to(root).parts
        # Skip if any ancestor directory (not the file itself) is hidden
        if any(part.startswith(".") for part in rel_parts[:-1]):
            continue
        if not p.is_file():
            continue
        if not any(p.name.lower().endswith(e) for e in EXTS):
            continue
        yield p


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
    log_lines.append(f"# search    = recursive (all subdirectories; hidden dirs skipped)")

    all_files = sorted(iter_fastq_files(in_dir))
    # Tally distinct subdirectories visited (for diagnostic clarity)
    subdirs_seen = sorted({str(f.parent.relative_to(in_dir)) or "." for f in all_files})
    log_lines.append(f"# subdirectories with FASTQs: {len(subdirs_seen)}")
    for sd in subdirs_seen:
        log_lines.append(f"#   {sd}")
    log_lines.append(f"# total FASTQ candidates: {len(all_files)}")

    pairs: List[Tuple[str, Path, Path, str]] = []
    seen_r2: set = set()        # set of resolved Path objects (NOT just names)
    skipped: List[str] = []
    warnings: List[str] = []

    def rel(p: Path) -> str:
        try:
            return str(p.relative_to(in_dir))
        except ValueError:
            return str(p)

    for f in all_files:
        if is_r2(f.name):
            continue
        if not is_r1(f.name):
            skipped.append(f"SKIP (no R1 pattern match): {rel(f)}")
            continue
        mate = find_mate(f)
        if mate is None:
            warnings.append(f"WARN: no R2 mate found for {rel(f)}")
            continue
        sample_id = derive_sample_id(f.name)
        pairs.append((sample_id, f, mate, "auto"))
        seen_r2.add(mate.resolve())

    # Detect ambiguous: R2 files not paired (by full path, since names can repeat across subdirs)
    for f in all_files:
        if is_r2(f.name) and f.resolve() not in seen_r2:
            warnings.append(f"WARN: orphan R2 file (no R1 partner): {rel(f)}")

    # Detect duplicate sample_ids — and list every file involved so the user can fix the layout
    by_sample: Dict[str, List[Path]] = {}
    for s, r1, _, _ in pairs:
        by_sample.setdefault(s, []).append(r1)
    duplicates = [k for k, v in by_sample.items() if len(v) > 1]
    for d in duplicates:
        files_listed = ", ".join(rel(p) for p in by_sample[d])
        warnings.append(f"ERROR: duplicate sample_id '{d}' derived from: {files_listed}")

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
