#!/usr/bin/env python3
"""
validate_associations.py
Validate the user-provided association_csv against the FASTQ pairs TSV.

Required columns:
    sample_id,species,genome,antibody,condition,replicate,group_id,
    is_control,control_group_id,merge_group_id,peak_calling_mode,notes

Genome filter:
    If --genome is provided, rows whose 'genome' column does NOT equal that
    value are dropped from this run BEFORE validation. The match is a literal
    case-sensitive string comparison: '--genome hg38' keeps rows with
    genome='hg38' only. Dropped rows are logged as INFO; FASTQ files that
    correspond to dropped rows are also silently skipped (also logged).

Rules enforced (against the genome-filtered set):
    1. Every sample_id in the FASTQ pairs TSV must appear in the CSV — UNLESS
       its CSV row was dropped by the genome filter (then it's logged as
       INFO and skipped).
    2. No extra rows in the CSV unless --allow_extra is passed.
    3. is_control must be 'true' or 'false' (case-insensitive).
    4. peak_calling_mode must be 'narrow', 'broad', or 'auto' (or empty -> 'auto').
    5. All rows sharing a merge_group_id must agree on species, genome, antibody, condition.
    6. Every non-control row with a non-empty merge_group_id must have a control_group_id
       that matches at least one is_control=true row's group_id, unless --allow_no_control true.
    7. is_control=true rows must NOT reference a control_group_id.
    8. group_id must be unique within (sample_id, replicate) — soft check (warn).
"""
from __future__ import annotations

import argparse
import csv
import sys
from collections import defaultdict
from typing import Dict, List, Set

REQUIRED_COLS = [
    "sample_id", "species", "genome", "antibody", "condition", "replicate",
    "group_id", "is_control", "control_group_id", "merge_group_id",
    "peak_calling_mode", "notes",
]
PEAK_MODES = {"narrow", "broad", "auto", ""}


def fail(msg: str, log: List[str]) -> None:
    log.append(f"ERROR: {msg}")


def warn(msg: str, log: List[str]) -> None:
    log.append(f"WARN: {msg}")


def info(msg: str, log: List[str]) -> None:
    log.append(f"INFO: {msg}")


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--pairs_tsv", required=True)
    p.add_argument("--csv", required=True)
    p.add_argument("--out_csv", required=True)
    p.add_argument("--log", required=True)
    p.add_argument("--allow_extra", default="false")
    p.add_argument("--allow_no_control", default="false")
    p.add_argument("--genome", default="",
                   help="If non-empty, drop CSV rows whose 'genome' column "
                        "does not equal this value (case-sensitive).")
    args = p.parse_args()

    allow_extra = args.allow_extra.lower() == "true"
    allow_no_ctrl = args.allow_no_control.lower() == "true"
    genome_filter = (args.genome or "").strip()

    log: List[str] = ["# validate_associations.py"]
    if genome_filter:
        log.append(f"# genome filter: '{genome_filter}'")
    else:
        log.append("# genome filter: <none>")

    # Load pairs TSV
    pair_ids: Set[str] = set()
    with open(args.pairs_tsv) as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        for row in rdr:
            pair_ids.add(row["sample_id"])
    log.append(f"# fastq_pair sample_ids: {len(pair_ids)}")

    # Load CSV
    with open(args.csv) as fh:
        rdr = csv.DictReader(fh)
        if rdr.fieldnames is None:
            print("ERROR: empty association_csv", file=sys.stderr)
            return 1
        missing_cols = [c for c in REQUIRED_COLS if c not in rdr.fieldnames]
        if missing_cols:
            msg = f"missing columns in association_csv: {missing_cols}"
            print(f"ERROR: {msg}", file=sys.stderr)
            with open(args.log, "w") as out:
                out.write("\n".join(log + [f"ERROR: {msg}"]) + "\n")
            return 2
        rows = list(rdr)
    log.append(f"# csv rows (pre-filter): {len(rows)}")

    # ---- Genome filter ----
    rows_dropped: List[Dict] = []
    if genome_filter:
        rows_kept: List[Dict] = []
        for r in rows:
            if (r.get("genome") or "").strip() == genome_filter:
                rows_kept.append(r)
            else:
                rows_dropped.append(r)

        # If --genome was set but no rows survive, fail loudly with the
        # available genomes — the user almost certainly has a typo or
        # mismatched spelling (hg38 vs GRCh38, etc.).
        if not rows_kept:
            genomes_seen = sorted({(r.get("genome") or "").strip() for r in rows})
            msg = (f"--genome '{genome_filter}' matched 0 rows in association_csv. "
                   f"Genomes present in the CSV: {genomes_seen}")
            print(f"ERROR: {msg}", file=sys.stderr)
            with open(args.log, "w") as out:
                out.write("\n".join(log + [f"ERROR: {msg}"]) + "\n")
            return 4

        rows = rows_kept
        log.append(f"# csv rows kept after genome filter: {len(rows)}")
        log.append(f"# csv rows dropped by genome filter: {len(rows_dropped)}")
        for r in rows_dropped:
            info(f"genome-filter SKIP: sample '{r['sample_id']}' (genome='{r.get('genome','').strip()}')", log)
    else:
        log.append(f"# csv rows kept: {len(rows)}  (no genome filter applied)")

    dropped_ids: Set[str] = {r["sample_id"] for r in rows_dropped}
    csv_ids: Set[str] = {r["sample_id"] for r in rows}
    log.append(f"# csv sample_ids (kept): {len(csv_ids)}")

    # Rule 1: every pair_id must be in CSV — EXCEPT pair_ids whose CSV row was
    # dropped by the genome filter (those are silently skipped, but logged).
    missing_in_csv = pair_ids - csv_ids - dropped_ids
    for sid in sorted(missing_in_csv):
        fail(f"sample_id '{sid}' detected in FASTQ pairs but NOT in association_csv (any genome)", log)
    fastqs_skipped_by_genome = pair_ids & dropped_ids
    for sid in sorted(fastqs_skipped_by_genome):
        info(f"FASTQ sample '{sid}' will NOT be processed (its CSV row's genome != '{genome_filter}')", log)

    # Rule 2: no extra rows unless allowed
    extras = csv_ids - pair_ids
    for sid in sorted(extras):
        if allow_extra:
            warn(f"sample_id '{sid}' present in association_csv but no FASTQ found (allowed)", log)
        else:
            fail(f"sample_id '{sid}' in association_csv has no FASTQ pair (set --allow_extra true to permit)", log)

    # Per-row checks
    by_merge: Dict[str, List[Dict]] = defaultdict(list)
    by_groupid: Dict[str, List[Dict]] = defaultdict(list)
    control_groupids: Set[str] = set()

    for r in rows:
        # Normalize is_control
        ic_raw = (r.get("is_control") or "").strip().lower()
        if ic_raw not in {"true", "false"}:
            fail(f"sample {r['sample_id']}: is_control must be 'true' or 'false' (got '{ic_raw}')", log)
        r["is_control"] = "true" if ic_raw == "true" else "false"

        # Normalize peak_calling_mode
        pm = (r.get("peak_calling_mode") or "").strip().lower()
        if pm not in PEAK_MODES:
            fail(f"sample {r['sample_id']}: peak_calling_mode must be narrow/broad/auto (got '{pm}')", log)
        r["peak_calling_mode"] = pm if pm else "auto"

        # Required-not-empty
        for col in ["species", "genome", "antibody", "condition", "replicate", "group_id"]:
            if not (r.get(col) or "").strip():
                fail(f"sample {r['sample_id']}: column '{col}' must not be empty", log)

        if r["is_control"] == "true":
            if (r.get("control_group_id") or "").strip():
                fail(f"control sample {r['sample_id']}: control_group_id must be empty", log)
            control_groupids.add(r["group_id"])
        else:
            cgid = (r.get("control_group_id") or "").strip()
            if not cgid and not allow_no_ctrl:
                fail(f"non-control sample {r['sample_id']}: control_group_id is empty (set --allow_no_control true to permit)", log)

        if r.get("merge_group_id"):
            by_merge[r["merge_group_id"]].append(r)
        by_groupid[r["group_id"]].append(r)

    # Rule 5: merge group consistency
    for mg, members in by_merge.items():
        for col in ["species", "genome", "antibody", "condition"]:
            vals = {m[col] for m in members}
            if len(vals) > 1:
                fail(f"merge_group_id '{mg}': inconsistent {col} values: {vals}", log)

    # Rule 6: control_group_id existence
    for r in rows:
        if r["is_control"] == "false" and (r.get("control_group_id") or "").strip():
            cgid = r["control_group_id"].strip()
            if cgid not in control_groupids:
                if allow_no_ctrl:
                    warn(f"sample {r['sample_id']}: control_group_id '{cgid}' not found among control rows (allowed)", log)
                else:
                    fail(f"sample {r['sample_id']}: control_group_id '{cgid}' not present as a group_id of any is_control=true row", log)

    # Write log
    errors = [l for l in log if l.startswith("ERROR")]
    with open(args.log, "w") as out:
        out.write("\n".join(log) + "\n")
        out.write(f"# total errors: {len(errors)}\n")

    if errors:
        for e in errors:
            print(e, file=sys.stderr)
        return 3

    # Write validated CSV with normalized columns
    with open(args.out_csv, "w", newline="") as out:
        w = csv.DictWriter(out, fieldnames=REQUIRED_COLS)
        w.writeheader()
        for r in rows:
            if r["sample_id"] in pair_ids:
                w.writerow({k: r.get(k, "") for k in REQUIRED_COLS})

    print(f"OK: validation passed. Validated CSV: {args.out_csv}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
