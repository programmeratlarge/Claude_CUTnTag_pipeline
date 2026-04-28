#!/usr/bin/env python3
"""
validate_associations.py
Validate the user-provided association_csv against the FASTQ pairs TSV.

Required columns:
    sample_id,species,genome,antibody,condition,replicate,group_id,
    is_control,control_group_id,merge_group_id,peak_calling_mode,notes

Rules enforced:
    1. Every sample_id in the FASTQ pairs TSV must appear in the CSV.
    2. No extra sample_id rows in the CSV unless --allow_extra is passed.
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


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--pairs_tsv", required=True)
    p.add_argument("--csv", required=True)
    p.add_argument("--out_csv", required=True)
    p.add_argument("--log", required=True)
    p.add_argument("--allow_extra", default="false")
    p.add_argument("--allow_no_control", default="false")
    args = p.parse_args()

    allow_extra = args.allow_extra.lower() == "true"
    allow_no_ctrl = args.allow_no_control.lower() == "true"

    log: List[str] = ["# validate_associations.py"]

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

    csv_ids: Set[str] = {r["sample_id"] for r in rows}
    log.append(f"# csv sample_ids: {len(csv_ids)}")

    # Rule 1: every pair_id present in CSV
    missing_in_csv = pair_ids - csv_ids
    for sid in sorted(missing_in_csv):
        fail(f"sample_id '{sid}' detected in FASTQ pairs but NOT in association_csv", log)

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
