#!/usr/bin/env python3
"""Compare frozen and recomputed MSH2 metric tables.

All numeric fields, including nested bootstrap metadata and confidence
intervals, are compared with an absolute tolerance. Cohort keys and other text
fields must match exactly.
"""
from __future__ import annotations

import argparse
import ast
import csv
import json
import math
from pathlib import Path
from typing import Any


TABLES = {
    "primary_metrics.csv": ("held_rep", "source_reps", "endpoint"),
    "sensitivity_metrics.csv": (
        "held_rep",
        "source_reps",
        "endpoint",
        "scenario",
        "pseudocount",
        "depth_frequency_threshold",
        "count_correction",
    ),
}


def read_csv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        return list(reader.fieldnames or []), list(reader)


def compare_objects(left: Any, right: Any, where: str, issues: list[dict], atol: float) -> None:
    if left is None or right is None:
        if left is not right:
            issues.append({"where": where, "left": left, "right": right, "reason": "null mismatch"})
        return
    if isinstance(left, bool) or isinstance(right, bool):
        if type(left) is not type(right) or left != right:
            issues.append({"where": where, "left": left, "right": right, "reason": "boolean mismatch"})
        return
    if isinstance(left, (int, float)) and isinstance(right, (int, float)):
        if not (math.isfinite(float(left)) and math.isfinite(float(right))):
            equal = left == right
        else:
            equal = abs(float(left) - float(right)) <= atol
        if not equal:
            issues.append({
                "where": where,
                "left": left,
                "right": right,
                "abs_diff": abs(float(left) - float(right)) if math.isfinite(float(left)) and math.isfinite(float(right)) else None,
                "reason": "numeric mismatch",
            })
        return
    if type(left) is not type(right):
        issues.append({"where": where, "left_type": type(left).__name__, "right_type": type(right).__name__, "reason": "type mismatch"})
        return
    if isinstance(left, dict):
        if set(left) != set(right):
            issues.append({"where": where, "left_keys": sorted(left), "right_keys": sorted(right), "reason": "dict key mismatch"})
        for key in sorted(set(left) & set(right)):
            compare_objects(left[key], right[key], f"{where}.{key}", issues, atol)
        return
    if isinstance(left, (list, tuple)):
        if len(left) != len(right):
            issues.append({"where": where, "left_len": len(left), "right_len": len(right), "reason": "sequence length mismatch"})
        for index, (a, b) in enumerate(zip(left, right)):
            compare_objects(a, b, f"{where}[{index}]", issues, atol)
        return
    if left != right:
        issues.append({"where": where, "left": left, "right": right, "reason": "text mismatch"})


def compare_cells(left: str, right: str, where: str, issues: list[dict], atol: float) -> None:
    if left == right:
        return
    if left == "" or right == "":
        issues.append({"where": where, "left": left, "right": right, "reason": "empty/value mismatch"})
        return
    try:
        left_num, right_num = float(left), float(right)
    except ValueError:
        left_num = right_num = None
    if left_num is not None and right_num is not None:
        compare_objects(left_num, right_num, where, issues, atol)
        return
    if left[:1] in "[{" and right[:1] in "[{":
        try:
            left_obj, right_obj = ast.literal_eval(left), ast.literal_eval(right)
        except (SyntaxError, ValueError):
            pass
        else:
            compare_objects(left_obj, right_obj, where, issues, atol)
            return
    issues.append({"where": where, "left": left, "right": right, "reason": "text mismatch"})


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-dir", required=True, type=Path, help="directory containing the frozen MSH2 metric CSVs")
    parser.add_argument("--recomputed-dir", required=True, type=Path, help="directory containing the recomputed MSH2 metric CSVs")
    parser.add_argument("--output-json", type=Path, help="optional path for the machine-readable comparison report")
    parser.add_argument("--atol", type=float, default=1e-10, help="absolute tolerance for numeric values (default: 1e-10)")
    args = parser.parse_args()
    if args.atol < 0:
        parser.error("--atol must be nonnegative")

    result: dict[str, Any] = {"tolerance_abs": args.atol, "files": {}, "issues": []}
    expected_rows = {"primary_metrics.csv": 9, "sensitivity_metrics.csv": 45}
    for filename, row_keys in TABLES.items():
        reference_header, reference_rows = read_csv(args.reference_dir / filename)
        recomputed_header, recomputed_rows = read_csv(args.recomputed_dir / filename)
        issues: list[dict] = []
        if reference_header != recomputed_header:
            issues.append({"reason": "header mismatch", "left": reference_header, "right": recomputed_header})
        reference_map: dict[tuple[str, ...], dict[str, str]] = {}
        recomputed_map: dict[tuple[str, ...], dict[str, str]] = {}
        for label, rows, target in (
            ("reference", reference_rows, reference_map),
            ("recomputed", recomputed_rows, recomputed_map),
        ):
            for row in rows:
                key = tuple(row.get(column, "") for column in row_keys)
                if key in target:
                    issues.append({"reason": f"{label} duplicate row key", "key": key})
                target[key] = row
        if set(reference_map) != set(recomputed_map):
            issues.append({
                "reason": "row key/cohort mismatch",
                "reference_only": sorted(set(reference_map) - set(recomputed_map)),
                "recomputed_only": sorted(set(recomputed_map) - set(reference_map)),
            })
        for key in sorted(set(reference_map) & set(recomputed_map)):
            for column in reference_header:
                compare_cells(
                    reference_map[key].get(column, ""),
                    recomputed_map[key].get(column, ""),
                    f"{filename}:{key!r}:{column}",
                    issues,
                    args.atol,
                )
        wanted = expected_rows[filename]
        if len(reference_rows) != wanted or len(recomputed_rows) != wanted:
            issues.append({
                "reason": "unexpected row count",
                "expected_each": wanted,
                "reference": len(reference_rows),
                "recomputed": len(recomputed_rows),
            })
        result["files"][filename] = {
            "reference_rows": len(reference_rows),
            "recomputed_rows": len(recomputed_rows),
            "columns": len(reference_header),
            "row_key": list(row_keys),
            "mismatch_count": len(issues),
            "mismatches": issues[:30],
        }
        result["issues"].extend({"file": filename, **issue} for issue in issues)

    _, primary_rows = read_csv(args.recomputed_dir / "primary_metrics.csv")
    ci_rows = []
    for row in primary_rows:
        raw = row.get("position_cluster_refit_bootstrap", "")
        if not raw:
            continue
        bootstrap = ast.literal_eval(raw)
        ci_rows.append({
            "held_rep": int(row["held_rep"]),
            "source_reps": row["source_reps"],
            "endpoint": row["endpoint"],
            "seed": bootstrap.get("seed"),
            "requested_draws": bootstrap.get("requested_draws"),
            "successful_draws": bootstrap.get("successful_draws"),
            "relative_rmse_gain_95pct_ci": bootstrap.get("relative_rmse_gain_95pct_ci"),
            "paired_mse_difference_full_minus_aa_95pct_ci": bootstrap.get("paired_mse_difference_full_minus_aa_95pct_ci"),
            "rmse_difference_full_minus_aa_95pct_ci": bootstrap.get("rmse_difference_full_minus_aa_95pct_ci"),
        })
    result["primary_bootstrap_ci_rows"] = ci_rows
    if len(ci_rows) != 3:
        result["issues"].append({"reason": "expected exactly 3 primary bootstrap CI rows", "observed": len(ci_rows)})
    for row in ci_rows:
        if row["requested_draws"] != 5000 or row["successful_draws"] != 5000:
            result["issues"].append({"reason": "unexpected primary bootstrap draw count", "row": row})
    result["outcome"] = "PASS" if not result["issues"] else "FAIL"

    if args.output_json:
        args.output_json.parent.mkdir(parents=True, exist_ok=True)
        args.output_json.write_text(json.dumps(result, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "outcome": result["outcome"],
        "primary_rows": result["files"]["primary_metrics.csv"]["recomputed_rows"],
        "sensitivity_rows": result["files"]["sensitivity_metrics.csv"]["recomputed_rows"],
        "primary_ci_rows": len(ci_rows),
        "mismatch_count": len(result["issues"]),
        "ci_seeds": [row["seed"] for row in ci_rows],
        "output_json": str(args.output_json) if args.output_json else None,
    }, indent=2))
    return 0 if result["outcome"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
