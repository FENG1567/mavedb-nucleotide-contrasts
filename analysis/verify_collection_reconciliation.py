#!/usr/bin/env python3
"""Reaggregate and verify the frozen 92/89/27 collection reconciliation tables.

This command consumes the packaged aggregate CSVs only. The optional raw-source
audit is a separate command because its 27 raw score exports are not bundled.
"""
from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_REVISION = REPO_ROOT / "data/derived/revision_20261123_review/collection"
DEFAULT_LEGACY = REPO_ROOT / "data/derived/revision_20261122/dna_audit/legacy_reproduction.json"
DEFAULT_OUT = REPO_ROOT / "results/collection_reconciliation.json"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def as_int(value: str | None) -> int:
    try:
        return int(float(value or 0))
    except (TypeError, ValueError):
        return 0


def state_counts(rows: list[dict[str, str]], field: str) -> dict[str, int]:
    return dict(sorted(Counter(row.get(field, "") for row in rows).items()))


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--collection-dir", type=Path, default=DEFAULT_REVISION,
                        help="directory with group, per-set and source-provenance CSVs")
    parser.add_argument("--legacy-summary", type=Path, default=DEFAULT_LEGACY)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args(argv)
    collection_dir = args.collection_dir.expanduser().resolve()
    paths = {
        "groups": collection_dir / "group_level_empirical_testability.csv",
        "per_set": collection_dir / "per_score_set_reconciliation.csv",
        "source": collection_dir / "source_provenance_per_score_set.csv",
    }
    missing = [str(path) for path in [*paths.values(), args.legacy_summary.resolve()] if not path.is_file()]
    if missing:
        raise SystemExit("Missing frozen collection input(s):\n  - " + "\n  - ".join(missing))
    groups, per_set, source = (read_csv(paths[name]) for name in ("groups", "per_set", "source"))
    legacy = json.loads(args.legacy_summary.expanduser().resolve().read_text(encoding="utf-8"))
    expected = legacy.get("expected", legacy.get("observed", {}))
    observed = {
        "manifest_score_sets": len(per_set),
        "analyzable_score_sets": sum(as_int(row.get("legacy_missense_rows")) > 0 for row in per_set),
        "source_score_sets_checked": len(source),
        "groups_with_multiple_legacy_keys": len(groups),
        "legacy_numeric_missense_rows": sum(as_int(row.get("legacy_missense_rows")) for row in per_set),
        "legacy_unique_keys_in_multikey_groups": sum(as_int(row.get("legacy_unique_keys_in_ge2_groups")) for row in per_set),
        "legacy_multikey_score_sets": sum(as_int(row.get("legacy_groups_ge2_unique_keys")) > 0 for row in per_set),
        "source_exact_numeric_multiset_matches": sum(str(row.get("exact_numeric_record_multiset_match", "")).lower() in {"true", "1", "yes"} for row in source),
        "source_numeric_rows_unmatched": sum(as_int(row.get("raw_numeric_rows_unmatched")) for row in source),
        "mapped_numeric_rows_without_raw_match": sum(as_int(row.get("mapped_numeric_rows_without_raw_match")) for row in source),
        "tpk1_unresolved_legacy_multikey_groups": sum(as_int(row.get("TPK1_full_construct_unresolved_groups")) for row in per_set if row.get("gene") == "TPK1"),
        "tpk1_unresolved_legacy_keys_in_multikey_groups": sum(as_int(row.get("legacy_unique_keys_in_ge2_groups")) for row in per_set if row.get("gene") == "TPK1"),
    }
    observed["accession_reference_group_states"] = state_counts(groups, "accession_reference_group_state")
    observed["versioned_request_group_states"] = state_counts(groups, "versioned_request_group_state")
    observed["target_sequence_group_states"] = state_counts(groups, "target_sequence_group_state")
    checks = {
        "historical_manifest_score_sets_92": observed["manifest_score_sets"] == 92,
        "historical_analyzable_score_sets_89": observed["analyzable_score_sets"] == 89,
        "source_score_sets_checked_27": observed["source_score_sets_checked"] == 27,
        "historical_multikey_groups_13041": observed["groups_with_multiple_legacy_keys"] == 13041,
        "historical_missense_rows_413090": observed["legacy_numeric_missense_rows"] == 413090,
        "historical_keys_31723": observed["legacy_unique_keys_in_multikey_groups"] == 31723,
        "multi_key_source_sets_27": observed["legacy_multikey_score_sets"] == 27,
        "all_27_numeric_source_matches": observed["source_exact_numeric_multiset_matches"] == 27 and observed["source_numeric_rows_unmatched"] == 0 and observed["mapped_numeric_rows_without_raw_match"] == 0,
        "tpk1_unresolved_multikey_groups_2474": observed["tpk1_unresolved_legacy_multikey_groups"] == 2474,
        "tpk1_unresolved_legacy_keys_6913": observed["tpk1_unresolved_legacy_keys_in_multikey_groups"] == 6913,
    }
    if not all(checks.values()):
        raise SystemExit("Collection frozen-table checks failed: " + json.dumps(checks, sort_keys=True))
    report = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "historical_expected": expected,
        "observed": observed,
        "reconciliation_checks": checks,
        "input_mode": "frozen aggregate CSV reaggregation; raw score records were not reopened",
        "independent_measurement_boundary": {
            "status": "not identified by these archived score and sequence qualification tables",
            "independent_groups_proven": None,
            "conservative_confirmed_lower_bound": 0,
        },
    }
    output = args.output.expanduser().resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, ensure_ascii=False), flush=True)


if __name__ == "__main__":
    main()
