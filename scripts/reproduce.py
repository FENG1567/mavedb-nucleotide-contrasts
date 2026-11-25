#!/usr/bin/env python3
"""Portable entry point for selected analysis and figure rebuilds.

Examples are listed in the repository README. Analysis commands forward their
remaining options to the corresponding self-written script; no data are fetched.
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
COMMANDS = {
    "msh2": REPO_ROOT / "analysis/analyze_msh2_v4.py",
    "brca2-map": REPO_ROOT / "analysis/rebuild_brca2_mapping.py",
    "brca2": REPO_ROOT / "analysis/brca2_core.py",
    "brca2-direction": REPO_ROOT / "analysis/brca2_direction_ties.py",
    "selection": REPO_ROOT / "analysis/selection_bias_null.py",
    "collection": REPO_ROOT / "analysis/verify_collection_reconciliation.py",
    "collection-raw": REPO_ROOT / "analysis/audit_collection_groups_precomputed.py",
    "splicing": REPO_ROOT / "analysis/splicing_fixed_pair_audit_v2.py",
}


def run_figures(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(prog="reproduce.py figures")
    parser.add_argument("--data-root", type=Path, default=REPO_ROOT / "data/derived")
    parser.add_argument("--output-dir", type=Path, default=REPO_ROOT / "results/figures")
    parser.add_argument("--source-dir", type=Path, default=REPO_ROOT / "results/figures/source_data")
    args = parser.parse_args(argv)
    builder = REPO_ROOT / "figures/build_display_figures.py"
    subprocess.run([sys.executable, str(builder), "--data-root", str(args.data_root),
                    "--output-dir", str(args.output_dir), "--source-dir", str(args.source_dir)], check=True)
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=["figures", *COMMANDS])
    # Parse only the command so child --help reaches the selected program.
    tokens = sys.argv[1:] if argv is None else argv
    args = parser.parse_args(tokens[:1])
    rest = tokens[1:]
    if args.command == "figures":
        return run_figures(rest)
    script = COMMANDS[args.command]
    return subprocess.run([sys.executable, str(script), *rest], check=False).returncode


if __name__ == "__main__":
    raise SystemExit(main())
