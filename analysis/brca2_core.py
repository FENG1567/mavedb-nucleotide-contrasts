#!/usr/bin/env python3
"""Recompute the paired BRCA2 primary Phi from the two frozen mapped score tables.

The primary estimate uses per-assay z scores over all mapped missense variants,
same-amino-acid groups with at least two shared HGVS DNA keys, sample covariance
within groups, and Phi = 1 - sum(group covariance) / sum(group mean variance).
The original 5,000-draw ordinary group-bootstrap stream uses seed 20261114.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT_DIR = REPO_ROOT / "data/external/brca2_cc0"
DEFAULT_OUT = REPO_ROOT / "results/brca2_core"
EXPECTED_PHI = 0.8750878921486126


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_mapped(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, dtype=str, keep_default_na=False)
    required = {"score", "protein_pos", "ref_aa", "alt_aa", "class", "hgvs_nt"}
    missing = sorted(required - set(frame.columns))
    if missing:
        raise ValueError(f"{path} is missing required columns: {missing}")
    frame["score"] = pd.to_numeric(frame["score"], errors="coerce")
    frame["protein_pos"] = pd.to_numeric(frame["protein_pos"], errors="coerce")
    frame = frame.dropna(subset=["score", "protein_pos"])
    frame = frame[frame["class"].astype(str).eq("missense")].copy()
    if frame.empty:
        raise ValueError(f"No mapped missense rows in {path}")
    frame["protein_pos"] = frame["protein_pos"].astype(int)
    sd = float(frame["score"].std(ddof=0))
    frame["z"] = (frame["score"] - frame["score"].mean()) / sd if sd else 0.0
    return frame


def group_values(frame: pd.DataFrame) -> dict[tuple[int, str, str], dict[str, float]]:
    grouped: dict[tuple[int, str, str], dict[str, float]] = defaultdict(dict)
    for row in frame.itertuples(index=False):
        grouped[(int(row.protein_pos), str(row.ref_aa), str(row.alt_aa))][str(row.hgvs_nt)] = float(row.z)
    return grouped


def paired_groups(hap1: pd.DataFrame, mes: pd.DataFrame) -> pd.DataFrame:
    a, b = group_values(hap1), group_values(mes)
    rows = []
    for key in a:
        if key not in b:
            continue
        dna_keys = sorted(set(a[key]) & set(b[key]))
        if len(dna_keys) < 2:
            continue
        x = np.asarray([a[key][dna] for dna in dna_keys], dtype=float)
        y = np.asarray([b[key][dna] for dna in dna_keys], dtype=float)
        xc, yc = x - x.mean(), y - y.mean()
        n = len(dna_keys)
        var_hap1 = float(np.dot(xc, xc) / (n - 1))
        var_mes = float(np.dot(yc, yc) / (n - 1))
        cov = float(np.dot(xc, yc) / (n - 1))
        if n == 2:
            d_hap1, d_mes = float(x[0] - x[1]), float(y[0] - y[1])
        else:
            order = np.argsort(x)
            d_hap1 = float(x[order[-1]] - x[order[0]])
            d_mes = float(y[order[-1]] - y[order[0]])
        rows.append({
            "group_id": f"BRCA2:{key[0]}:{key[1]}>{key[2]}",
            "protein_pos": key[0], "ref_aa": key[1], "alt_aa": key[2],
            "n_variants": n, "dna_keys": ";".join(dna_keys),
            "var_hap1": var_hap1, "var_mes": var_mes,
            "mean_var": (var_hap1 + var_mes) / 2.0, "cov": cov,
            "d_hap1": d_hap1, "d_mes": d_mes,
            "sign_concordant": bool(d_hap1 * d_mes > 0),
        })
    return pd.DataFrame(rows)


def bootstrap_summary(groups: pd.DataFrame, draws: int = 5000) -> dict:
    if len(groups) != 388 or int(groups["n_variants"].sum()) != 823:
        raise ValueError(
            "Expected frozen BRCA2 core of 388 groups / 823 shared DNA keys; "
            f"observed {len(groups)} / {int(groups['n_variants'].sum())}"
        )
    cov = groups["cov"].to_numpy(dtype=float)
    mean_var = groups["mean_var"].to_numpy(dtype=float)
    phi = float(1.0 - cov.sum() / mean_var.sum())

    # Phase A called the ratio bootstrap first; retain its RNG sequence exactly.
    rng = np.random.default_rng(20261114)
    idx = rng.integers(0, len(groups), size=(draws, len(groups)))
    ratios = cov[idx].sum(axis=1) / np.maximum(mean_var[idx].sum(axis=1), 1e-12)
    phi_ci = [float(1 - np.quantile(ratios, 0.975)), float(1 - np.quantile(ratios, 0.025))]

    by_position = groups.groupby("protein_pos", sort=False).agg(
        group_n=("cov", "size"), cov_sum=("cov", "sum"),
        mean_var_sum=("mean_var", "sum"), concordant_n=("sign_concordant", "sum")
    ).reset_index()
    n_positions = len(by_position)
    rng_cluster = np.random.default_rng(20261122)
    weights = rng_cluster.multinomial(n_positions, np.full(n_positions, 1 / n_positions), size=draws)
    sampled_n = weights @ by_position["group_n"].to_numpy(dtype=float)
    cov_sum = weights @ by_position["cov_sum"].to_numpy(dtype=float)
    mean_var_sum = weights @ by_position["mean_var_sum"].to_numpy(dtype=float)
    concordant_sum = weights @ by_position["concordant_n"].to_numpy(dtype=float)
    q = lambda x: [float(np.quantile(x, .025)), float(np.quantile(x, .975))]
    return {
        "groups": int(len(groups)), "protein_positions": int(n_positions),
        "shared_dna_variants": int(groups["n_variants"].sum()),
        "mean_cov": float(cov.mean()), "phi_A_primary": phi,
        "phi_A_primary_ci_original_group_bootstrap": phi_ci,
        "direction_R": float(groups["sign_concordant"].mean()),
        "position_cluster_bootstrap_5000_seed_20261122": {
            "phi_ci": q(1 - cov_sum / mean_var_sum),
            "direction_R_ci": q(concordant_sum / sampled_n),
            "bootstrap_group_count_range": [int(sampled_n.min()), int(sampled_n.max())],
        },
    }


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--hap1", type=Path, default=DEFAULT_INPUT_DIR / "mapped_BRCA2_00001225-a-1.csv")
    parser.add_argument("--mes", type=Path, default=DEFAULT_INPUT_DIR / "mapped_BRCA2_00001242-a-1.csv")
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--bootstrap-draws", type=int, default=5000)
    args = parser.parse_args(argv)
    hap1_path, mes_path = args.hap1.expanduser().resolve(), args.mes.expanduser().resolve()
    for path in (hap1_path, mes_path):
        if not path.is_file():
            raise FileNotFoundError(f"Mapped BRCA2 input not found: {path}")
    out = args.output_dir.expanduser().resolve()
    out.mkdir(parents=True, exist_ok=True)
    groups = paired_groups(load_mapped(hap1_path), load_mapped(mes_path))
    summary = bootstrap_summary(groups, args.bootstrap_draws)
    summary["primary_phi_matches_frozen_value"] = abs(summary["phi_A_primary"] - EXPECTED_PHI) <= 1e-12
    summary["inputs"] = {
        "hap1": {"path": str(hap1_path), "sha256": sha256(hap1_path)},
        "mes": {"path": str(mes_path), "sha256": sha256(mes_path)},
    }
    groups.to_csv(out / "brca2_paired_groups.csv", index=False)
    (out / "brca2_core_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
