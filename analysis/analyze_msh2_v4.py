#!/usr/bin/env python3
"""Source-only MSH2 replicate-held-out analysis of deposited GEO counts.

This reconstructs per-DNA count-table frequency ratios. It does not reproduce
Jia et al.'s amino-acid aggregation/sign filters or assert an author-equivalent
clinical LoF score.
"""
from __future__ import annotations

import gzip
import argparse
import hashlib
import json
import re
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[1]
INPUT = REPO_ROOT / "data/external/GSE162130_MSH2_HAP1_raw_count.tsv.gz"
OUTPUT = REPO_ROOT / "results/msh2_v4"
REPS = (1, 2, 3)
ENDPOINTS = {
    "6TG_mock": ("db6", "mock"),
    "mock_P0": ("mock", "p0"),
    "DB6_P0": ("db6", "p0"),
}
COLS = {
    r: {"p0": f"R{r}D_P0", "db6": f"R{r}DB6_P2", "mock": f"R{r}DB_P2"}
    for r in REPS
}
WT_COL = "WT_plasmid"
ALL_COUNT_COLS = [COLS[r][c] for r in REPS for c in ("p0", "db6", "mock")] + [WT_COL]
AA_RE = re.compile(r"^([A-Z*])_(\d+)_([A-Z*])$")
NT_RE = re.compile(r"^([ACGT]{3})_(\d+)_([ACGT]{3})$")
TILE_RE = re.compile(r"^Tile(\d+)_nt_wt$")
FOLDS = 5
FOLD_SEED = 20261123
BOOT_N = 5000
BOOT_SEED = 20261124
PRIMARY_FREQ = 1.0 / 200000.0
DEPTH_SENSITIVITIES = (1.0 / 400000.0, 1.0 / 100000.0)
PRIMARY_PC = 0.5
PC_SENSITIVITIES = (0.1, 1.0)
MOCK_QC_CUTOFF = -2.0


def log(message: str) -> None:
    stamp = datetime.now(timezone.utc).isoformat()
    print(f"[{stamp}] {message}", flush=True)


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def file_rows(path: Path) -> int:
    n = 0
    with gzip.open(path, "rb") as f:
        for _ in f:
            n += 1
    return n


def assign_tiles(raw: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    raw = raw.copy()
    marker = raw["var_nt"].astype("string").str.extract(TILE_RE, expand=False)
    raw["tile_num"] = pd.to_numeric(marker, errors="coerce").ffill()
    if raw["tile_num"].isna().any():
        raise ValueError("Found raw rows before first TileXX_nt_wt marker")
    raw["tile_num"] = raw["tile_num"].astype(int)
    marker_rows = marker.notna()
    tile_order = raw.loc[marker_rows, "tile_num"].astype(int).tolist()
    if tile_order != list(range(1, 22)):
        raise ValueError(f"Expected ordered tile markers 1..21; observed {tile_order}")

    aa = raw["var_AA"].astype("string").str.extract(AA_RE)
    nt = raw["var_nt"].astype("string").str.extract(NT_RE)
    valid = (aa[0].notna() & nt[0].notna() & (aa[0] != aa[2]) &
             (aa[2] != "*") & (aa[1] == nt[1]) & ~marker_rows)
    miss = raw.loc[valid].copy().reset_index(drop=True)
    miss["ref"] = aa.loc[valid, 0].astype(str).to_numpy()
    miss["pos"] = aa.loc[valid, 1].astype(int).to_numpy()
    miss["alt"] = aa.loc[valid, 2].astype(str).to_numpy()
    miss["nt"] = miss["var_nt"].astype(str)
    miss["group_id"] = (miss["pos"].astype(str) + ":" + miss["ref"] + ">" + miss["alt"])
    if miss["nt"].duplicated().any():
        raise ValueError("Duplicate missense DNA key in the deposited count table")
    group_tile_n = miss.groupby("group_id")["tile_num"].nunique()
    if (group_tile_n != 1).any():
        raise ValueError("At least one amino-acid group crosses tile blocks")

    tile_records = []
    tile_total: dict[int, dict[str, float]] = {}
    tile_wt_total: dict[int, float] = {}
    for tile, block in raw.groupby("tile_num", sort=True):
        tile = int(tile)
        block_marker = block["var_nt"].astype(str).str.match(TILE_RE)
        if int(block_marker.sum()) != 1 or not bool(block_marker.iloc[0]):
            raise ValueError(f"Tile {tile} block must begin with one WT marker")
        count_totals = {}
        for col in ALL_COUNT_COLS:
            values = pd.to_numeric(block[col], errors="coerce")
            count_totals[col] = float(values.sum()) if values.notna().all() else float("nan")
        wt_values = pd.to_numeric(block[WT_COL], errors="coerce")
        tile_wt_total[tile] = float(wt_values.sum()) if wt_values.notna().all() else float("nan")
        tile_total[tile] = count_totals
        tile_miss = miss.loc[miss["tile_num"].eq(tile)]
        pos = tile_miss["pos"]
        row = {
            "tile": tile,
            "count_table_rows_including_tile_wt": int(len(block)),
            "missense_dna_variants": int(len(tile_miss)),
            "missense_groups": int(tile_miss["group_id"].nunique()),
            "position_min": int(pos.min()) if len(pos) else None,
            "position_max": int(pos.max()) if len(pos) else None,
            "wt_plasmid_count_table_total": tile_wt_total[tile],
        }
        for r in REPS:
            for condition in ("p0", "db6", "mock"):
                row[f"R{r}_{condition}_count_table_total"] = count_totals[COLS[r][condition]]
        tile_records.append(row)
    tile_totals = pd.DataFrame(tile_records)
    if tile_totals["wt_plasmid_count_table_total"].isna().any() or \
       (tile_totals["wt_plasmid_count_table_total"] <= 0).any():
        raise ValueError("WT-plasmid tile count-table totals are missing or non-positive")
    for tile, totals in tile_total.items():
        if any((not np.isfinite(v) or v <= 0) for k, v in totals.items() if k != WT_COL):
            raise ValueError(f"Missing/non-positive sample count-table total in tile {tile}")
    provenance = {
        "tile_marker_rows": int(marker_rows.sum()),
        "tile_order": tile_order,
        "tile_assignment": "contiguous rows after TileXX_nt_wt through the next tile marker",
        "denominator_definition": (
            "sum of deposited count-table rows within tile, including the tile WT row; "
            "not a raw sequencer read total"
        ),
        "wt_denominator_definition": (
            "sum of WT_plasmid counts in the same tile block and deposited row universe"
        ),
    }
    return raw, miss, {"tile_totals": tile_totals, "tile_total": tile_total,
                       "tile_wt_total": tile_wt_total, "provenance": provenance}


def build_count_arrays(miss: pd.DataFrame, tile_meta: dict) -> tuple[dict, dict]:
    tile_total = tile_meta["tile_total"]
    tile_wt_total = tile_meta["tile_wt_total"]
    tile_ids = miss["tile_num"].to_numpy(int)
    wt_variant = pd.to_numeric(miss[WT_COL], errors="coerce").to_numpy(float)
    corrected: dict[int, dict[str, np.ndarray]] = {r: {} for r in REPS}
    raw_counts: dict[int, dict[str, np.ndarray]] = {r: {} for r in REPS}
    for r in REPS:
        for condition in ("p0", "db6", "mock"):
            col = COLS[r][condition]
            raw_value = pd.to_numeric(miss[col], errors="coerce").to_numpy(float)
            corr = raw_value.copy()
            for tile in sorted(tile_total):
                idx = tile_ids == tile
                n_sample = tile_total[tile][col]
                n_wt = tile_wt_total[tile]
                good = idx & np.isfinite(raw_value) & np.isfinite(wt_variant)
                adjusted = raw_value[good] - wt_variant[good] * n_sample / n_wt
                adjusted[adjusted < 0] = 0.9
                corr[good] = adjusted
                corr[idx & ~good] = np.nan
            raw_counts[r][condition] = raw_value
            corrected[r][condition] = corr
    return corrected, raw_counts


def totals_for_rows(miss: pd.DataFrame, tile_meta: dict, rep: int, condition: str) -> np.ndarray:
    col = COLS[rep][condition]
    totals = tile_meta["tile_total"]
    return np.array([totals[int(tile)][col] for tile in miss["tile_num"]], dtype=float)


def score_matrix(miss: pd.DataFrame, counts: dict, tile_meta: dict,
                 pseudocount: float) -> dict[str, dict[int, np.ndarray]]:
    result = {endpoint: {} for endpoint in ENDPOINTS}
    for endpoint, (num_cond, den_cond) in ENDPOINTS.items():
        for rep in REPS:
            num = counts[rep][num_cond]
            den = counts[rep][den_cond]
            n_num = totals_for_rows(miss, tile_meta, rep, num_cond)
            n_den = totals_for_rows(miss, tile_meta, rep, den_cond)
            good = (np.isfinite(num) & np.isfinite(den) & np.isfinite(n_num) &
                    np.isfinite(n_den) & (n_num > 0) & (n_den > 0))
            score = np.full(len(miss), np.nan, dtype=float)
            score[good] = (np.log2((num[good] + pseudocount) / n_num[good]) -
                           np.log2((den[good] + pseudocount) / n_den[good]))
            result[endpoint][rep] = score
    return result


def frequency_arrays(miss: pd.DataFrame, counts: dict, tile_meta: dict) -> dict[int, dict[str, np.ndarray]]:
    result = {r: {} for r in REPS}
    for r in REPS:
        for condition in ("p0", "db6", "mock"):
            values = counts[r][condition]
            total = totals_for_rows(miss, tile_meta, r, condition)
            freq = np.full(len(miss), np.nan, dtype=float)
            good = np.isfinite(values) & np.isfinite(total) & (total > 0)
            freq[good] = values[good] / total[good]
            result[r][condition] = freq
    return result


def make_fold_map(positions: np.ndarray) -> dict[int, int]:
    unique = np.array(sorted(set(map(int, positions))), dtype=int)
    rng = np.random.default_rng(FOLD_SEED)
    shuffled = rng.permutation(unique)
    return {int(pos): fold for fold, part in enumerate(np.array_split(shuffled, FOLDS))
            for pos in part}


def source_eligibility(miss: pd.DataFrame, corrected: dict, raw_counts: dict, tile_meta: dict,
                       frequencies: dict, held: int, threshold: float,
                       primary_scores: dict) -> tuple[np.ndarray, dict]:
    sources = [r for r in REPS if r != held]
    p0_pass = np.ones(len(miss), dtype=bool)
    all_depth_pass = np.ones(len(miss), dtype=bool)
    raw_nonzero_pass = np.ones(len(miss), dtype=bool)
    mock_qc_pass = np.ones(len(miss), dtype=bool)
    for r in sources:
        p0_pass &= np.isfinite(frequencies[r]["p0"]) & (frequencies[r]["p0"] > threshold)
        for condition in ("p0", "mock", "db6"):
            f = frequencies[r][condition]
            all_depth_pass &= np.isfinite(f) & (f > threshold)
            raw_nonzero_pass &= np.isfinite(raw_counts[r][condition]) & (raw_counts[r][condition] > 0)
        mock_p0 = primary_scores["mock_P0"][r]
        mock_qc_pass &= np.isfinite(mock_p0) & (mock_p0 > MOCK_QC_CUTOFF)
    source_measurable = all_depth_pass & raw_nonzero_pass
    eligible = source_measurable & mock_qc_pass
    floor_only_threshold_pass = np.zeros(len(miss), dtype=bool)
    floor_source_total_min = np.inf
    final_zero_detail = {}
    for r in sources:
        for condition in ("p0", "mock", "db6"):
            raw_value = raw_counts[r][condition]
            corr_value = corrected[r][condition]
            f = frequencies[r][condition]
            raw_zero = ~np.isfinite(raw_value) | (raw_value == 0)
            floor_only = raw_zero & np.isfinite(corr_value) & np.isclose(corr_value, 0.9)
            floor_only_threshold_pass |= floor_only & np.isfinite(f) & (f > threshold)
            final_zero_detail[f"final_eligible_raw_zero_R{r}_{condition}"] = int(np.sum(eligible & raw_zero))
            final_zero_detail[f"thresholdpass_raw_zero_R{r}_{condition}"] = int(np.sum(all_depth_pass & raw_zero))
            totals = tile_meta["tile_total"]
            min_total = min(totals[int(tile)][COLS[r][condition]] for tile in miss["tile_num"].unique())
            floor_source_total_min = min(floor_source_total_min, min_total)
    detail = {
        "source_reps": ",".join(map(str, sources)),
        "candidate_missense_variants": int(len(miss)),
        "candidate_amino_acid_groups": int(miss["group_id"].nunique()),
        "observed_nonzero_any_source_condition": int(sum(
            np.any(np.column_stack([raw_counts[r][c] > 0 for r in sources
                                    for c in ("p0", "mock", "db6")]), axis=1))),
        "p0_frequency_pass_variants": int(p0_pass.sum()),
        "all_three_condition_frequency_pass_variants": int(all_depth_pass.sum()),
        "all_source_raw_conditions_nonzero_pass_variants": int(raw_nonzero_pass.sum()),
        "source_measurable_after_frequency_and_raw_nonzero": int(source_measurable.sum()),
        "source_mock_qc_pass_variants": int(mock_qc_pass.sum()),
        "source_eligible_after_mock_qc_variants": int(eligible.sum()),
        "source_eligible_variants": int(eligible.sum()),
        "frequency_qualified_with_raw_zero_in_source": int(np.sum(all_depth_pass & ~raw_nonzero_pass)),
        "frequency_qualified_floor_only_in_source": int(floor_only_threshold_pass.sum()),
        "minimum_source_tile_count_table_total": float(floor_source_total_min),
    }
    detail.update(final_zero_detail)
    return eligible, detail


def add_group_size_filter(miss: pd.DataFrame, eligible: np.ndarray) -> tuple[np.ndarray, dict]:
    eligible_groups = miss.loc[eligible].groupby("group_id")["nt"].nunique()
    multi_group_ids = eligible_groups[eligible_groups >= 2].index
    multi = eligible & miss["group_id"].isin(multi_group_ids).to_numpy()
    detail = {
        "eligible_groups_any_size": int(len(eligible_groups)),
        "eligible_groups_multidna": int(len(multi_group_ids)),
        "source_eligible_multidna_variants": int(multi.sum()),
        "source_eligible_multidna_positions": int(miss.loc[multi, "pos"].nunique()),
    }
    return multi, detail


def build_source_rows(miss: pd.DataFrame, source_mask: np.ndarray, held: int,
                      scores: dict, fold_map: dict[int, int]) -> pd.DataFrame:
    sources = [r for r in REPS if r != held]
    a, b = sources
    frame = miss.loc[source_mask, ["var_nt", "group_id", "pos", "tile_num", "nt"]].copy()
    frame["source_row_index"] = np.flatnonzero(source_mask)
    frame = frame.reset_index(drop=True)
    frame["score_a"] = scores["_current"][a][source_mask]
    frame["score_b"] = scores["_current"][b][source_mask]
    frame["pos_fold"] = frame["pos"].map(fold_map).astype(int)
    frame["group_mean_a"] = frame.groupby("group_id")["score_a"].transform("mean")
    frame["group_mean_b"] = frame.groupby("group_id")["score_b"].transform("mean")
    frame["resid_a"] = frame["score_a"] - frame["group_mean_a"]
    frame["resid_b"] = frame["score_b"] - frame["group_mean_b"]
    frame["x_source"] = (frame["resid_a"] + frame["resid_b"]) / 2.0
    frame["pred_aa"] = (frame["group_mean_a"] + frame["group_mean_b"]) / 2.0
    frame["source_cross_variant"] = frame["resid_a"] * frame["resid_b"]
    frame["source_x2_variant"] = frame["x_source"] ** 2
    return frame


def source_betas(frame: pd.DataFrame) -> tuple[dict[int, float], pd.DataFrame, pd.DataFrame]:
    by_group = frame.groupby("group_id", sort=False).agg(
        pos=("pos", "first"), tile=("tile_num", "first"),
        c_source=("source_cross_variant", "mean"),
        v_source=("source_x2_variant", "mean"),
    ).reset_index()
    pos_folds = frame.groupby("pos", sort=False)["pos_fold"].first()
    by_group["fold"] = by_group["pos"].map(pos_folds).astype(int)
    betas = {}
    rows = []
    for fold in range(FOLDS):
        train = by_group[by_group["fold"] != fold]
        c = float(train["c_source"].mean()) if len(train) else float("nan")
        v = float(train["v_source"].mean()) if len(train) else float("nan")
        beta = float(np.clip(c / v, 0.0, 1.0)) if v > 0 and np.isfinite(c / v) else 0.0
        betas[fold] = beta
        rows.append({"fold": fold, "beta": beta, "C_source_group_equal": c,
                     "V_source_group_equal": v, "train_groups": int(len(train)),
                     "train_positions": int(train["pos"].nunique()),
                     "validation_groups": int((by_group["fold"] == fold).sum()),
                     "validation_positions": int(by_group.loc[by_group["fold"] == fold, "pos"].nunique())})
    return betas, pd.DataFrame(rows), by_group


def bootstrap_refit(groups: pd.DataFrame, source_groups: pd.DataFrame,
                    n_boot: int = BOOT_N,
                    seed: int = BOOT_SEED) -> dict:
    """Position-cluster bootstrap; refit source beta within every resampled fold."""
    # Resample the entire source-eligible position set. A position with no
    # complete target group can still contribute to source-only beta refitting.
    positions = np.array(sorted(source_groups["pos"].unique()), dtype=int)
    if not len(positions):
        return {"n_draws": 0}
    pos_to_i = {int(p): i for i, p in enumerate(positions)}
    p_index = np.array([pos_to_i[int(p)] for p in groups["pos"]], dtype=int)
    n_eval_group_pos = np.bincount(p_index, minlength=len(positions)).astype(float)
    sum_mse_aa = np.bincount(p_index, weights=groups["mse_aa"].to_numpy(float), minlength=len(positions))
    sum_cross = np.bincount(p_index, weights=groups["cross"].to_numpy(float), minlength=len(positions))
    sum_x2 = np.bincount(p_index, weights=groups["x2"].to_numpy(float), minlength=len(positions))
    source_subset = source_groups.copy()
    source_subset["_pidx"] = source_subset["pos"].map(pos_to_i).astype(int)
    n_source_group_pos = np.bincount(source_subset["_pidx"].to_numpy(int), minlength=len(positions)).astype(float)
    sum_c = np.bincount(source_subset["_pidx"].to_numpy(int),
                         weights=source_subset["c_source"].to_numpy(float), minlength=len(positions))
    sum_v = np.bincount(source_subset["_pidx"].to_numpy(int),
                         weights=source_subset["v_source"].to_numpy(float), minlength=len(positions))
    pos_fold = np.array([int(source_subset.loc[source_subset["pos"].eq(p), "fold"].iloc[0]) for p in positions])
    rng = np.random.default_rng(seed)
    out_gain, out_delta, out_rmse_diff = [], [], []
    for draw in range(n_boot):
        sampled = rng.integers(0, len(positions), size=len(positions))
        w_pos = np.bincount(sampled, minlength=len(positions)).astype(float)
        aa_num = full_num = denom = 0.0
        for f in range(FOLDS):
            train_pos = pos_fold != f
            train_w = w_pos[train_pos] * n_source_group_pos[train_pos]
            if train_w.sum() <= 0:
                continue
            c = float(np.sum(w_pos[train_pos] * sum_c[train_pos]) / train_w.sum())
            v = float(np.sum(w_pos[train_pos] * sum_v[train_pos]) / train_w.sum())
            beta = float(np.clip(c / v, 0.0, 1.0)) if v > 0 and np.isfinite(c / v) else 0.0
            test_pos = pos_fold == f
            wt = w_pos[test_pos]
            if wt.sum() <= 0:
                continue
            aa_num += float(np.sum(wt * sum_mse_aa[test_pos]))
            full_num += float(np.sum(wt * (sum_mse_aa[test_pos] -
                                           2.0 * beta * sum_cross[test_pos] +
                                           beta * beta * sum_x2[test_pos])))
            denom += float(np.sum(wt * n_eval_group_pos[test_pos]))
        if denom <= 0 or aa_num <= 0 or full_num < 0:
            continue
        mse_aa, mse_full = aa_num / denom, full_num / denom
        out_delta.append(mse_full - mse_aa)
        out_gain.append(1.0 - np.sqrt(mse_full / mse_aa))
        out_rmse_diff.append(np.sqrt(mse_full) - np.sqrt(mse_aa))
        if (draw + 1) % 500 == 0:
            log(f"bootstrap progress {draw + 1}/{n_boot}")
    def ci(values):
        return [float(x) for x in np.quantile(np.asarray(values), [0.025, 0.975])]
    return {
        "requested_draws": int(n_boot),
        "successful_draws": int(len(out_gain)),
        "seed": int(seed),
        "refit_beta_each_fold_each_draw": True,
        "cluster": "protein position, sampled with replacement; all AA groups at a position retained",
        "conditional_scope": "within-study positions and fixed count table; not independent-pool uncertainty",
        "relative_rmse_gain_95pct_ci": ci(out_gain) if out_gain else None,
        "paired_mse_difference_full_minus_aa_95pct_ci": ci(out_delta) if out_delta else None,
        "rmse_difference_full_minus_aa_95pct_ci": ci(out_rmse_diff) if out_rmse_diff else None,
    }


def evaluate(miss: pd.DataFrame, eligible: np.ndarray, held: int, endpoint: str,
             endpoint_scores: dict, all_endpoint_scores: dict, fold_map: dict[int, int],
             complete_mask: np.ndarray, run_bootstrap: bool = False,
             bootstrap_seed: int = BOOT_SEED) -> tuple[dict, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    sources = [r for r in REPS if r != held]
    # All eligible DNA records define the source model. Target completeness only
    # gates evaluation rows and never affects fitting or beta calibration.
    source_scores = {"_current": {r: endpoint_scores[r] for r in sources}}
    source_frame = build_source_rows(miss, eligible, held, source_scores, fold_map)
    beta_by_fold, beta_table, source_group_stats = source_betas(source_frame)
    source_frame["beta"] = source_frame["pos_fold"].map(beta_by_fold).astype(float)
    source_frame["pred_full"] = source_frame["pred_aa"] + source_frame["beta"] * source_frame["x_source"]
    source_frame["target"] = all_endpoint_scores[endpoint][held][eligible]
    source_frame["target_complete"] = complete_mask[eligible]
    eval_rows = source_frame.loc[source_frame["target_complete"]].copy()
    if eval_rows.empty:
        raise ValueError(f"No complete target groups for held replicate {held}, endpoint {endpoint}")
    group_rows = []
    for group_id, g in eval_rows.groupby("group_id", sort=False):
        y = g["target"].to_numpy(float)
        aa = g["pred_aa"].to_numpy(float)
        full = g["pred_full"].to_numpy(float)
        x = g["x_source"].to_numpy(float)
        ea = g["resid_a"].to_numpy(float)
        eb = g["resid_b"].to_numpy(float)
        group_rows.append({
            "group_id": group_id,
            "pos": int(g["pos"].iloc[0]),
            "tile": int(g["tile_num"].iloc[0]),
            "fold": int(g["pos_fold"].iloc[0]),
            "n_variants": int(len(g)),
            "mse_aa": float(np.mean((y - aa) ** 2)),
            "mae_aa": float(np.mean(np.abs(y - aa))),
            "mse_full": float(np.mean((y - full) ** 2)),
            "mae_full": float(np.mean(np.abs(y - full))),
            "cross": float(np.mean((y - aa) * x)),
            "x2": float(np.mean(x ** 2)),
            "c_source": float(np.mean(ea * eb)),
            "v_source": float(np.mean(x ** 2)),
        })
    groups = pd.DataFrame(group_rows)
    mse_aa = float(groups["mse_aa"].mean())
    mse_full = float(groups["mse_full"].mean())
    rmse_aa, rmse_full = float(np.sqrt(mse_aa)), float(np.sqrt(mse_full))
    summary = {
        "held_rep": int(held),
        "source_reps": ",".join(map(str, sources)),
        "endpoint": endpoint,
        "n_groups": int(len(groups)),
        "n_variants": int(eval_rows.shape[0]),
        "n_positions": int(groups["pos"].nunique()),
        "n_tiles": int(groups["tile"].nunique()),
        "mse_aa": mse_aa,
        "mse_full": mse_full,
        "rmse_aa": rmse_aa,
        "rmse_full": rmse_full,
        "rmse_difference_full_minus_aa": rmse_full - rmse_aa,
        "mae_aa_group_equal": float(groups["mae_aa"].mean()),
        "mae_full_group_equal": float(groups["mae_full"].mean()),
        "paired_mse_difference_full_minus_aa": mse_full - mse_aa,
        "relative_rmse_gain": 1.0 - rmse_full / rmse_aa if rmse_aa else None,
        "relative_mse_gain": 1.0 - mse_full / mse_aa if mse_aa else None,
        "beta_min": float(beta_table["beta"].min()),
        "beta_max": float(beta_table["beta"].max()),
        "beta_mean_unweighted_fold_descriptive": float(beta_table["beta"].mean()),
    }
    ci = None
    if run_bootstrap:
        ci = bootstrap_refit(groups, source_group_stats, n_boot=BOOT_N, seed=bootstrap_seed)
        summary["position_cluster_refit_bootstrap"] = ci
        summary["relative_rmse_gain_95pct_ci"] = ci.get("relative_rmse_gain_95pct_ci")
        summary["paired_mse_difference_full_minus_aa_95pct_ci"] = ci.get(
            "paired_mse_difference_full_minus_aa_95pct_ci")
        summary["rmse_difference_full_minus_aa_95pct_ci"] = ci.get(
            "rmse_difference_full_minus_aa_95pct_ci")
    tile_rows = []
    for tile, g in groups.groupby("tile", sort=True):
        mse_a_tile, mse_f_tile = float(g["mse_aa"].mean()), float(g["mse_full"].mean())
        tile_rows.append({
            "held_rep": held, "source_reps": ",".join(map(str, sources)),
            "endpoint": endpoint, "tile": int(tile), "n_groups": int(len(g)),
            "n_positions": int(g["pos"].nunique()), "n_variants": int(g["n_variants"].sum()),
            "rmse_aa": float(np.sqrt(mse_a_tile)), "rmse_full": float(np.sqrt(mse_f_tile)),
            "paired_mse_difference_full_minus_aa": mse_f_tile - mse_a_tile,
            "relative_rmse_gain": 1.0 - np.sqrt(mse_f_tile / mse_a_tile) if mse_a_tile else None,
        })
    eval_rows["endpoint"] = endpoint
    eval_rows["held_rep"] = held
    eval_rows["source_reps"] = ",".join(map(str, sources))
    eval_rows["target_zero_num"] = False
    eval_rows["target_zero_den"] = False
    num_cond, den_cond = ENDPOINTS[endpoint]
    eval_source_idx = eval_rows["source_row_index"].to_numpy(int)
    eval_rows["target_zero_num"] = pd.to_numeric(
        miss.iloc[eval_source_idx][COLS[held][num_cond]], errors="coerce").to_numpy() == 0
    eval_rows["target_zero_den"] = pd.to_numeric(
        miss.iloc[eval_source_idx][COLS[held][den_cond]], errors="coerce").to_numpy() == 0
    eval_rows["target_score"] = eval_rows["target"].to_numpy(float)
    beta_table["held_rep"] = held
    beta_table["source_reps"] = ",".join(map(str, sources))
    beta_table["endpoint"] = endpoint
    return summary, groups, pd.DataFrame(tile_rows), eval_rows, beta_table


def write_markdown(metrics: pd.DataFrame, attrition: pd.DataFrame, zero: pd.DataFrame,
                   source: dict, output_path: Path) -> None:
    lines = [
        "# MSH2 source-only reanalysis (revision v4)", "",
        "## Methods",
        "The analyzed count table is the public GEO GSE162130 supplement `GSE162130_MSH2_HAP1_raw_count.tsv.gz` (SHA256 recorded in `provenance.json`). It is byte-identical to the prior incoming count table, but the analysis now uses the deposited R1–R3 P0, passage-2 mock (DB), passage-2 6-TG (DB6), and WT_plasmid columns. The deposited file contains 21 `TileXX_nt_wt` marker rows followed by the tile's enumerated codon rows. For every tile and sample, the denominator was defined as the sum of the deposited count-table rows in that tile block, including its WT row; the WT control denominator was the sum of WT_plasmid rows in the same tile block. These are count-table totals, not raw sequencer read totals.",
        "",
        "For each tile/sample, WT error counts were reconstructed by subtracting the mutation-specific WT_plasmid count multiplied by sample count-table total / WT_plasmid count-table total. Negative adjusted counts were set to 0.9, following the paper's stated floor; frequency-ratio scores used a 0.5-count pseudocount. The primary score is the per-DNA log2 ratio of corrected count-table frequencies for 6-TG/mock. Mock/P0 and DB6/P0 were evaluated in parallel. This is a reconstruction from deposited tile-count totals, not a claim to reproduce Jia et al.'s full amino-acid LoF pipeline.",
        "",
        f"In each replicate-held-out round, source-only eligibility required both source replicates to have corrected P0, mock, and 6-TG frequencies > {PRIMARY_FREQ:.8g} (1/200,000), raw counts >0 in all six source condition columns (so the 0.9 error-correction floor alone could not qualify a DNA key), and source mock/P0 score > {MOCK_QC_CUTOFF:g}; eligible amino-acid groups needed at least two distinct DNA variants. These criteria were fixed before examining held-out outcomes. Target raw zero counts were retained; missing counts were not imputed. Only groups with all source-eligible DNA members observed in all three endpoint conditions contributed to a matched evaluation set. The same-AA predictor was the mean source group score; the DNA predictor added a cross-replicate residual mean multiplied by a source-only beta. Five deterministic protein-position folds (seed {FOLD_SEED}) held all groups at a position together. Metrics equally weighted AA groups, with variants equally weighted within group. The primary 6-TG/mock gain interval used {BOOT_N:,} position-cluster draws with beta refit inside each fold/draw.",
        "",
        "## Results",
        "The per-held-replicate estimates below are descriptive internal transfer within the same published study. The three held-out rounds are not independent external validations. Positive relative RMSE gain favors adding within-AA DNA information; negative gain favors the AA-only group-mean baseline.",
        "",
        "### Primary and parallel endpoints",
        "| Held replicate | Endpoint | Groups | Variants | AA RMSE | DNA RMSE | AA MAE | DNA MAE | Paired MSE difference (DNA−AA) | Relative RMSE gain | 95% position-cluster refit CI |",
        "|---:|:---|---:|---:|---:|---:|---:|---:|---:|---:|:---|",
    ]
    for _, row in metrics.sort_values(["held_rep", "endpoint"]).iterrows():
        ci = row.get("relative_rmse_gain_95pct_ci")
        if isinstance(ci, list):
            ci_text = f"[{ci[0]:.4f}, {ci[1]:.4f}]"
        else:
            ci_text = "—"
        lines.append(
            f"| {int(row.held_rep)} | {row.endpoint} | {int(row.n_groups)} | {int(row.n_variants)} | "
            f"{row.rmse_aa:.4f} | {row.rmse_full:.4f} | {row.mae_aa_group_equal:.4f} | "
            f"{row.mae_full_group_equal:.4f} | {row.paired_mse_difference_full_minus_aa:.4f} | "
            f"{row.relative_rmse_gain:.4%} | {ci_text} |"
        )
    lines.extend(["", "### Source-to-target attrition", "",
                  "See `attrition.csv` for candidate enumeration, source frequency/mock-QC attrition, multi-DNA groups, and common target-complete evaluation counts for every held replicate. See `target_zero_missing.csv` for measured target zeros separately from missing counts.",
                  "", "### Interpretation limits", "",
                  "This reanalysis estimates within-study, replicate-held-out transfer among variants already enumerated in the same MSH2 library. The source dishes share a tile-specific starting library; the uncertainty intervals resample protein positions conditional on this fixed dataset and do not represent independent transduction pools. WT_plasmid is a shared reference across replicate samples, so the WT-correction step does not isolate biological function from common technical signal. The 6-TG/mock endpoint is primary; mock/P0 and DB6/P0 are parallel endpoints using the same source-defined cohort. No external prediction, unseen-variant performance, clinical interpretation, or isolated biological mechanism is established.",
                  ""])
    output_path.write_text("\n".join(lines), encoding="utf-8")


def main(argv: list[str] | None = None) -> None:
    global INPUT, OUTPUT
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=INPUT,
                        help="GSE162130 count table (.tsv.gz); no network download is attempted")
    parser.add_argument("--output", type=Path, default=OUTPUT,
                        help="directory for computed output (default: repo/results/msh2_v4)")
    args = parser.parse_args(argv)
    INPUT = args.input.expanduser().resolve()
    OUTPUT = args.output.expanduser().resolve()
    if not INPUT.is_file():
        raise FileNotFoundError(f"MSH2 source not found: {INPUT}; provide --input")
    OUTPUT.mkdir(parents=True, exist_ok=True)
    log("loading published GEO count table")
    if not INPUT.is_file():
        raise FileNotFoundError(INPUT)
    raw = pd.read_csv(INPUT, sep="\t", compression="gzip", low_memory=False)
    required = ["var_nt", "var_AA", *ALL_COUNT_COLS]
    missing_columns = [c for c in required if c not in raw.columns]
    if missing_columns:
        raise ValueError(f"Missing columns: {missing_columns}")
    for col in ALL_COUNT_COLS:
        raw[col] = pd.to_numeric(raw[col], errors="coerce")
    raw, miss, tile_meta = assign_tiles(raw)
    tile_meta["tile_totals"].to_csv(OUTPUT / "tile_totals.csv", index=False)
    log(f"loaded {len(raw):,} rows; {len(miss):,} missense DNA rows; {miss.group_id.nunique():,} AA groups")
    corrected, raw_counts = build_count_arrays(miss, tile_meta)
    frequencies = frequency_arrays(miss, corrected, tile_meta)
    primary_scores = score_matrix(miss, corrected, tile_meta, PRIMARY_PC)
    fold_map = make_fold_map(miss["pos"].unique())
    pd.DataFrame([{"position": p, "fold": f} for p, f in fold_map.items()]).sort_values("position").to_csv(
        OUTPUT / "position_folds.csv", index=False)

    # Base source-only cohort counts and a common target completeness mask.
    attr_rows, source_masks, common_complete = [], {}, {}
    for held in REPS:
        eligible, detail = source_eligibility(miss, corrected, raw_counts, tile_meta, frequencies,
                                             held, PRIMARY_FREQ, primary_scores)
        multi, group_detail = add_group_size_filter(miss, eligible)
        detail.update(group_detail)
        source_masks[held] = multi
        all_scores = {ep: {r: primary_scores[ep][r] for r in REPS} for ep in ENDPOINTS}
        perrow_complete = np.ones(len(miss), dtype=bool)
        for ep in ENDPOINTS:
            perrow_complete &= np.isfinite(all_scores[ep][held])
        group_complete = pd.DataFrame({"group_id": miss.loc[multi, "group_id"],
                                       "complete": perrow_complete[multi]}).groupby("group_id")["complete"].all()
        comp = miss["group_id"].map(group_complete).fillna(False).to_numpy(bool)
        common_complete[held] = comp
        detail.update({"held_rep": held, "source_reps": ",".join(str(r) for r in REPS if r != held),
                       "depth_frequency_threshold": PRIMARY_FREQ,
                       "threshold_per_200k_reads": 1.0 / PRIMARY_FREQ,
                       "mock_qc_log2_cutoff": MOCK_QC_CUTOFF,
                       "candidate_positions": int(miss["pos"].nunique()),
                       "target_zero_counts_retained": True,
                       "target_missing_dna_rows_not_imputed": True,
                       "common_target_complete_groups": int(miss.loc[multi & comp, "group_id"].nunique()),
                       "common_target_complete_variants": int((multi & comp).sum()),
                       "common_target_complete_positions": int(miss.loc[multi & comp, "pos"].nunique())})
        attr_rows.append(detail)

    metrics, tile_rows, prediction_rows, beta_rows = [], [], [], []
    group_rows_all = []
    for held in REPS:
        for endpoint in ENDPOINTS:
            log(f"primary held-replicate analysis: held={held}, endpoint={endpoint}")
            result, groups, tiles, preds, betas = evaluate(
                miss, source_masks[held], held, endpoint, primary_scores[endpoint],
                primary_scores, fold_map, common_complete[held],
                run_bootstrap=(endpoint == "6TG_mock"), bootstrap_seed=BOOT_SEED + held)
            metrics.append(result)
            groups["held_rep"] = held
            groups["endpoint"] = endpoint
            groups_all_row = groups.to_dict(orient="records")
            group_rows_all.extend(groups_all_row)
            tiles["scenario"] = "primary_WT_corrected_pc0.5"
            tile_rows.extend(tiles.to_dict(orient="records"))
            preds["scenario"] = "primary_WT_corrected_pc0.5"
            prediction_rows.extend(preds.to_dict(orient="records"))
            betas["scenario"] = "primary_WT_corrected_pc0.5"
            beta_rows.extend(betas.to_dict(orient="records"))

    metrics_df = pd.DataFrame(metrics)
    attrition_df = pd.DataFrame(attr_rows)
    tile_metrics_df = pd.DataFrame(tile_rows)
    predictions_df = pd.DataFrame(prediction_rows)
    group_errors_df = pd.DataFrame(group_rows_all)
    betas_df = pd.DataFrame(beta_rows)
    metrics_df.to_csv(OUTPUT / "primary_metrics.csv", index=False)
    attrition_df.to_csv(OUTPUT / "attrition.csv", index=False)
    tile_metrics_df.to_csv(OUTPUT / "tile_metrics.csv", index=False)
    predictions_df.to_csv(OUTPUT / "primary_variant_predictions.csv.gz", index=False, compression="gzip")
    group_errors_df.to_csv(OUTPUT / "primary_group_errors.csv.gz", index=False, compression="gzip")
    betas_df.to_csv(OUTPUT / "fold_betas.csv", index=False)

    # Prespecified point-estimate sensitivities only; the bootstrap belongs to
    # the primary WT-corrected, 0.5-count, 1/200,000 cohort.
    sensitivity_rows = []
    for pc in PC_SENSITIVITIES:
        sens_scores = score_matrix(miss, corrected, tile_meta, pc)
        for held in REPS:
            for endpoint in ENDPOINTS:
                result, *_ = evaluate(miss, source_masks[held], held, endpoint,
                                      sens_scores[endpoint], sens_scores, fold_map,
                                      common_complete[held], run_bootstrap=False)
                result.update({"scenario": "pseudocount_sensitivity", "pseudocount": pc,
                               "depth_frequency_threshold": PRIMARY_FREQ,
                               "count_correction": "WT_error_subtracted"})
                sensitivity_rows.append(result)
    for threshold in DEPTH_SENSITIVITIES:
        sens_primary_scores = score_matrix(miss, corrected, tile_meta, PRIMARY_PC)
        for held in REPS:
            eligible, _ = source_eligibility(miss, corrected, raw_counts, tile_meta, frequencies,
                                             held, threshold, sens_primary_scores)
            multi, _ = add_group_size_filter(miss, eligible)
            perrow_complete = np.ones(len(miss), dtype=bool)
            for ep in ENDPOINTS:
                perrow_complete &= np.isfinite(sens_primary_scores[ep][held])
            group_complete = pd.DataFrame({"group_id": miss.loc[multi, "group_id"],
                                           "complete": perrow_complete[multi]}).groupby("group_id")["complete"].all()
            sens_complete = miss["group_id"].map(group_complete).fillna(False).to_numpy(bool)
            for endpoint in ENDPOINTS:
                result, *_ = evaluate(miss, multi, held, endpoint,
                                      sens_primary_scores[endpoint], sens_primary_scores,
                                      fold_map, sens_complete, run_bootstrap=False)
                result.update({"scenario": "depth_threshold_sensitivity", "pseudocount": PRIMARY_PC,
                               "depth_frequency_threshold": threshold,
                               "count_correction": "WT_error_subtracted"})
                sensitivity_rows.append(result)
    # WT-error-correction sensitivity: same primary source eligibility and
    # evaluation membership, but scores are calculated from uncorrected counts.
    uncorrected_scores = score_matrix(miss, raw_counts, tile_meta, PRIMARY_PC)
    for held in REPS:
        for endpoint in ENDPOINTS:
            result, *_ = evaluate(miss, source_masks[held], held, endpoint,
                                  uncorrected_scores[endpoint], uncorrected_scores,
                                  fold_map, common_complete[held], run_bootstrap=False)
            result.update({"scenario": "uncorrected_count_sensitivity_same_primary_cohort",
                           "pseudocount": PRIMARY_PC, "depth_frequency_threshold": PRIMARY_FREQ,
                           "count_correction": "uncorrected"})
            sensitivity_rows.append(result)
    sensitivity_df = pd.DataFrame(sensitivity_rows)
    sensitivity_df.to_csv(OUTPUT / "sensitivity_metrics.csv", index=False)

    # Explicit target-zero/missing audit over source-eligible multi-DNA rows.
    zero_rows = []
    for held in REPS:
        idx = source_masks[held]
        for endpoint, (num_cond, den_cond) in ENDPOINTS.items():
            num = raw_counts[held][num_cond][idx]
            den = raw_counts[held][den_cond][idx]
            target_score = primary_scores[endpoint][held][idx]
            zero_rows.append({
                "held_rep": held, "source_reps": ",".join(str(r) for r in REPS if r != held),
                "endpoint": endpoint, "source_eligible_multidna_variants": int(idx.sum()),
                "target_zero_numerator_variants": int(np.sum(np.isfinite(num) & (num == 0))),
                "target_zero_denominator_variants": int(np.sum(np.isfinite(den) & (den == 0))),
                "target_zero_both_variants": int(np.sum(np.isfinite(num) & np.isfinite(den) & (num == 0) & (den == 0))),
                "target_missing_numerator_variants": int(np.sum(~np.isfinite(num))),
                "target_missing_denominator_variants": int(np.sum(~np.isfinite(den))),
                "target_missing_score_variants": int(np.sum(~np.isfinite(target_score))),
                "target_complete_common_groups": int(miss.loc[idx & common_complete[held], "group_id"].nunique()),
                "target_complete_common_variants": int(np.sum(idx & common_complete[held])),
            })
    zero_df = pd.DataFrame(zero_rows)
    zero_df.to_csv(OUTPUT / "target_zero_missing.csv", index=False)

    source_hash = sha256_file(INPUT)
    script_path = Path(__file__)
    provenance = {
        "analysis": "MSH2 source-only replicate-held-out endpoint reconstruction",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "input": {
            "geo_accession": "GSE162130",
            "filename": INPUT.name,
            "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE162nnn/GSE162130/suppl/GSE162130_MSH2_HAP1_raw_count.tsv.gz",
            "size_bytes": INPUT.stat().st_size,
            "sha256": source_hash,
            "decompressed_rows_including_header": file_rows(INPUT),
            "header_columns": required,
            "rows_incoming_precheck_sha256_match": source_hash == "9e315efaf715db4f906aaf5b82b4175847af819d5f65659f03e5070d710e2299",
        },
        "script": {"path": str(script_path), "sha256": sha256_file(script_path)},
        "score_and_selection": {
            "primary_endpoint": "log2(corrected 6TG count-table frequency / corrected mock count-table frequency)",
            "parallel_endpoints": list(ENDPOINTS),
            "frequency_threshold": PRIMARY_FREQ,
            "frequency_threshold_sensitivities": list(DEPTH_SENSITIVITIES),
            "pseudocount_primary_corrected_reads": PRIMARY_PC,
            "pseudocount_sensitivities_corrected_reads": list(PC_SENSITIVITIES),
            "source_mock_qc": f"source mock/P0 log2 score > {MOCK_QC_CUTOFF:g}",
            "raw_source_nonzero_rule": "raw P0, mock, and 6-TG counts >0 in both source replicates; prevents the 0.9 WT correction floor alone from qualifying a source variant",
            "group_rule": "at least two distinct DNA keys after source-only eligibility",
            "target_count_used_for_eligibility_or_beta": False,
            "target_raw_zero_retained": True,
            "target_missing_imputed": False,
            "jia_opposite_sign_or_AA_aggregation_filters_applied": False,
            "WT_error_correction": "sample count minus mutation-specific WT_plasmid count scaled by same-tile sample/WT count-table total; negative adjusted counts set to 0.9",
            "WT_correction_scope": "shared WT_plasmid reference, conditional reconstruction; not a causal/biological isolation",
            "denominator": tile_meta["provenance"]["denominator_definition"],
            "WT_denominator": tile_meta["provenance"]["wt_denominator_definition"],
        },
        "validation": {
            "outer": "leave one of three culture replicates out; source eligibility from the other two",
            "position_folds": FOLDS,
            "position_fold_seed": FOLD_SEED,
            "source_beta": "cross-replicate within-AA residual covariance / source residual variance, clipped to [0,1], fit outside each validation-position fold",
            "baseline": "same-budget mean of two source replicate AA group means",
            "weighting": "AA groups equally weighted; DNA variants equally weighted within group",
            "bootstrap": {"draws": BOOT_N, "seed_base": BOOT_SEED, "cluster": "protein position", "beta_refit": True},
            "inference_limit": "conditional within-study position uncertainty; three held-replicate rounds are overlapping internal transfer, not independent external validation",
        },
        "tile_layout_and_denominator": tile_meta["provenance"],
        "input_policy": "provided through --input; the analysis does not download source data",
    }
    (OUTPUT / "provenance.json").write_text(json.dumps(provenance, indent=2), encoding="utf-8")
    write_markdown(metrics_df, attrition_df, zero_df, provenance, OUTPUT / "manuscript_ready_methods_results.md")
    json_tables = {
        "primary_metrics": json.loads(metrics_df.to_json(orient="records", double_precision=15)),
        "sensitivity_metrics": json.loads(sensitivity_df.to_json(orient="records", double_precision=15)),
        "attrition": json.loads(attrition_df.to_json(orient="records", double_precision=15)),
        "target_zero_missing": json.loads(zero_df.to_json(orient="records", double_precision=15)),
        "tile_metrics": json.loads(tile_metrics_df.to_json(orient="records", double_precision=15)),
    }
    (OUTPUT / "summary.json").write_text(json.dumps({
        "analysis": "MSH2 source-only held-replicate endpoint reconstruction",
        "provenance_file": "provenance.json",
        "primary_metrics_csv": "primary_metrics.csv",
        "summary_tables": json_tables,
    }, indent=2, allow_nan=False), encoding="utf-8")
    log(f"complete: {len(metrics_df)} primary estimates, {len(sensitivity_df)} sensitivity rows")
    log(f"outputs: {OUTPUT}")


if __name__ == "__main__":
    main()
