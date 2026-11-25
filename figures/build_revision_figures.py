#!/usr/bin/env python3
"""Build the four GENETICS v4 main figures from locked revision outputs.

Run on the analysis server only with the revision's system Python, one thread,
and user-site packages disabled. The script stops before drawing if any required
input or required field is absent. It never synthesizes observations.
"""
from __future__ import annotations

import csv
import argparse
import hashlib
import json
import math
import os
import sys
from pathlib import Path
from typing import Any

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DATA_ROOT = REPO_ROOT / "data/derived"
OUT = REPO_ROOT / "results/figures"
SOURCE = REPO_ROOT / "results/figure_source"
INPUTS: dict[str, Path] = {}
OPTIONAL_INPUTS: dict[str, Path] = {}


def configure_paths(data_root: Path, output_dir: Path, source_dir: Path) -> None:
    global INPUTS, OPTIONAL_INPUTS, OUT, SOURCE
    data_root = data_root.expanduser().resolve()
    revision = data_root / "revision_20261123_review"
    previous = data_root / "revision_20261122"
    INPUTS = {
        "collection_summary": revision / "collection/collection_group_reconciliation.json",
        "collection_groups": revision / "collection/group_level_empirical_testability.csv",
        "collection_per_set": revision / "collection/per_score_set_reconciliation.csv",
        "historical_reproduction": previous / "dna_audit/legacy_reproduction.json",
        "brca2_summary": revision / "pairs/brca2_direction_summary_20261123.json",
        "brca2_groups": revision / "pairs/brca2_direction_group_sensitivities_20261123.csv",
        "brca2_covariance_history": previous / "statistics/brca2_covariance_summary.json",
        "msh2_summary": revision / "msh2/output/summary.json",
        "msh2_metrics": revision / "msh2/output/primary_metrics.csv",
        "msh2_attrition": revision / "msh2/output/attrition.csv",
        "msh2_group_errors": revision / "msh2/output/primary_group_errors.csv.gz",
        "selection_null_spec": REPO_ROOT / "analysis/selection_bias_null.py",
        "selection_null": data_root / "outputs/selection_bias_null.csv",
        "real_workflow": data_root / "outputs/real_workflow_comparison.csv",
    }
    OPTIONAL_INPUTS = {
        "msh2_sensitivity": revision / "msh2/output/sensitivity_metrics.csv",
    }
    OUT = output_dir.expanduser().resolve()
    SOURCE = source_dir.expanduser().resolve()

BLUE = "#3B6F9C"
GREEN = "#5B8E75"
ORANGE = "#C7864D"
GREY = "#6D7378"
LIGHT = "#D7DEE4"
PALE_BLUE = "#E7EEF4"
PALE_ORANGE = "#F5EADF"
PALE_GREEN = "#E7EFEA"

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
    "font.size": 7,
    "axes.labelsize": 7,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 7,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 0.65,
    "xtick.major.width": 0.65,
    "ytick.major.width": 0.65,
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "savefig.facecolor": "white",
    "figure.facecolor": "white",
})


def fail(message: str) -> None:
    raise SystemExit(f"FIGURE PREFLIGHT FAILED: {message}")


def read_json(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        fail(f"expected a JSON object: {path}")
    return value


def require_columns(frame: pd.DataFrame, columns: list[str], label: str) -> None:
    missing = [column for column in columns if column not in frame.columns]
    if missing:
        fail(f"{label} is missing required columns: {missing}")


def require_keys(value: dict[str, Any], keys: list[str], label: str) -> None:
    missing = [key for key in keys if key not in value]
    if missing:
        fail(f"{label} is missing required keys: {missing}")


def number(value: Any, label: str) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError):
        fail(f"expected numeric {label}; observed {value!r}")
    if not math.isfinite(result):
        fail(f"expected finite {label}; observed {value!r}")
    return result


def bool_series(series: pd.Series) -> pd.Series:
    return series.astype(str).str.strip().str.lower().isin({"true", "1", "yes"})


def csv_list(value: Any, label: str) -> list[float] | None:
    if value is None or (isinstance(value, float) and math.isnan(value)) or str(value).strip() == "":
        return None
    if isinstance(value, (list, tuple)):
        items = value
    else:
        text = str(value).strip().strip("[]")
        items = [part.strip() for part in text.split(",") if part.strip()]
    if len(items) != 2:
        fail(f"expected two interval bounds for {label}; observed {value!r}")
    return [number(items[0], f"{label} lower bound"), number(items[1], f"{label} upper bound")]


def check_inputs() -> None:
    missing = [f"{name}: {path}" for name, path in INPUTS.items() if not path.is_file()]
    if missing:
        print("Required figure inputs are not finalized; no graphics were rendered:", file=sys.stderr)
        for item in missing:
            print(f"  - {item}", file=sys.stderr)
        raise SystemExit(2)


def load_inputs() -> dict[str, Any]:
    check_inputs()
    data: dict[str, Any] = {
        "collection_summary": read_json(INPUTS["collection_summary"]),
        "historical_reproduction": read_json(INPUTS["historical_reproduction"]),
        "collection_groups": pd.read_csv(INPUTS["collection_groups"], dtype=str, keep_default_na=False),
        "collection_per_set": pd.read_csv(INPUTS["collection_per_set"], dtype=str, keep_default_na=False),
        "brca2_summary": read_json(INPUTS["brca2_summary"]),
        "brca2_groups": pd.read_csv(INPUTS["brca2_groups"]),
        "brca2_covariance_history": read_json(INPUTS["brca2_covariance_history"]),
        "msh2_summary": read_json(INPUTS["msh2_summary"]),
        "msh2_metrics": pd.read_csv(INPUTS["msh2_metrics"]),
        "msh2_attrition": pd.read_csv(INPUTS["msh2_attrition"]),
        "msh2_group_errors": pd.read_csv(INPUTS["msh2_group_errors"], compression="gzip"),
        "selection_null": pd.read_csv(INPUTS["selection_null"]),
        "real_workflow": pd.read_csv(INPUTS["real_workflow"]),
    }
    optional = OPTIONAL_INPUTS["msh2_sensitivity"]
    data["msh2_sensitivity"] = pd.read_csv(optional) if optional.is_file() else None
    validate_inputs(data)
    return data


def validate_inputs(data: dict[str, Any]) -> None:
    legacy = data["historical_reproduction"].get("observed", {})
    require_keys(legacy, ["manifest_score_sets", "missense_rows", "unique_keys_in_ge2_groups",
                          "groups_ge2", "group_size_distribution"], "historical reproduction")
    if int(legacy["manifest_score_sets"]) != 92 or int(legacy["missense_rows"]) != 413090 or int(legacy["unique_keys_in_ge2_groups"]) != 31723:
        fail("historical denominators do not match the locked 92 / 413090 / 31723 contract")

    collection_summary = data["collection_summary"]
    require_keys(collection_summary, ["observed", "reconciliation_checks"], "collection reconciliation")
    require_keys(collection_summary["observed"], ["historical_manifest_score_sets", "historical_multi_key_groups",
        "historical_keys", "accession_reference_group_states_all", "versioned_request_group_states_all"],
        "collection reconciliation observed")
    require_columns(data["collection_groups"], ["gene", "accession_reference_group_state",
        "versioned_request_group_state", "target_sequence_group_state",
        "empirical_independent_DNA_measurement_provenance"], "collection group table")
    if len(data["collection_groups"]) != int(collection_summary["observed"]["historical_multi_key_groups"]):
        fail("collection group table row count does not match the reconciliation summary")

    brca = data["brca2_summary"]
    require_keys(brca, ["group_universe", "preserved_primary", "estimators", "diagnostics"], "BRCA2 summary")
    require_columns(data["brca2_groups"], ["protein_pos", "n_variants", "primary_d_Hap1_z", "primary_d_mES_z",
        "primary_concordant", "Hap1_any_exact_score_tie", "mES_any_exact_score_tie",
        "primary_tie_break_sign_ambiguous", "reverse_concordant", "fixed_concordant"], "BRCA2 group table")
    if len(data["brca2_groups"]) != int(brca["group_universe"]["groups"]):
        fail("BRCA2 group table row count does not match the summary")
    if len(data["brca2_groups"]) != 388 or int(pd.to_numeric(data["brca2_groups"]["n_variants"]).sum()) != 823:
        fail("BRCA2 frozen core must contain 388 same-amino-acid groups and 823 paired DNA keys")
    old_stats = data["brca2_covariance_history"]
    require_keys(old_stats, ["phi_A_primary_preserved", "phi_A_primary_ci_preserved", "rank_scale_sensitivity",
                             "class_control_missense_phi_sensitivity_only"], "BRCA2 scale sensitivity summary")
    if not math.isclose(float(old_stats["phi_A_primary_preserved"]), 0.8750878921486126, rel_tol=0, abs_tol=1e-12):
        fail("BRCA2 primary Phi differs from the frozen estimate")

    msh2 = data["msh2_summary"].get("summary_tables", {})
    require_keys(msh2, ["primary_metrics", "attrition"], "MSH2 summary tables")
    require_columns(data["msh2_metrics"], ["held_rep", "endpoint", "n_groups", "relative_rmse_gain",
        "rmse_aa", "rmse_full", "mae_aa_group_equal", "mae_full_group_equal"], "MSH2 primary metrics")
    require_columns(data["msh2_attrition"], ["held_rep", "candidate_missense_variants",
        "candidate_amino_acid_groups", "observed_nonzero_any_source_condition",
        "source_eligible_variants", "eligible_groups_multidna", "common_target_complete_groups"], "MSH2 attrition")
    require_columns(data["msh2_group_errors"], ["held_rep", "endpoint", "group_id", "mse_aa", "mae_aa",
        "mse_full", "mae_full"], "MSH2 group errors")
    if set(data["msh2_metrics"]["endpoint"].astype(str)) != {"6TG_mock", "mock_P0", "DB6_P0"}:
        fail("MSH2 primary metrics do not contain exactly the three locked endpoints")

    require_columns(data["selection_null"], ["tau_over_sigma", "pipeline", "mean_R", "ci_low", "ci_high",
        "P_R_ge_60", "mean_selected"], "selection-null simulation")
    require_columns(data["real_workflow"], ["pipeline", "n_selected", "R", "R_ci_low", "R_ci_high"],
        "real-workflow comparison")
    for path in (INPUTS["selection_null_spec"],):
        if not path.is_file():
            fail(f"missing simulation provenance script: {path}")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def save_figure(fig: mpl.figure.Figure, stem: str) -> None:
    fig.savefig(OUT / f"{stem}.pdf", bbox_inches="tight", metadata={"Creator": "GENETICS revision v4 figure builder"})
    fig.savefig(OUT / f"{stem}.png", dpi=600, bbox_inches="tight")
    fig.savefig(OUT / f"{stem}.svg", bbox_inches="tight")
    plt.close(fig)


def label_panel(ax: plt.Axes, letter: str) -> None:
    ax.text(-0.10, 1.05, letter, transform=ax.transAxes, fontsize=8, fontweight="bold", va="bottom")


def clean_axis(ax: plt.Axes) -> None:
    ax.spines["left"].set_color(GREY)
    ax.spines["bottom"].set_color(GREY)
    ax.tick_params(colors=GREY, length=2.5)
    ax.xaxis.label.set_color(GREY)
    ax.yaxis.label.set_color(GREY)


def source_export(data: dict[str, Any]) -> None:
    SOURCE.mkdir(parents=True, exist_ok=True)
    groups = data["collection_groups"]
    group_tiers = []
    for column, tier in (("accession_reference_group_state", "accession_reference"),
                         ("versioned_request_group_state", "versioned_request"),
                         ("target_sequence_group_state", "target_sequence_conditional")):
        counts = groups[column].value_counts(dropna=False)
        for state, count in counts.items():
            group_tiers.append({"tier": tier, "group_state": state, "n_groups": int(count)})
    pd.DataFrame(group_tiers).to_csv(SOURCE / "figure1_group_tiers.csv", index=False)

    brca_cols = ["group_id", "protein_pos", "n_variants", "primary_d_Hap1_z", "primary_d_mES_z",
        "primary_concordant", "Hap1_any_exact_score_tie", "mES_any_exact_score_tie",
        "primary_tie_break_sign_ambiguous", "reverse_concordant", "fixed_concordant"]
    data["brca2_groups"][brca_cols].to_csv(SOURCE / "figure2_brca2_contrasts.csv", index=False)
    data["msh2_metrics"].to_csv(SOURCE / "figure3_msh2_metrics.csv", index=False)
    data["msh2_attrition"].to_csv(SOURCE / "figure3_msh2_attrition.csv", index=False)
    data["msh2_group_errors"].to_csv(SOURCE / "figure3_msh2_group_errors.csv.gz", index=False, compression="gzip")
    data["selection_null"].to_csv(SOURCE / "figure4_selection_null.csv", index=False)
    data["real_workflow"].to_csv(SOURCE / "figure4_real_workflow.csv", index=False)
    manifest = {
        "generated_utc": __import__("datetime").datetime.now(__import__("datetime").timezone.utc).isoformat(),
        "inputs": {name: {"path": str(path), "sha256": sha256(path)} for name, path in INPUTS.items()},
        "optional_inputs": {name: {"path": str(path), "present": path.is_file(),
            "sha256": sha256(path) if path.is_file() else None} for name, path in OPTIONAL_INPUTS.items()},
        "panel_source_map": {
            "figure1": ["historical_reproduction.observed", "collection_group_reconciliation.json", "group_level_empirical_testability.csv"],
            "figure2": ["brca2_direction_summary_20261123.json", "brca2_direction_group_sensitivities_20261123.csv", "brca2_covariance_summary.json"],
            "figure3": ["msh2/output/summary.json", "primary_metrics.csv", "attrition.csv", "primary_group_errors.csv.gz"],
            "figure4": ["outputs/selection_bias_null.csv", "outputs/real_workflow_comparison.csv", "scripts/selection_bias_null.py"],
        },
    }
    (SOURCE / "input_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    rule = {
        "label": "Rule A two-variant selection acceptance region",
        "predicate_exact": "abs(D1) > 1.96 * sqrt(2 * max((D1 - D2)^2 / 4, 1e-12))",
        "predicate_away_from_origin": "abs(D1 - D2) < (sqrt(2) / 1.96) * abs(D1)",
        "grid_floor_inside_max": 1e-12,
        "display": "exact conceptual acceptance grid; no empirical observations",
    }
    (SOURCE / "figure4_analytic_rule_definition.json").write_text(json.dumps(rule, indent=2) + "\n", encoding="utf-8")


def figure1(data: dict[str, Any]) -> None:
    historical = data["historical_reproduction"]["observed"]
    recon = data["collection_summary"]["observed"]
    groups = data["collection_groups"]
    fig = plt.figure(figsize=(7.2, 4.0), constrained_layout=False)
    gs = fig.add_gridspec(1, 2, left=0.10, right=0.99, bottom=0.29, top=0.92,
                          wspace=0.35, width_ratios=[1.18, 1.55])
    ax0 = fig.add_subplot(gs[0, 0])
    ax1 = fig.add_subplot(gs[0, 1])
    label_panel(ax0, "a")
    label_panel(ax1, "b")

    cards = [
        ("92", "score sets in manifest", BLUE),
        (f"{int(historical['missense_rows']):,}", "missense score rows", GREY),
        (f"{int(historical['unique_keys_in_ge2_groups']):,}", "distinct keys in historical ≥2-key groups", ORANGE),
        (f"{int(historical['groups_ge2']):,}", "historical groups with ≥2 keys", GREEN),
    ]
    ax0.set_axis_off()
    for i, (value, text, color) in enumerate(cards):
        y = 0.94 - i * 0.20
        ax0.text(0.02, y, value, transform=ax0.transAxes, fontsize=12, fontweight="bold", color=color, va="center")
        ax0.text(0.02, y - 0.07, text, transform=ax0.transAxes, fontsize=7, color=GREY, va="center", wrap=True)
    dist = historical["group_size_distribution"]
    sizes = sorted((int(k), int(v)) for k, v in dist.items())
    inset = ax0.inset_axes([0.05, 0.02, 0.90, 0.18])
    size_x = np.arange(len(sizes))
    inset.bar(size_x, [v for _, v in sizes], color=LIGHT, edgecolor="white", linewidth=0.4)
    inset.set_xticks(size_x, [str(k) for k, _ in sizes])
    inset.set_xlabel("historical keys per group", labelpad=1)
    inset.set_ylabel("groups", labelpad=1)
    inset.tick_params(length=2, pad=1)
    for spine in inset.spines.values():
        spine.set_color(LIGHT)

    states = [
        ("Accession reference", "accession_reference_group_state", ["ge2", "one", "zero", "unknown_0_or_1"]),
        ("Versioned request", "versioned_request_group_state", ["ge2", "one", "zero", "unknown_0_or_1"]),
        ("Target-sequence conditional", "target_sequence_group_state",
         ["conditional_ge2", "conditional_one", "zero_target_compatible_under_assumption",
          "not_assessed_target_sequence", "tpk1_full_construct_unresolved"]),
    ]
    colors = [BLUE, GREEN, ORANGE, GREY, LIGHT]
    y_positions = np.arange(len(states))
    left = np.zeros(len(states), dtype=float)
    state_labels = ["≥2", "1", "0 qualifying", "unknown / not assessed", "TPK1 unresolved"]
    for index in range(5):
        vals = []
        for _, col, cats in states:
            state = cats[index] if index < len(cats) else None
            vals.append(int((groups[col] == state).sum()) if state else 0)
        ax1.barh(y_positions, vals, left=left, color=colors[index], edgecolor="white", linewidth=0.45,
                 label=state_labels[index])
        left += np.asarray(vals)
    ax1.set_yticks(y_positions, [item[0] for item in states])
    ax1.invert_yaxis()
    ax1.set_xlabel("candidate score-table groups")
    ax1.set_xlim(0, max(left) * 1.07)
    ax1.grid(axis="x", color=LIGHT, linewidth=0.45, zorder=0)
    ax1.set_axisbelow(True)
    handles, legend_labels = ax1.get_legend_handles_labels()
    fig.legend(handles, legend_labels, frameon=False, ncol=3, loc="center",
               bbox_to_anchor=(0.72, 0.12), fontsize=7, handlelength=1.2, columnspacing=0.9)
    ref = groups["empirical_independent_DNA_measurement_provenance"].astype(str)
    provenance = "not established" if ref.str.contains("not established", case=False).all() else "see source rows"
    fig.text(0.72, 0.025,
             f"Physical measurement independence: {provenance}.",
             ha="center", fontsize=7, color=GREY)
    for ax in (ax1,):
        clean_axis(ax)
    fig.suptitle("")
    save_figure(fig, "figure1_collection_testability")


def figure2(data: dict[str, Any]) -> None:
    groups = data["brca2_groups"].copy()
    summary = data["brca2_summary"]
    history = data["brca2_covariance_history"]
    groups["primary_concordant_bool"] = bool_series(groups["primary_concordant"])
    groups["tie_any"] = bool_series(groups["Hap1_any_exact_score_tie"]) | bool_series(groups["mES_any_exact_score_tie"])
    groups["n_variants"] = pd.to_numeric(groups["n_variants"], errors="raise")
    groups["primary_d_Hap1_z"] = pd.to_numeric(groups["primary_d_Hap1_z"], errors="raise")
    groups["primary_d_mES_z"] = pd.to_numeric(groups["primary_d_mES_z"], errors="raise")
    fig, axs = plt.subplots(1, 3, figsize=(7.2, 3.15), gridspec_kw={"width_ratios": [1.35, 1.05, 0.95]},
                            constrained_layout=True)
    ax0, ax1, ax2 = axs
    for ax, letter in zip(axs, "abc"):
        label_panel(ax, letter)

    for nclass, color, marker, name in [("two", BLUE, "o", "2 variants"), ("multi", ORANGE, "^", "≥3 variants")]:
        mask = groups["n_variants"].eq(2) if nclass == "two" else groups["n_variants"].ge(3)
        concordant = mask & groups["primary_concordant_bool"] & ~groups["tie_any"]
        discordant = mask & ~groups["primary_concordant_bool"] & ~groups["tie_any"]
        tied = mask & groups["tie_any"]
        ax0.scatter(groups.loc[concordant, "primary_d_Hap1_z"], groups.loc[concordant, "primary_d_mES_z"],
                    s=13, c=color, marker=marker, alpha=0.70, linewidths=0, label=f"{name}; concordant")
        ax0.scatter(groups.loc[discordant, "primary_d_Hap1_z"], groups.loc[discordant, "primary_d_mES_z"],
                    s=14, c=GREY, marker=marker, alpha=0.75, linewidths=0, label=f"{name}; discordant")
        if tied.any():
            ax0.scatter(groups.loc[tied, "primary_d_Hap1_z"], groups.loc[tied, "primary_d_mES_z"],
                        s=28, facecolors="none", edgecolors=ORANGE if nclass == "multi" else BLUE,
                        marker=marker, linewidths=0.8, label=f"{name}; group has exact score tie")
    ax0.axhline(0, color=LIGHT, linewidth=0.6, zorder=0)
    ax0.axvline(0, color=LIGHT, linewidth=0.6, zorder=0)
    ax0.set_xlabel("HAP1 contrast (z-score units)")
    ax0.set_ylabel("mES contrast (z-score units)")
    ax0.legend(frameon=False, fontsize=7, loc="best", handletextpad=0.3, labelspacing=0.25)
    clean_axis(ax0)

    estimators = [
        ("primary_Hap1_extremes", "HAP1 extremes"),
        ("reverse_mES_extremes", "mES extremes"),
        ("fixed_first_two_HGVS_keys", "fixed first two"),
    ]
    markers = {"all": "o", "n2": "s", "n3plus": "^"}
    colors_by_size = {"all": BLUE, "n2": GREEN, "n3plus": ORANGE}
    offsets = {"all": -0.12, "n2": 0.0, "n3plus": 0.12}
    for xi, (est, _) in enumerate(estimators):
        for size in ("all", "n2", "n3plus"):
            row = summary["estimators"][est][size]
            interval = row["position_cluster_bootstrap_95"]["ci95_percentile"]
            r = float(row["R"])
            lo, hi = [float(x) for x in interval]
            ax1.errorbar(xi + offsets[size], r * 100, yerr=[[max(0, (r - lo) * 100)], [max(0, (hi - r) * 100)]],
                         fmt=markers[size], ms=3.8, color=colors_by_size[size], capsize=1.6, linewidth=0.75,
                         label=f"{size} (n={int(row['groups'])})" if xi == 0 else None, zorder=3)
    ax1.set_xticks(range(len(estimators)), [row[1] for row in estimators], rotation=20, ha="right")
    ax1.set_ylabel("sign-concordant groups, R (%)")
    ax1.set_ylim(0, 100)
    ax1.grid(axis="y", color=LIGHT, linewidth=0.45)
    ax1.set_axisbelow(True)
    ax1.legend(frameon=False, fontsize=7, loc="upper left", handletextpad=0.3, labelspacing=0.3)
    clean_axis(ax1)

    phi = number(history["phi_A_primary_preserved"], "primary Phi")
    phi_ci = [number(x, "primary Phi CI") for x in history["phi_A_primary_ci_preserved"]]
    rank_phi = number(history["rank_scale_sensitivity"]["phi_rank_scale"], "rank-scale Phi sensitivity")
    class_phi = number(history["class_control_missense_phi_sensitivity_only"], "class-control Phi sensitivity")
    ylabels = ["primary Phi", "rank-scale sensitivity", "class-control sensitivity"]
    vals = [phi, rank_phi, class_phi]
    ax2.errorbar(phi * 100, 0, xerr=[[max(0, (phi - phi_ci[0]) * 100)], [max(0, (phi_ci[1] - phi) * 100)]],
                 fmt="o", ms=4.2, color=BLUE, capsize=2, linewidth=0.9, label="original group bootstrap interval")
    ax2.scatter(rank_phi * 100, 1, s=18, marker="s", color=GREEN, label="rank units")
    ax2.scatter(class_phi * 100, 2, s=18, marker="^", color=ORANGE, label="alternate class control")
    ax2.set_yticks(range(3), ylabels)
    ax2.invert_yaxis()
    ax2.set_xlim(0, 105)
    ax2.set_xlabel("Phi (%)")
    ax2.grid(axis="x", color=LIGHT, linewidth=0.45)
    ax2.set_axisbelow(True)
    clean_axis(ax2)
    save_figure(fig, "figure2_brca2_paired_contrasts")


def figure3(data: dict[str, Any]) -> None:
    metrics = data["msh2_metrics"].copy()
    attrition = data["msh2_attrition"].copy()
    metrics["held_rep"] = pd.to_numeric(metrics["held_rep"], errors="raise").astype(int)
    attrition["held_rep"] = pd.to_numeric(attrition["held_rep"], errors="raise").astype(int)
    fig = plt.figure(figsize=(7.2, 4.10), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, height_ratios=[0.95, 1.05], width_ratios=[1.12, 1.0])
    axv = fig.add_subplot(gs[0, 0])
    axg = fig.add_subplot(gs[1, 0])
    axgain = fig.add_subplot(gs[0, 1])
    axerr = fig.add_subplot(gs[1, 1])
    label_panel(axv, "a")
    label_panel(axgain, "b")
    label_panel(axerr, "c")

    stages_v = ["candidate", "any source read", "source gates", "target complete"]
    cols_v = ["candidate_missense_variants", "observed_nonzero_any_source_condition", "source_eligible_variants", "common_target_complete_variants"]
    rep_colors = {1: BLUE, 2: GREEN, 3: ORANGE}
    x = np.arange(len(stages_v))
    width = 0.22
    for offset, (_, row) in enumerate(attrition.sort_values("held_rep").iterrows()):
        held = int(row["held_rep"])
        values = [number(row[c], c) for c in cols_v]
        axv.bar(x + (offset - 1) * width, values, width=width, color=rep_colors[held], label=f"held R{held}")
    axv.set_xticks(x, stages_v, rotation=18, ha="right")
    axv.set_ylabel("DNA variants")
    axv.legend(frameon=False, fontsize=7, ncol=3, loc="upper right", handlelength=1)
    axv.grid(axis="y", color=LIGHT, linewidth=0.4)
    axv.set_axisbelow(True)
    clean_axis(axv)

    stages_g = ["candidate AA groups", "source-eligible ≥2 DNA", "target-complete"]
    cols_g = ["candidate_amino_acid_groups", "eligible_groups_multidna", "common_target_complete_groups"]
    x = np.arange(len(stages_g))
    for offset, (_, row) in enumerate(attrition.sort_values("held_rep").iterrows()):
        held = int(row["held_rep"])
        values = [number(row[c], c) for c in cols_g]
        axg.bar(x + (offset - 1) * width, values, width=width, color=rep_colors[held])
    axg.set_xticks(x, stages_g, rotation=18, ha="right")
    axg.set_ylabel("amino-acid groups")
    axg.grid(axis="y", color=LIGHT, linewidth=0.4)
    axg.set_axisbelow(True)
    clean_axis(axg)

    endpoint_order = ["6TG_mock", "mock_P0", "DB6_P0"]
    endpoint_labels = {"6TG_mock": "6-TG / mock", "mock_P0": "mock / P0", "DB6_P0": "6-TG / P0"}
    endpoint_colors = {"6TG_mock": BLUE, "mock_P0": GREEN, "DB6_P0": ORANGE}
    for yi, endpoint in enumerate(endpoint_order):
        sub = metrics[metrics.endpoint.astype(str).eq(endpoint)].sort_values("held_rep")
        for _, row in sub.iterrows():
            val = number(row["relative_rmse_gain"], f"{endpoint} relative RMSE gain") * 100
            lo = hi = None
            if endpoint == "6TG_mock":
                ci = csv_list(row.get("relative_rmse_gain_95pct_ci"), "6-TG/mock gain CI")
                if ci:
                    lo, hi = ci
            held = int(row["held_rep"])
            yj = yi + (held - 2) * 0.13
            if lo is None:
                axgain.scatter(val, yj, s=16, color=endpoint_colors[endpoint], marker="o", zorder=3)
            else:
                axgain.errorbar(val, yj, xerr=[[max(0, (val/100 - lo) * 100)], [max(0, (hi - val/100) * 100)]],
                                fmt="o", ms=3.8, color=endpoint_colors[endpoint], capsize=1.7, linewidth=0.8, zorder=3)
    axgain.axvline(0, color=GREY, linewidth=0.7, linestyle=(0, (2, 2)))
    axgain.set_yticks(range(3), [endpoint_labels[e] for e in endpoint_order])
    axgain.set_xlabel("relative RMSE gain (%)")
    axgain.grid(axis="x", color=LIGHT, linewidth=0.4)
    axgain.set_axisbelow(True)
    clean_axis(axgain)

    # Absolute-error scales are reported directly; each linked pair is one held replicate.
    error_axes = [axerr]
    metric_specs = [("rmse_aa", "rmse_full", "RMSE"), ("mae_aa_group_equal", "mae_full_group_equal", "MAE")]
    # Compact endpoint-by-error small multiples inside the bottom-right panel.
    axerr.remove()
    inner = gs[1, 1].subgridspec(2, 3, hspace=0.38, wspace=0.30)
    axes = []
    for r, (a_col, f_col, metric_label) in enumerate(metric_specs):
        for c, endpoint in enumerate(endpoint_order):
            ax = fig.add_subplot(inner[r, c])
            axes.append(ax)
            sub = metrics[metrics.endpoint.astype(str).eq(endpoint)].sort_values("held_rep")
            for _, row in sub.iterrows():
                held = int(row["held_rep"])
                yvals = [number(row[a_col], a_col), number(row[f_col], f_col)]
                ax.plot([0, 1], yvals, marker="o", markersize=2.8, linewidth=0.7,
                        color=rep_colors[held], alpha=0.9)
            if c == 0:
                ax.set_ylabel(metric_label, fontsize=7)
            if r == 0:
                ax.set_title(endpoint_labels[endpoint], fontsize=7, pad=2)
            if r == 1:
                ax.set_xticks([0, 1], ["AA", "DNA"], rotation=25)
            else:
                ax.set_xticks([0, 1], [])
            ax.tick_params(length=1.5, pad=1)
            ax.grid(axis="y", color=LIGHT, linewidth=0.35)
            ax.set_axisbelow(True)
            ax.spines["left"].set_color(GREY)
            ax.spines["bottom"].set_color(GREY)
    if axes:
        axes[0].text(-0.55, 1.28, "c", transform=axes[0].transAxes, fontsize=8, fontweight="bold", va="bottom")
    save_figure(fig, "figure3_msh2_controlled_endpoints")


def figure4(data: dict[str, Any]) -> None:
    null = data["selection_null"].copy()
    ratio = pd.to_numeric(null["tau_over_sigma"], errors="raise")
    null0 = null[np.isclose(ratio, 0.0)].copy()
    wanted = ["R_all", "R_A", "R_B", "R_C"]
    null0 = null0.set_index("pipeline").reindex(wanted).reset_index()
    if null0["mean_R"].isna().any():
        fail("selection-null file lacks the prespecified tau/sigma=0 pipelines R_all/R_A/R_B/R_C")
    real = data["real_workflow"].copy()
    fig = plt.figure(figsize=(7.2, 3.25), constrained_layout=True)
    gs = fig.add_gridspec(1, 3, width_ratios=[0.9, 1.1, 1.15])
    ax0, ax1, ax2 = [fig.add_subplot(gs[0, i]) for i in range(3)]
    for ax, letter in zip((ax0, ax1, ax2), "abc"):
        label_panel(ax, letter)

    lim = 2.6
    grid_axis = np.linspace(-lim, lim, 401)
    d1_grid, d2_grid = np.meshgrid(grid_axis, grid_axis)
    accepted = np.abs(d1_grid) > 1.96 * np.sqrt(2 * np.maximum((d1_grid - d2_grid) ** 2 / 4, 1e-12))
    ax0.contourf(d1_grid, d2_grid, accepted.astype(int), levels=[0.5, 1.5],
                 colors=[PALE_GREEN], alpha=0.95, zorder=0)
    ax0.contour(d1_grid, d2_grid, accepted.astype(float), levels=[0.5], colors=[GREEN], linewidths=0.7)
    ax0.plot(grid_axis, grid_axis, color=BLUE, linewidth=0.7, linestyle=(0, (2, 2)))
    ax0.axhline(0, color=GREY, linewidth=0.65)
    ax0.axvline(0, color=GREY, linewidth=0.65)
    ax0.set_xlim(-lim, lim)
    ax0.set_ylim(-lim, lim)
    ax0.set_aspect("equal", adjustable="box")
    ax0.set_xlabel("D1")
    ax0.set_ylabel("D2")
    ax0.text(1.75, 1.75, "accepted", fontsize=7, color=GREEN, ha="center")
    ax0.text(-1.75, -1.75, "accepted", fontsize=7, color=GREEN, ha="center")
    clean_axis(ax0)

    y = np.arange(len(null0))
    mean = pd.to_numeric(null0["mean_R"], errors="raise").to_numpy(float) * 100
    low = pd.to_numeric(null0["ci_low"], errors="raise").to_numpy(float) * 100
    high = pd.to_numeric(null0["ci_high"], errors="raise").to_numpy(float) * 100
    colors = [GREY, BLUE, GREEN, ORANGE]
    nsel_labels = []
    for i in range(len(null0)):
        ax1.plot([low[i], high[i]], [y[i], y[i]], color=colors[i], linewidth=1.5, solid_capstyle="round")
        ax1.scatter(mean[i], y[i], color=colors[i], s=18, zorder=3)
        nsel = number(null0.iloc[i]["mean_selected"], "null mean selected")
        label = ["all", "rule A", "rule B", "rule C"][i]
        nsel_labels.append(f"{label} (mean n={nsel:.1f})")
    ax1.set_yticks(y, nsel_labels)
    ax1.invert_yaxis()
    ax1.set_xlim(0, 108)
    ax1.set_xlabel("mean dataset-level R (%)")
    ax1.grid(axis="x", color=LIGHT, linewidth=0.4)
    ax1.set_axisbelow(True)
    clean_axis(ax1)

    # This panel preserves each workflow's own selected population and interval availability.
    real["R"] = pd.to_numeric(real["R"], errors="raise")
    real["n_selected"] = pd.to_numeric(real["n_selected"], errors="raise")
    real["R_ci_low"] = pd.to_numeric(real["R_ci_low"], errors="coerce")
    real["R_ci_high"] = pd.to_numeric(real["R_ci_high"], errors="coerce")
    y2 = np.arange(len(real))
    for i, row in real.reset_index(drop=True).iterrows():
        r = float(row["R"])
        lo, hi = row["R_ci_low"], row["R_ci_high"]
        if pd.notna(lo) and pd.notna(hi):
            ax2.errorbar(r * 100, y2[i], xerr=[[max(0, (r - float(lo)) * 100)],
                [max(0, (float(hi) - r) * 100)]], fmt="o", ms=4, color=BLUE, capsize=1.8, linewidth=0.8)
        else:
            ax2.scatter(r * 100, y2[i], s=18, color=BLUE, marker="o")
    label_map = {
        "A_validation_dependent": "A validation-dependent",
        "B_external_lab": "B fixed (lab label)",
        "B_external_tech": "B fixed (technical label)",
        "C_position_crossfit": "C position crossfit",
    }
    labels = [f"{label_map.get(str(row.pipeline), str(row.pipeline).replace('_', ' '))} ({'mean n' if str(row.pipeline) == 'C_position_crossfit' else 'n'}={float(row.n_selected):g})"
              for row in real.itertuples()]
    ax2.set_yticks(y2, labels)
    ax2.invert_yaxis()
    ax2.set_xlim(0, 108)
    ax2.set_xlabel("workflow-specific R (%)")
    ax2.grid(axis="x", color=LIGHT, linewidth=0.4)
    ax2.set_axisbelow(True)
    clean_axis(ax2)
    save_figure(fig, "figure4_structural_selection_counterexample")


def optional_supplement(data: dict[str, Any]) -> None:
    sensitivity = data.get("msh2_sensitivity")
    if sensitivity is None:
        return
    require_columns(sensitivity, ["held_rep", "endpoint", "relative_rmse_gain", "scenario",
        "pseudocount", "depth_frequency_threshold", "count_correction"], "MSH2 sensitivity metrics")
    subset = sensitivity[sensitivity.endpoint.astype(str).eq("6TG_mock")].copy()
    if subset.empty:
        return
    fig, ax = plt.subplots(figsize=(7.2, 3.0), constrained_layout=True)
    labels = []
    for yi, (setting, group) in enumerate(subset.groupby(["scenario", "pseudocount", "depth_frequency_threshold", "count_correction"], sort=True)):
        scenario, pseudo, depth, correction = setting
        name = (f"source frequency > {float(depth):g}" if str(scenario).startswith("depth") else
                f"pseudocount {float(pseudo):g}" if str(scenario).startswith("pseudocount") else
                "uncorrected counts; primary cohort")
        labels.append(name)
        for row in group.itertuples():
            held = int(row.held_rep)
            ax.scatter(float(row.relative_rmse_gain) * 100, yi + (held - 2) * .12, s=18,
                       color=[BLUE, GREEN, ORANGE][held - 1], marker=["o", "s", "^"][held - 1],
                       label=f"held replicate {held}" if yi == 0 else None)
    ax.axvline(0, color=GREY, linewidth=0.65, linestyle=(0, (2, 2)))
    ax.set_yticks(range(len(labels)), labels)
    ax.set_xlabel("6-TG/mock relative RMSE gain (%)")
    ax.set_ylabel("prespecified sensitivity")
    ax.legend(frameon=False, fontsize=7, loc="upper left")
    clean_axis(ax)
    save_figure(fig, "figureS1_msh2_sensitivity")


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, default=DEFAULT_DATA_ROOT,
                        help="root containing the frozen data/derived tree")
    parser.add_argument("--output-dir", type=Path, default=REPO_ROOT / "results/figures",
                        help="directory for rendered main and S1 figures")
    parser.add_argument("--source-dir", type=Path, default=REPO_ROOT / "results/figure_source",
                        help="directory for exported source tables and input hash manifest")
    args = parser.parse_args(argv)
    configure_paths(args.data_root, args.output_dir, args.source_dir)
    check_inputs()
    data = load_inputs()
    OUT.mkdir(parents=True, exist_ok=True)
    source_export(data)
    figure1(data)
    figure2(data)
    figure3(data)
    figure4(data)
    optional_supplement(data)
    print("Main figure PDFs and PNG previews written to:", OUT)
    print("Source-data tables and input hashes written to:", SOURCE)


if __name__ == "__main__":
    main()


