#!/usr/bin/env python3
"""Render the display revision from frozen, verified figure inputs."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
from typing import Any

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.lines import Line2D
from matplotlib.text import Text
import numpy as np
import pandas as pd

import build_revision_figures as base


ROOT = Path(__file__).resolve().parent
TASK_ROOT = ROOT.parent
DEFAULT_DATA_ROOT = TASK_ROOT / "data" / "derived"
DEFAULT_OUT = TASK_ROOT / "results" / "figures"
DEFAULT_SOURCE = DEFAULT_OUT / "source_data"

BLUE = "#3B6F9C"
GREEN = "#5B8E75"
ORANGE = "#C7864D"
GREY = "#6D7378"
LIGHT = "#D7DEE4"
PALE_BLUE = "#E7EEF4"
PALE_GREEN = "#E7EFEA"
PALE_ORANGE = "#F5EADF"
REP_COLORS = {1: BLUE, 2: GREEN, 3: ORANGE}
REP_MARKERS = {1: "o", 2: "s", 3: "^"}
INK = "#192C3B"
FIGURE_WIDTH = 170 / 25.4

if os.environ.get("FIGURE_ARIAL_BOLD"):
    font_manager.fontManager.addfont(os.environ["FIGURE_ARIAL_BOLD"])

try:
    font_manager.findfont(font_manager.FontProperties(family="Arial", weight="bold"),
                          fallback_to_default=False)
    families = ["Arial", "DejaVu Sans"]
except ValueError:
    families = ["DejaVu Sans"]

mpl.rcParams.update({
    "font.family": families,
    "font.sans-serif": ["Arial", "DejaVu Sans"],
    "font.size": 9.5,
    "font.weight": "bold",
    "text.color": INK,
    "axes.labelsize": 10,
    "axes.labelweight": "bold",
    "axes.titlesize": 11,
    "axes.titleweight": "bold",
    "xtick.labelsize": 9,
    "ytick.labelsize": 9.5,
    "legend.fontsize": 9,
    "legend.frameon": False,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 0.7,
    "xtick.major.width": 0.7,
    "ytick.major.width": 0.7,
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "savefig.facecolor": "white",
    "figure.facecolor": "white",
})


def num(value: Any, label: str) -> float:
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"non-finite {label}: {value!r}")
    return result


def save(fig: mpl.figure.Figure, stem: str) -> None:
    DEFAULT_OUT.mkdir(parents=True, exist_ok=True)
    for label in fig.findobj(Text):
        label.set_fontweight("bold")
    fig.savefig(DEFAULT_OUT / f"{stem}.pdf")
    fig.savefig(DEFAULT_OUT / f"{stem}.png", dpi=600)
    fig.savefig(DEFAULT_OUT / f"{stem}.svg")
    plt.close(fig)


def panel(fig: mpl.figure.Figure, letter: str, x: float, y: float, title: str = "") -> None:
    fig.text(x, y, letter, fontsize=18, va="top", color=INK)
    if title:
        offset = 25 / (72 * fig.get_figwidth())
        title_y = y - 3.6 / (72 * fig.get_figheight())
        fig.text(x + offset, title_y, title, fontsize=11.5, va="top", color=INK)


def tidy(ax: plt.Axes) -> None:
    ax.spines["left"].set_color(GREY)
    ax.spines["bottom"].set_color(GREY)
    ax.tick_params(colors=INK, length=2.8, pad=3)
    ax.xaxis.label.set_color(INK)
    ax.yaxis.label.set_color(INK)


def tier_data(data: dict[str, Any]) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    historical = data["historical_reproduction"]["observed"]
    observed = data["collection_summary"]["observed"]
    counts = pd.DataFrame([
        {"measure": "manifest score sets", "value": int(historical["manifest_score_sets"]), "unit": "score sets"},
        {"measure": "analyzable score sets", "value": int(observed["historical_analyzable_score_sets"]), "unit": "score sets"},
        {"measure": "sets with multi-key groups", "value": int(observed["eligible_score_sets_with_multi_key_groups"]), "unit": "score sets"},
        {"measure": "missense rows", "value": int(historical["missense_rows"]), "unit": "score rows"},
        {"measure": "historical group–key entries", "value": int(historical["unique_keys_in_ge2_groups"]), "unit": "group–key entries"},
        {"measure": "historical multi-key groups", "value": int(historical["groups_ge2"]), "unit": "amino-acid groups"},
    ])

    group_table = data["collection_groups"]
    specs = [
        ("Accession reference", "accession_reference_group_state", "ge2", "accession_reference_groups_ge2_all", 3836, 19),
        ("Versioned request", "versioned_request_group_state", "ge2", "versioned_request_groups_ge2_all", 2679, 13),
        ("Conditional target sequence", "target_sequence_group_state", "conditional_ge2", "target_conditional_groups_ge2_non_TPK1", 1013, 7),
    ]
    denominator = int(observed["historical_multi_key_groups"])
    per_set = data["collection_per_set"]
    hero = []
    for tier, group_col, eligible_state, set_col, expected_groups, expected_sets in specs:
        n_groups = int(group_table[group_col].astype(str).eq(eligible_state).sum())
        n_sets = int(pd.to_numeric(per_set[set_col], errors="raise").gt(0).sum())
        if (n_groups, n_sets) != (expected_groups, expected_sets):
            raise ValueError(f"frozen collection tier mismatch for {tier}: {(n_groups, n_sets)}")
        hero.append({"tier": tier, "qualified_groups": n_groups, "all_historical_groups": denominator,
                     "sets_with_qualified_groups": n_sets, "eligible_state": eligible_state})
    hero_df = pd.DataFrame(hero)

    states = []
    state_specs = [
        ("Accession reference", "accession_reference_group_state",
         [("ge2", "≥2 supported"), ("one", "1 supported"), ("zero", "0 supported"), ("unknown_0_or_1", "unknown (0 or 1)")]),
        ("Versioned request", "versioned_request_group_state",
         [("ge2", "≥2 supported"), ("one", "1 supported"), ("zero", "0 supported"), ("unknown_0_or_1", "unknown (0 or 1)")]),
        ("Target conditional", "target_sequence_group_state",
         [("conditional_ge2", "≥2 target-compatible*"), ("conditional_one", "1 target-compatible*"),
          ("zero_target_compatible_under_assumption", "0 target-compatible*"),
          ("not_assessed_target_sequence", "not assessed"),
          ("tpk1_full_construct_unresolved", "TPK1 unresolved")]),
    ]
    for tier, column, expected in state_specs:
        observed_states = group_table[column].astype(str)
        all_states = list(dict.fromkeys([state for state, _ in expected] + sorted(set(observed_states))))
        labels = dict(expected)
        for state in all_states:
            states.append({"tier": tier, "state": state,
                           "state_label": labels.get(state, state.replace("_", " ")),
                           "groups": int(observed_states.eq(state).sum())})
    state_df = pd.DataFrame(states)
    if state_df.groupby("tier")["groups"].sum().ne(denominator).any():
        raise ValueError("collection state bars do not each sum to the historical group denominator")
    return counts, hero_df, state_df


def figure1(data: dict[str, Any]) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    counts, hero, _ = tier_data(data)
    fig = plt.figure(figsize=(FIGURE_WIDTH, 3.15))
    ax0 = fig.add_axes([0.11, 0.47, 0.86, 0.46])
    ax1 = fig.add_axes([0.28, 0.18, 0.69, 0.27])
    panel(fig, "A", 0.02, 0.96)
    panel(fig, "B", 0.02, 0.46)
    ax0.set_axis_off()
    for i, row in counts.iterrows():
        col, r = i % 3, i // 3
        x = 0.01 + col * 0.335
        y = 0.84 - r * 0.48
        color = [BLUE, "#376750", "#975521", INK, BLUE, "#376750"][i]
        value = f"{int(row.value):,}"
        ax0.text(x, y, value, transform=ax0.transAxes, fontsize=17, color=color, va="center")
        ax0.text(x, y - 0.165, str(row.measure), transform=ax0.transAxes, fontsize=9.2, va="center")
        ax0.text(x, y - 0.278, str(row.unit), transform=ax0.transAxes, fontsize=8.5, va="center")
    colors = [BLUE, GREEN, ORANGE]
    y = np.arange(len(hero))
    ax1.barh(y, hero["qualified_groups"], color=colors, height=0.56, edgecolor="white", linewidth=0.5)
    ax1.set_yticks(y, ["Reference compatible", "Versioned request", "Conditional target\nsequence"])
    ax1.invert_yaxis()
    ax1.set_xlim(0, 7200)
    ax1.set_xlabel("Qualified amino-acid groups (of 13,041 historical multi-key groups)")
    ax1.grid(axis="x", color=LIGHT, linewidth=0.45)
    ax1.set_axisbelow(True)
    for yi, row in hero.iterrows():
        ax1.text(int(row.qualified_groups) + 90, yi,
                 f"{int(row.qualified_groups):,} / {int(row.all_historical_groups):,}  ({int(row.sets_with_qualified_groups)} sets)",
                 va="center", ha="left", fontsize=9.3)
    tidy(ax1)
    save(fig, "figure1_collection_testability")
    return counts, hero, tier_data(data)[2]


def figure_s3(data: dict[str, Any], state_df: pd.DataFrame) -> pd.DataFrame:
    historical = data["historical_reproduction"]["observed"]
    dist = historical["group_size_distribution"]
    hist = pd.DataFrame([{"keys_per_group": int(k), "groups": int(v)} for k, v in dist.items()]).sort_values("keys_per_group")
    fig = plt.figure(figsize=(FIGURE_WIDTH, 5.40))
    ax0 = fig.add_axes([0.145, 0.69, 0.82, 0.22])
    ax1 = fig.add_axes([0.44, 0.105, 0.525, 0.425])
    panel(fig, "A", 0.02, 0.96)
    panel(fig, "B", 0.02, 0.565)
    x = np.arange(len(hist))
    bars = ax0.bar(x, hist["groups"], color=BLUE, width=0.68, edgecolor="white", linewidth=0.5)
    ax0.set_xticks(x, hist["keys_per_group"].astype(int).astype(str))
    ax0.set_xlabel("Historical DNA keys per amino-acid group")
    ax0.set_ylabel("Groups")
    ax0.set_ylim(0, hist.groups.max() * 1.24)
    ax0.bar_label(bars, labels=[f"{int(v):,}" for v in hist.groups], padding=3, fontsize=9.5)
    ax0.grid(axis="y", color=LIGHT, linewidth=0.45)
    ax0.set_axisbelow(True)
    tidy(ax0)

    tier_order = ["Accession reference", "Versioned request", "Target conditional"]
    prefix = {"Accession reference": "Reference", "Versioned request": "Versioned request", "Target conditional": "Target conditional"}
    state_colors = {"ge2": BLUE, "one": GREEN, "zero": ORANGE, "unknown_0_or_1": GREY,
                    "conditional_ge2": BLUE, "conditional_one": GREEN,
                    "zero_target_compatible_under_assumption": ORANGE,
                    "not_assessed_target_sequence": GREY, "tpk1_full_construct_unresolved": LIGHT}
    rows = []
    for tier in tier_order:
        sub = state_df[state_df.tier.eq(tier)]
        rows.extend(sub.to_dict("records"))
    plot_df = pd.DataFrame(rows).reset_index(drop=True)
    ypos = np.arange(len(plot_df))
    ax1.barh(ypos, plot_df.groups, color=[state_colors.get(str(s), GREY) for s in plot_df.state],
             height=0.68, edgecolor="white", linewidth=0.45)
    labels = [f"{prefix[t]} — {label}" for t, label in zip(plot_df.tier, plot_df.state_label)]
    ax1.set_yticks(ypos, labels)
    ax1.invert_yaxis()
    ax1.set_xlabel("Groups (denominator 13,041 for each tier)")
    ax1.grid(axis="x", color=LIGHT, linewidth=0.45)
    ax1.set_axisbelow(True)
    ax1.set_xlim(0, 14200)
    ax1.axhspan(3.5, 7.5, color="#F1F4F7", zorder=0)
    for yi, n in enumerate(plot_df.groups.astype(int)):
        ax1.text(max(n + 100, 130), yi, f"{n:,}", va="center", ha="left", fontsize=9.2)
    tidy(ax1)
    save(fig, "Supplementary_Figure_S3_collection_detail")
    return hist


def figure2(data: dict[str, Any]) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    groups = data["brca2_groups"].copy()
    summary = data["brca2_summary"]
    history = data["brca2_covariance_history"]
    groups["n_variants"] = pd.to_numeric(groups["n_variants"], errors="raise").astype(int)
    groups["primary_d_Hap1_z"] = pd.to_numeric(groups["primary_d_Hap1_z"], errors="raise")
    groups["primary_d_mES_z"] = pd.to_numeric(groups["primary_d_mES_z"], errors="raise")
    groups["primary_concordant_bool"] = base.bool_series(groups["primary_concordant"])
    if len(groups) != 388:
        raise ValueError(f"expected 388 BRCA2 groups, found {len(groups)}")

    fig = plt.figure(figsize=(FIGURE_WIDTH, 4.88))
    ax0 = fig.add_axes([0.115, 0.49, 0.36, 0.38])
    ax1 = fig.add_axes([0.635, 0.49, 0.335, 0.38])
    ax2 = fig.add_axes([0.247, 0.093, 0.497, 0.155])
    panel(fig, "A", 0.02, 0.972, "388 paired groups")
    panel(fig, "B", 0.532, 0.972, "Direction checks")
    panel(fig, "C", 0.02, 0.302, "Scale sensitivity")
    for nclass, marker, label in [(2, "o", "2 variants"), (3, "^", "≥3 variants")]:
        mask = groups.n_variants.eq(2) if nclass == 2 else groups.n_variants.ge(3)
        colors = np.where(groups.primary_concordant_bool[mask], BLUE, ORANGE)
        ax0.scatter(groups.loc[mask, "primary_d_Hap1_z"], groups.loc[mask, "primary_d_mES_z"],
                    s=13, c=colors, marker=marker, alpha=0.75, linewidths=0, label=label)
    ax0.axhline(0, color=LIGHT, linewidth=0.65, zorder=0)
    ax0.axvline(0, color=LIGHT, linewidth=0.65, zorder=0)
    ax0.set_xlabel("HAP1 contrast (z-score units)")
    ax0.set_ylabel("mES contrast (z-score units)")
    scatter_key = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor=BLUE, markeredgecolor=BLUE,
               markersize=4.2, label="Concordant"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor=ORANGE, markeredgecolor=ORANGE,
               markersize=4.2, label="Discordant"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor=GREY, markeredgecolor=GREY,
               markersize=4.2, label="2 variants"),
        Line2D([0], [0], marker="^", color="none", markerfacecolor=GREY, markeredgecolor=GREY,
               markersize=4.2, label="≥3 variants"),
    ]
    ax0.legend(handles=scatter_key, loc="upper center", bbox_to_anchor=(0.5, -0.245),
               ncol=2, handletextpad=0.35, labelspacing=0.22, columnspacing=1.05)
    tidy(ax0)

    estimators = [
        ("primary_Hap1_extremes", "HAP1\nextremes"),
        ("reverse_mES_extremes", "mES\nextremes"),
        ("fixed_first_two_HGVS_keys", "Fixed\nfirst two"),
    ]
    strata = [("all", BLUE, "o", -0.17, "All"), ("n2", GREEN, "s", 0.0, "2 variants"),
              ("n3plus", ORANGE, "^", 0.17, "≥3 variants")]
    direction_rows = []
    for xi, (method, short_label) in enumerate(estimators):
        for size, color, marker, offset, legend in strata:
            row = summary["estimators"][method][size]
            r = num(row["R"], "BRCA2 R") * 100
            lo, hi = [num(v, "BRCA2 interval") * 100 for v in row["position_cluster_bootstrap_95"]["ci95_percentile"]]
            direction_rows.append({"method": method, "stratum": size, "groups": int(row["groups"]),
                                   "R_percent": r, "ci_low_percent": lo, "ci_high_percent": hi})
            ax1.errorbar(r, xi + offset, xerr=[[r - lo], [hi - r]], fmt=marker, ms=4.0,
                         color=color, capsize=1.7, linewidth=0.75, label=legend if xi == 0 else None, zorder=3)
    ax1.axvline(50, color=GREY, linewidth=0.85, linestyle=(0, (3, 2)), zorder=2)
    ax1.set_yticks(range(len(estimators)), [item[1] for item in estimators])
    ax1.invert_yaxis()
    ax1.set_xlim(30, 70)
    ax1.set_xticks([30, 40, 50, 60, 70])
    ax1.set_xlabel("Sign-concordant groups, R (%)")
    ax1.legend(loc="upper center", bbox_to_anchor=(0.43, -0.245),
               ncol=3, handlelength=1.1, handletextpad=0.35, columnspacing=0.8)
    ax1.grid(axis="x", color=LIGHT, linewidth=0.45)
    ax1.set_axisbelow(True)
    tidy(ax1)

    phi = num(history["phi_A_primary_preserved"], "primary Phi") * 100
    phi_ci = [num(v, "primary Phi CI") * 100 for v in history["phi_A_primary_ci_preserved"]]
    rank = num(history["rank_scale_sensitivity"]["phi_rank_scale"], "rank-scale Phi") * 100
    alt = num(history["class_control_missense_phi_sensitivity_only"], "class-control Phi") * 100
    phi_rows = [
        {"estimate": "Primary", "phi_percent": phi, "ci_low_percent": phi_ci[0], "ci_high_percent": phi_ci[1]},
        {"estimate": "Rank-scale sensitivity", "phi_percent": rank, "ci_low_percent": np.nan, "ci_high_percent": np.nan},
        {"estimate": "Class-control sensitivity", "phi_percent": alt, "ci_low_percent": np.nan, "ci_high_percent": np.nan},
    ]
    vals = [phi, rank, alt]
    for yi, value in enumerate(vals):
        if yi == 0:
            ax2.errorbar(value, yi, xerr=[[value - phi_ci[0]], [phi_ci[1] - value]], fmt="o", ms=4.0,
                         color=BLUE, capsize=2, linewidth=0.85, zorder=3)
        else:
            ax2.scatter(value, yi, s=22, color=GREEN if yi == 1 else ORANGE,
                        marker="s" if yi == 1 else "^", zorder=3)
    ax2.set_yticks(range(3), ["Primary", "Rank scale", "Class control"])
    ax2.invert_yaxis()
    ax2.set_ylim(2.55, -0.35)
    ax2.set_xlim(65, 105)
    ax2.set_xticks([70, 80, 90, 100])
    ax2.set_xlabel("Phi (%)")
    labels = [f"{phi:.2f} ({phi_ci[0]:.2f}–{phi_ci[1]:.2f})", f"{rank:.3f}", f"{alt:.2f}"]
    for yi, label in enumerate(labels):
        ax2.text(1.045, yi, label, transform=ax2.get_yaxis_transform(),
                 fontsize=9.2, va="center", ha="left")
    ax2.grid(axis="x", color=LIGHT, linewidth=0.45)
    ax2.set_axisbelow(True)
    tidy(ax2)
    save(fig, "figure2_brca2_paired_contrasts")
    return groups, pd.DataFrame(direction_rows), pd.DataFrame(phi_rows)


def figure3(data: dict[str, Any]) -> tuple[pd.DataFrame, pd.DataFrame]:
    metrics = data["msh2_metrics"].copy()
    attrition = data["msh2_attrition"].copy()
    metrics["held_rep"] = pd.to_numeric(metrics["held_rep"], errors="raise").astype(int)
    metrics["relative_rmse_gain"] = pd.to_numeric(metrics["relative_rmse_gain"], errors="raise")
    endpoint_order = ["6TG_mock", "mock_P0", "DB6_P0"]
    endpoint_labels = {"6TG_mock": "6-TG / mock", "mock_P0": "mock / P0", "DB6_P0": "6-TG / P0"}
    fig = plt.figure(figsize=(FIGURE_WIDTH, 3.13))
    ax0 = fig.add_axes([0.112, 0.235, 0.346, 0.605])
    ax1 = fig.add_axes([0.632, 0.235, 0.333, 0.605])
    panel(fig, "A", 0.02, 0.965, "Primary endpoint")
    panel(fig, "B", 0.52, 0.965, "All endpoint gains")

    primary = metrics[metrics.endpoint.astype(str).eq("6TG_mock")].sort_values("held_rep")
    if len(primary) != 3:
        raise ValueError("primary MSH2 endpoint must include all three held replicates")
    for row in primary.itertuples():
        held = int(row.held_rep)
        values = [num(row.rmse_aa, "AA RMSE"), num(row.rmse_full, "DNA RMSE")]
        ax0.plot([0, 1], values, color=REP_COLORS[held], linewidth=1.0, marker=REP_MARKERS[held],
                 markersize=4.5, label=f"Held R{held}", zorder=3)
        for x, value in enumerate(values):
            offset = (-7, 0) if x == 0 else (7, 0)
            ax0.annotate(f"{value:.2f}", (x, value), xytext=offset,
                         textcoords="offset points", ha="right" if x == 0 else "left",
                         va="center", fontsize=9.4)
    ax0.set_xticks([0, 1], ["AA model", "DNA model"])
    ax0.set_ylabel("6-TG / mock RMSE")
    ax0.set_xlim(-0.38, 1.38)
    ax0.set_ylim(1.97, 2.35)
    ax0.grid(axis="y", color=LIGHT, linewidth=0.45)
    ax0.set_axisbelow(True)
    tidy(ax0)

    offsets = {1: -0.15, 2: 0.0, 3: 0.15}
    row_positions = {endpoint: i for i, endpoint in enumerate(endpoint_order)}
    for endpoint, yi in row_positions.items():
        sub = metrics[metrics.endpoint.astype(str).eq(endpoint)].sort_values("held_rep")
        if len(sub) != 3:
            raise ValueError(f"MSH2 endpoint {endpoint} does not include all three held replicates")
        for row in sub.itertuples():
            held = int(row.held_rep)
            value = num(row.relative_rmse_gain, f"{endpoint} relative RMSE gain") * 100
            yy = yi + offsets[held]
            ci = base.csv_list(getattr(row, "relative_rmse_gain_95pct_ci", None), f"{endpoint} primary CI") if endpoint == "6TG_mock" else None
            if ci:
                lo, hi = [num(x, "gain CI") * 100 for x in ci]
                ax1.errorbar(value, yy, xerr=[[max(0, value - lo)], [max(0, hi - value)]],
                             fmt=REP_MARKERS[held], ms=4.0, color=REP_COLORS[held], capsize=1.7,
                             linewidth=0.8, zorder=3)
            else:
                ax1.scatter(value, yy, s=19, color=REP_COLORS[held], marker=REP_MARKERS[held], zorder=3)
    ax1.axhspan(-0.42, 0.42, color=PALE_BLUE, zorder=0)
    ax1.axvline(0, color=GREY, linewidth=0.85, linestyle=(0, (3, 2)), zorder=2)
    ax1.set_yticks(range(3), [endpoint_labels[e] for e in endpoint_order])
    ax1.get_yticklabels()[0].set_fontweight("bold")
    ax1.invert_yaxis()
    ax1.set_xlabel("Relative RMSE gain (%)")
    ax1.set_ylim(2.48, -0.48)
    ax1.grid(axis="x", color=LIGHT, linewidth=0.45)
    ax1.set_axisbelow(True)
    handles, names = ax0.get_legend_handles_labels()
    fig.legend(handles, names, loc="lower center", bbox_to_anchor=(0.5, 0.027),
               ncol=3, handlelength=1.4, handletextpad=0.55, columnspacing=1.3)
    tidy(ax1)
    save(fig, "figure3_msh2_controlled_endpoints")
    return metrics, attrition


def figure4(data: dict[str, Any]) -> tuple[pd.DataFrame, pd.DataFrame]:
    null = data["selection_null"].copy()
    null["tau_over_sigma"] = pd.to_numeric(null.tau_over_sigma, errors="raise")
    null0 = null[np.isclose(null.tau_over_sigma, 0.0)].copy()
    wanted = ["R_all", "R_A", "R_B", "R_C"]
    null0 = null0.set_index("pipeline").reindex(wanted).reset_index()
    if null0.mean_R.isna().any():
        raise ValueError("frozen null simulation lacks tau/sigma=0 rules")
    real = data["real_workflow"].copy()
    fig = plt.figure(figsize=(FIGURE_WIDTH, 3.95))
    ax0 = fig.add_axes([0.105, 0.335, 0.345, 0.54])
    ax1 = fig.add_axes([0.638, 0.635, 0.33, 0.23])
    ax2 = fig.add_axes([0.638, 0.169, 0.33, 0.23])
    panel(fig, "A", 0.02, 0.97, "Analytic two-variant rule")
    panel(fig, "B", 0.493, 0.97, "Null simulation (τ/σ = 0)")
    panel(fig, "C", 0.493, 0.505, "Observed workflows")

    lim = 2.6
    grid = np.linspace(-lim, lim, 401)
    d1, d2 = np.meshgrid(grid, grid)
    accepted = np.abs(d1) > 1.96 * np.sqrt(2 * np.maximum((d1 - d2) ** 2 / 4, 1e-12))
    ax0.contourf(d1, d2, accepted.astype(int), levels=[0.5, 1.5], colors=[PALE_GREEN], alpha=0.95, zorder=0)
    ax0.contour(d1, d2, accepted.astype(float), levels=[0.5], colors=[GREEN], linewidths=0.75)
    ax0.plot(grid, grid, color=BLUE, linewidth=0.8, linestyle=(0, (3, 2)), label="D1 = D2")
    ax0.axhline(0, color=GREY, linewidth=0.55)
    ax0.axvline(0, color=GREY, linewidth=0.55)
    ax0.set_xlim(-lim, lim)
    ax0.set_ylim(-lim, lim)
    ax0.set_aspect("equal", adjustable="box")
    ax0.set_xlabel("D1")
    ax0.set_ylabel("D2")
    fig.text(0.105, 0.197, "Shading: accepted region\nNo observations plotted",
             fontsize=9.3, va="top", linespacing=1.3)
    ax0.legend(loc="upper left", handlelength=1.4)
    tidy(ax0)

    names = ["All", "A", "B", "C"]
    colors = [GREY, BLUE, GREEN, ORANGE]
    y = np.arange(4)
    for i, row in null0.iterrows():
        mean = num(row.mean_R, "null mean R") * 100
        lo, hi = num(row.ci_low, "null simulation interval") * 100, num(row.ci_high, "null simulation interval") * 100
        if lo < 0 or hi > 110:
            raise ValueError(f"null simulation interval falls outside display axis: {lo}–{hi}")
        ax1.plot([lo, hi], [y[i], y[i]], color=colors[i], linewidth=1.6, solid_capstyle="round")
        ax1.scatter(mean, y[i], color=colors[i], s=20, zorder=3)
        selected = num(row.mean_selected, "mean selected groups")
        names[i] = f"{names[i]} (n̄={selected:.1f})"
    ax1.axvline(50, color=GREY, linewidth=0.8, linestyle=(0, (3, 2)), zorder=2)
    ax1.set_yticks(y, names)
    ax1.invert_yaxis()
    ax1.set_xlim(0, 110)
    ax1.set_xticks([0, 25, 50, 75, 100])
    ax1.set_xlabel("Mean R across null simulations (%)")
    ax1.grid(axis="x", color=LIGHT, linewidth=0.45)
    ax1.set_axisbelow(True)
    tidy(ax1)

    real["R"] = pd.to_numeric(real.R, errors="raise")
    real["n_selected"] = pd.to_numeric(real.n_selected, errors="raise")
    real["R_ci_low"] = pd.to_numeric(real.R_ci_low, errors="coerce")
    real["R_ci_high"] = pd.to_numeric(real.R_ci_high, errors="coerce")
    label_map = {
        "A_validation_dependent": "A",
        "B_external_lab": "B-lab",
        "B_external_tech": "B-tech",
        "C_position_crossfit": "C",
    }
    y2 = np.arange(len(real))
    labels = []
    workflow_colors = [BLUE, GREEN, ORANGE, GREY]
    for i, row in real.reset_index(drop=True).iterrows():
        r = float(row.R) * 100
        lo, hi = row.R_ci_low, row.R_ci_high
        if pd.notna(lo) and pd.notna(hi):
            if float(lo) * 100 < 0 or float(hi) * 100 > 108:
                raise ValueError(f"observed-workflow interval falls outside display axis: {lo}–{hi}")
            ax2.errorbar(r, y2[i], xerr=[[max(0, (float(row.R) - float(lo)) * 100)],
                                         [max(0, (float(hi) - float(row.R)) * 100)]],
                         fmt="o", ms=4.2, color=workflow_colors[i], capsize=1.8, linewidth=0.8, zorder=3)
        else:
            ax2.scatter(r, y2[i], s=19, color=workflow_colors[i], marker="o", zorder=3)
        n_label = f"mean {float(row.n_selected):.1f}" if str(row.pipeline) == "C_position_crossfit" else f"n={float(row.n_selected):g}"
        labels.append(f"{label_map.get(str(row.pipeline), str(row.pipeline))} ({n_label})")
    ax2.set_yticks(y2, labels)
    ax2.invert_yaxis()
    ax2.set_xlim(0, 108)
    ax2.set_xlabel("Workflow-specific R (%)")
    ax2.grid(axis="x", color=LIGHT, linewidth=0.45)
    ax2.set_axisbelow(True)
    tidy(ax2)
    save(fig, "figure4_structural_selection_counterexample")
    return null0, real


def figure_s1(data: dict[str, Any]) -> tuple[pd.DataFrame, pd.DataFrame]:
    sensitivity = data["msh2_sensitivity"].copy()
    if sensitivity is None:
        raise ValueError("frozen MSH2 sensitivity table is missing")
    subset = sensitivity[sensitivity.endpoint.astype(str).eq("6TG_mock")].copy()
    group_cols = ["scenario", "pseudocount", "depth_frequency_threshold", "count_correction"]
    settings = subset[group_cols].drop_duplicates().copy()
    if len(settings) != 5:
        raise ValueError(f"expected five MSH2 sensitivity settings, found {len(settings)}")

    def setting_label(row: Any) -> str:
        scenario = str(row.scenario)
        if scenario.startswith("depth"):
            depth = float(row.depth_frequency_threshold)
            denom = round(1.0 / depth)
            return f"Source frequency >1/{denom:,}"
        if scenario.startswith("pseudocount"):
            return f"Pseudocount {float(row.pseudocount):g}"
        return "Uncorrected counts"

    labeled = [(setting_label(row), row) for row in settings.itertuples(index=False)]
    order_key = {"Pseudocount 0.1": 0, "Pseudocount 1": 1, "Source frequency >1/400,000": 2,
                 "Source frequency >1/100,000": 3, "Uncorrected counts": 4}
    labeled.sort(key=lambda item: order_key.get(item[0], 99))
    metrics = data["msh2_metrics"].copy()
    metrics["held_rep"] = pd.to_numeric(metrics.held_rep, errors="raise").astype(int)
    primary = metrics[metrics.endpoint.astype(str).eq("6TG_mock")].sort_values("held_rep")
    fig = plt.figure(figsize=(FIGURE_WIDTH, 3.20))
    ax = fig.add_axes([0.324, 0.177, 0.643, 0.68])
    labels = [name for name, _ in labeled] + ["Primary reference"]
    for yi, (name, setting) in enumerate(labeled):
        mask = np.ones(len(subset), dtype=bool)
        for col in group_cols:
            mask &= subset[col].astype(str).eq(str(getattr(setting, col)))
        rows = subset[mask].sort_values("held_rep")
        if len(rows) != 3:
            raise ValueError(f"sensitivity {name} does not include three held replicates")
        for row in rows.itertuples():
            rep = int(row.held_rep)
            ax.scatter(num(row.relative_rmse_gain, "sensitivity RMSE gain") * 100,
                       yi + {1: -0.14, 2: 0.0, 3: 0.14}[rep], s=20,
                       color=REP_COLORS[rep], marker=REP_MARKERS[rep], zorder=3)
    yi = len(labeled)
    for row in primary.itertuples():
        rep = int(row.held_rep)
        value = num(row.relative_rmse_gain, "primary RMSE gain") * 100
        ci = base.csv_list(getattr(row, "relative_rmse_gain_95pct_ci", None), "primary gain CI")
        if not ci:
            raise ValueError(f"primary held replicate {rep} lacks its frozen CI")
        lo, hi = [num(v, "primary CI bound") * 100 for v in ci]
        ax.errorbar(value, yi + {1: -0.14, 2: 0.0, 3: 0.14}[rep], xerr=[[max(0, value - lo)], [max(0, hi - value)]],
                    fmt=REP_MARKERS[rep], ms=4.1, color=REP_COLORS[rep], capsize=1.8, linewidth=0.8, zorder=3)
    ax.axhline(yi - 0.5, color=LIGHT, linewidth=0.8)
    ax.axvline(0, color=GREY, linewidth=0.8, linestyle=(0, (3, 2)), zorder=2)
    ax.set_yticks(range(len(labels)), labels)
    ax.invert_yaxis()
    ax.set_xlabel("6-TG / mock relative RMSE gain (%)")
    ax.set_xlim(-1, 8)
    ax.grid(axis="x", color=LIGHT, linewidth=0.45)
    ax.set_axisbelow(True)
    handles = [Line2D([0], [0], marker=REP_MARKERS[r], color="none", markerfacecolor=REP_COLORS[r],
                      markeredgecolor=REP_COLORS[r], markersize=4.8, label=f"Held R{r}") for r in (1, 2, 3)]
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.645, 0.982),
               ncol=3, handletextpad=0.35, columnspacing=1.25)
    tidy(ax)
    save(fig, "Supplementary_Figure_S1_MSH2_sensitivity")
    source = pd.concat([subset.assign(source_row="sensitivity"), primary.assign(source_row="primary_reference")],
                       ignore_index=True, sort=False)
    return subset, source


def figure_s2(data_root: Path) -> pd.DataFrame:
    bard_path = data_root / "revision_20261122" / "external_data" / "bard1_analysis_20261122" / "direction_summary_and_position_cluster_CI.csv"
    vhl_path = data_root / "revision_20261122" / "secondary_sources" / "vhl_rna_crossrep" / "qc" / "same_region_sensitivity.csv"
    bard, vhl = pd.read_csv(bard_path), pd.read_csv(vhl_path)
    if len(bard) != 9 or len(vhl) != 4:
        raise ValueError(f"expected 9 BARD1 and 4 VHL endpoints; found {len(bard)} and {len(vhl)}")
    fig = plt.figure(figsize=(FIGURE_WIDTH, 3.31))
    axes = [fig.add_axes([0.108, 0.222, 0.355, 0.62]),
            fig.add_axes([0.607, 0.222, 0.355, 0.62])]
    rows = []
    for panel_letter, ax, table, gene, n_groups in zip("AB", axes, [bard, vhl], ["BARD1", "VHL"], [580, 85]):
        panel(fig, panel_letter, 0.02 if panel_letter == "A" else 0.52, 0.966, f"{gene} ({n_groups} groups)")
        y = np.arange(len(table))[::-1]
        labels = []
        for i, row in table.reset_index(drop=True).iterrows():
            if gene == "BARD1":
                label = str(row.direction).replace(":", " / ")
                relation = str(row.replicate_relation)
                effect = num(row.mean_group_sample_covariance, "BARD1 covariance")
                low, high = num(row.cov_ci_lower, "BARD1 CI"), num(row.cov_ci_upper, "BARD1 CI")
                groups, positions, variants = 580, 445, 1240
            else:
                label = f"{row.rna_measure} / {row.function_measure}"
                relation = "within" if str(row.pairing).startswith("within") else "cross"
                effect = num(row.mean_group_covariance, "VHL covariance")
                low = num(row.mean_group_covariance_ci95_low, "VHL CI")
                high = num(row.mean_group_covariance_ci95_high, "VHL CI")
                groups, positions, variants = 85, 82, 179
            color = GREY if relation == "within" else BLUE
            labels.append(label)
            ax.errorbar(effect, y[i], xerr=np.array([[effect - low], [high - effect]]), fmt="o",
                        color=color, markersize=3.7, capsize=2.0, elinewidth=0.85, zorder=3)
            rows.append({"gene": gene, "direction": label, "replicate_relation": relation,
                         "mean_group_covariance_log2_units_squared": effect,
                         "ci_lower": low, "ci_upper": high, "groups": groups,
                         "protein_positions": positions, "DNA_target_observations": variants})
        ax.axvline(0, lw=0.7, linestyle=(0, (3, 2)), color=GREY, zorder=0)
        ax.set_yticks(y, labels)
        ax.set_xlim(-0.10, 0.115)
        ax.set_ylim(-0.6, len(table) - 0.3)
        ax.set_xticks([-0.10, -0.05, 0, 0.05, 0.10])
        ax.set_xlabel("Mean group covariance (log₂ units²)")
        tidy(ax)
    legend = [Line2D([0], [0], marker="o", color=GREY, linestyle="none", markersize=4, label="Same repeat"),
              Line2D([0], [0], marker="o", color=BLUE, linestyle="none", markersize=4, label="Different repeats")]
    fig.legend(handles=legend, loc="lower center", bbox_to_anchor=(0.5, 0.014),
               ncol=2, handletextpad=0.45, columnspacing=1.55)
    save(fig, "Supplementary_Figure_S2_RNA_cross_repeat_checks")
    return pd.DataFrame(rows)


def export_sources(data: dict[str, Any], data_root: Path, frames: dict[str, pd.DataFrame]) -> None:
    DEFAULT_SOURCE.mkdir(parents=True, exist_ok=True)
    counts, hero, states = frames["f1"]
    hist = frames["s3_hist"]
    counts.to_csv(DEFAULT_SOURCE / "figure1_collection_scope.csv", index=False)
    hero.to_csv(DEFAULT_SOURCE / "figure1_qualified_group_counts.csv", index=False)
    states.to_csv(DEFAULT_SOURCE / "figure1_full_group_states.csv", index=False)
    hist.to_csv(DEFAULT_SOURCE / "figure1_group_size_distribution.csv", index=False)
    frames["f2_scatter"].to_csv(DEFAULT_SOURCE / "figure2_brca2_scatter_points.csv", index=False)
    frames["f2_directions"].to_csv(DEFAULT_SOURCE / "figure2_direction_interval_points.csv", index=False)
    frames["f2_phi"].to_csv(DEFAULT_SOURCE / "figure2_phi_estimates.csv", index=False)
    frames["f3_metrics"].to_csv(DEFAULT_SOURCE / "figure3_msh2_all_endpoint_metrics.csv", index=False)
    frames["f3_attrition"].to_csv(DEFAULT_SOURCE / "figure3_msh2_attrition.csv", index=False)
    data["msh2_group_errors"].to_csv(DEFAULT_SOURCE / "figure3_msh2_group_errors.csv.gz", index=False, compression="gzip")
    frames["f4_null"].to_csv(DEFAULT_SOURCE / "figure4_null_simulation_tau0.csv", index=False)
    frames["f4_workflow"].to_csv(DEFAULT_SOURCE / "figure4_observed_workflows.csv", index=False)
    frames["s1_source"].to_csv(DEFAULT_SOURCE / "supplementary_s1_sensitivity_and_primary.csv", index=False)
    frames["s2_points"].to_csv(DEFAULT_SOURCE / "supplementary_s2_covariance_points.csv", index=False)

    rule = {
        "label": "Rule A two-variant analytic acceptance region",
        "predicate_exact": "abs(D1) > 1.96 * sqrt(2 * max((D1 - D2)^2 / 4, 1e-12))",
        "predicate_away_from_origin": "abs(D1 - D2) < (sqrt(2) / 1.96) * abs(D1)",
        "grid_floor_inside_max": 1e-12,
        "display": "analytic decision region; no empirical observations plotted",
    }
    (DEFAULT_SOURCE / "figure4_analytic_rule_definition.json").write_text(json.dumps(rule, indent=2) + "\n", encoding="utf-8")

    input_manifest = {}
    for name, path in {**base.INPUTS, **base.OPTIONAL_INPUTS}.items():
        if name == "selection_null_spec":
            rel = Path("analysis") / path.name
        else:
            rel = Path("data/derived") / path.relative_to(data_root)
        input_manifest[name] = {"relative_path": str(rel).replace("\\", "/"), "sha256": base.sha256(path)}
    for name, path in {
        "BARD1_covariance": data_root / "revision_20261122/external_data/bard1_analysis_20261122/direction_summary_and_position_cluster_CI.csv",
        "VHL_covariance": data_root / "revision_20261122/secondary_sources/vhl_rna_crossrep/qc/same_region_sensitivity.csv",
    }.items():
        rel = Path("data/derived") / path.relative_to(data_root)
        input_manifest[name] = {"relative_path": str(rel).replace("\\", "/"), "sha256": base.sha256(path)}
    manifest = {
        "backend": "Python / matplotlib",
        "input_root": "data/derived",
        "input_sha256": input_manifest,
        "outputs": [
            "figure1_collection_testability", "figure2_brca2_paired_contrasts",
            "figure3_msh2_controlled_endpoints", "figure4_structural_selection_counterexample",
            "Supplementary_Figure_S1_MSH2_sensitivity", "Supplementary_Figure_S2_RNA_cross_repeat_checks",
            "Supplementary_Figure_S3_collection_detail",
        ],
        "point_exports": sorted(p.name for p in DEFAULT_SOURCE.glob("*.csv*")),
    }
    (DEFAULT_SOURCE / "source_hash_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")

    preservation = {
        "status": "rendered from frozen verified inputs; no observations reanalyzed or altered",
        "main_figures": {
            "figure1": {"displayed_scope_counts": [92, 89, 27, 413090, 31723, 13041],
                        "hero_counts": hero.to_dict("records"),
                        "moved_to_supplementary_s3": ["historical group-size histogram", "full category/state counts for all three identifier tiers"],
                        "excluded_observations": 0},
            "figure2": {"brca2_scatter_groups": int(len(frames["f2_scatter"])),
                        "direction_interval_points": int(len(frames["f2_directions"])),
                        "phi_values": frames["f2_phi"].to_dict("records"), "excluded_observations": 0},
            "figure3": {"endpoint_gain_points": int(len(frames["f3_metrics"])),
                        "main_panel_absolute_rmse_scope": "6TG_mock AA vs DNA for each of 3 held replicates",
                        "all_endpoint_rmse_mae_values": "figure3_msh2_all_endpoint_metrics.csv; all 3 endpoints and 3 held replicates",
                        "attrition": "figure3_msh2_attrition.csv and original frozen attrition table",
                        "excluded_metrics": 0},
            "figure4": {"null_tau_over_sigma_zero_rules": int(len(frames["f4_null"])),
                        "observed_workflows": int(len(frames["f4_workflow"])),
                        "excluded_points": 0},
        },
        "supplementary": {
            "S1": {"sensitivity_rows": int(frames["s1_sensitivity"].groupby(["scenario", "pseudocount", "depth_frequency_threshold", "count_correction"]).ngroups),
                   "primary_reference_rows": int(frames["s1_source"]["source_row"].eq("primary_reference").sum()),
                   "primary_ci_count": int(frames["s1_source"].query("source_row == 'primary_reference'")["relative_rmse_gain_95pct_ci"].notna().sum())},
            "S2": {"BARD1_points": int(frames["s2_points"].gene.eq("BARD1").sum()),
                   "VHL_points": int(frames["s2_points"].gene.eq("VHL").sum())},
            "S3": {"group_size_bins": int(len(hist)), "state_rows": int(len(states)),
                   "all_three_tiers_sum_to_historical_group_denominator": True},
        },
        "source_files": input_manifest,
    }
    (DEFAULT_SOURCE / "figure_preservation.json").write_text(json.dumps(preservation, indent=2) + "\n", encoding="utf-8")


def main() -> None:
    global DEFAULT_OUT, DEFAULT_SOURCE
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, default=DEFAULT_DATA_ROOT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--source-dir", type=Path, default=DEFAULT_SOURCE)
    args = parser.parse_args()
    DEFAULT_OUT = args.output_dir.expanduser().resolve()
    DEFAULT_SOURCE = args.source_dir.expanduser().resolve()
    data_root = args.data_root.expanduser().resolve()
    base.configure_paths(data_root, DEFAULT_OUT, DEFAULT_SOURCE)
    base.INPUTS["selection_null_spec"] = TASK_ROOT / "analysis" / "selection_bias_null.py"
    base.check_inputs()
    data = base.load_inputs()

    scope, hero, states = figure1(data)
    hist = figure_s3(data, states)
    scatter, directions, phi = figure2(data)
    metrics, attrition = figure3(data)
    null0, workflow = figure4(data)
    sensitivity, s1_source = figure_s1(data)
    s2_points = figure_s2(data_root)
    export_sources(data, data_root, {
        "f1": (scope, hero, states), "s3_hist": hist, "f2_scatter": scatter,
        "f2_directions": directions, "f2_phi": phi, "f3_metrics": metrics,
        "f3_attrition": attrition, "f4_null": null0, "f4_workflow": workflow,
        "s1_sensitivity": sensitivity, "s1_source": s1_source, "s2_points": s2_points,
    })
    print(json.dumps({"figures_out": str(DEFAULT_OUT), "source_data": str(DEFAULT_SOURCE),
                      "figures": 7, "brca2_points": len(scatter), "msh2_metrics": len(metrics)}), flush=True)


if __name__ == "__main__":
    main()
