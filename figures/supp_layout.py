"""Layout-only redraws of the three BRCA2 supplementary figures.

Each function reads the already-exported source-data CSVs in ``inputs`` and
returns an open Matplotlib Figure. The observations, estimates, and confidence
intervals are not recomputed or filtered.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import pandas as pd
from matplotlib.lines import Line2D

import layout_common as lc


_REPLICATE_COLORS = {1: lc.BLUE, 2: lc.GREEN, 3: lc.ORANGE}
_REPLICATE_MARKERS = {1: "o", 2: "s", 3: "^"}
_REPLICATE_OFFSETS = {1: -0.15, 2: 0.0, 3: 0.15}
_S1_ORDER = [
    "Pseudocount 0.1",
    "Pseudocount 1",
    "Source frequency >1/400,000",
    "Source frequency >1/100,000",
    "Uncorrected counts",
]
_S1_ORDER_INDEX = {label: i for i, label in enumerate(_S1_ORDER)}
_TIER_ORDER = ["Accession reference", "Versioned request", "Target conditional"]
_STATE_COLORS = {
    "ge2": lc.BLUE,
    "one": lc.GREEN,
    "zero": lc.ORANGE,
    "unknown_0_or_1": lc.GREY,
    "conditional_ge2": lc.BLUE,
    "conditional_one": lc.GREEN,
    "zero_target_compatible_under_assumption": lc.ORANGE,
    "not_assessed_target_sequence": lc.GREY,
    "tpk1_full_construct_unresolved": lc.LIGHT,
}


def _read(inputs: Path, filename: str) -> pd.DataFrame:
    return pd.read_csv(Path(inputs) / filename)


def _require_columns(table: pd.DataFrame, columns: tuple[str, ...], filename: str) -> None:
    missing = [column for column in columns if column not in table.columns]
    if missing:
        raise ValueError(f"{filename} is missing required columns: {', '.join(missing)}")


def _num(value: object, label: str) -> float:
    number = float(value)
    if not math.isfinite(number):
        raise ValueError(f"non-finite {label}: {value!r}")
    return number


def _s1_setting_label(row: pd.Series) -> str:
    scenario = str(row["scenario"])
    if scenario.startswith("depth"):
        depth = _num(row["depth_frequency_threshold"], "depth-frequency threshold")
        if depth <= 0:
            raise ValueError(f"depth-frequency threshold must be positive: {depth}")
        denominator = round(1.0 / depth)
        return f"Source frequency >1/{denominator:,}"
    if scenario.startswith("pseudocount"):
        return f"Pseudocount {_num(row['pseudocount'], 'pseudocount'):g}"
    return "Uncorrected counts"


def _ci_pair(value: object, label: str) -> tuple[float, float]:
    try:
        pair = json.loads(str(value))
    except (TypeError, json.JSONDecodeError) as error:
        raise ValueError(f"could not read {label} as a two-value JSON array: {value!r}") from error
    if not isinstance(pair, list) or len(pair) != 2:
        raise ValueError(f"{label} must contain exactly two bounds: {value!r}")
    low, high = (_num(pair[0], f"{label} lower bound"), _num(pair[1], f"{label} upper bound"))
    if low > high:
        raise ValueError(f"{label} lower bound exceeds upper bound: {low} > {high}")
    return low, high


def draw_s1(inputs: Path):
    """Return the S1 sensitivity plot with all 15 sensitivity and 3 primary points."""
    filename = "supplementary_s1_sensitivity_and_primary.csv"
    table = _read(inputs, filename)
    _require_columns(
        table,
        (
            "endpoint",
            "source_row",
            "scenario",
            "pseudocount",
            "depth_frequency_threshold",
            "count_correction",
            "held_rep",
            "relative_rmse_gain",
            "relative_rmse_gain_95pct_ci",
        ),
        filename,
    )
    if set(table["endpoint"].astype(str)) != {"6TG_mock"}:
        raise ValueError("S1 source must contain only the original 6TG_mock endpoint")

    sensitivity = table.loc[table["source_row"].astype(str).eq("sensitivity")].copy()
    primary = table.loc[table["source_row"].astype(str).eq("primary_reference")].copy()
    if len(sensitivity) != 15 or len(primary) != 3:
        raise ValueError(f"expected 15 sensitivity and 3 primary rows; found {len(sensitivity)} and {len(primary)}")

    by_setting: dict[str, list[pd.Series]] = {label: [] for label in _S1_ORDER}
    for _, row in sensitivity.iterrows():
        label = _s1_setting_label(row)
        if label not in by_setting:
            raise ValueError(f"unexpected S1 setting: {label}")
        by_setting[label].append(row)
    if any(len(by_setting[label]) != 3 for label in _S1_ORDER):
        raise ValueError("each of the five S1 sensitivity settings must contain three held replicates")
    if set(int(_num(value, "held replicate")) for value in primary["held_rep"]) != {1, 2, 3}:
        raise ValueError("S1 primary reference must contain held replicates 1, 2, and 3 exactly once")

    fig = lc.newfig(width=7.2, height=3.65)
    ax = lc.axes(fig, [0.385, 0.17, 0.57, 0.66])
    lc.style(ax)
    labels = _S1_ORDER + ["Primary analysis"]
    labels[4] = "Without WT subtraction"
    ax.set_yticks(range(len(labels)), labels)
    ax.set_ylim(len(labels) - 0.42, -0.42)
    ax.set_xlabel("6-TG / mock relative RMSE gain (%)")
    ax.set_xlim(-1, 8)
    ax.set_xticks([-1, 0, 2, 4, 6, 8])
    ax.grid(axis="x", color=lc.LIGHT, linewidth=0.6)
    ax.axvline(0, color=lc.GREY, linewidth=0.9, linestyle=(0, (3, 2)), zorder=0)
    ax.axhline(len(_S1_ORDER) - 0.5, color=lc.LIGHT, linewidth=0.9, zorder=1)

    for y, label in enumerate(_S1_ORDER):
        rows = sorted(by_setting[label], key=lambda row: int(_num(row["held_rep"], "held replicate")))
        reps = [int(_num(row["held_rep"], "held replicate")) for row in rows]
        if sorted(reps) != [1, 2, 3]:
            raise ValueError(f"S1 setting {label} does not contain held replicates 1, 2, and 3 exactly once")
        for row in rows:
            rep = int(_num(row["held_rep"], "held replicate"))
            value = 100 * _num(row["relative_rmse_gain"], "S1 relative RMSE gain")
            ax.scatter(
                value,
                y + _REPLICATE_OFFSETS[rep],
                s=34,
                color=_REPLICATE_COLORS[rep],
                marker=_REPLICATE_MARKERS[rep],
                zorder=3,
            )

    for _, row in primary.sort_values("held_rep").iterrows():
        rep = int(_num(row["held_rep"], "held replicate"))
        value = 100 * _num(row["relative_rmse_gain"], "primary relative RMSE gain")
        low, high = _ci_pair(row["relative_rmse_gain_95pct_ci"], "primary relative RMSE gain CI")
        low, high = 100 * low, 100 * high
        if low > value or value > high:
            raise ValueError(f"primary held replicate {rep} estimate falls outside its original CI")
        ax.errorbar(
            value,
            len(_S1_ORDER) + _REPLICATE_OFFSETS[rep],
            xerr=[[value - low], [high - value]],
            fmt=_REPLICATE_MARKERS[rep],
            markersize=5.0,
            color=_REPLICATE_COLORS[rep],
            capsize=2.5,
            elinewidth=1.0,
            zorder=3,
        )

    fig.legend(
        handles=lc.legend_handles(),
        frameon=False,
        loc="upper center",
        bbox_to_anchor=(0.59, 0.96),
        ncol=3,
        handletextpad=0.45,
        columnspacing=1.25,
        prop={"family": "Arial", "weight": "bold", "size": 9.5},
    )
    return fig


def draw_s2(inputs: Path):
    """Return side-by-side BARD1 and VHL covariance plots with every source CI."""
    filename = "supplementary_s2_covariance_points.csv"
    table = _read(inputs, filename)
    _require_columns(
        table,
        (
            "gene",
            "direction",
            "replicate_relation",
            "mean_group_covariance_log2_units_squared",
            "ci_lower",
            "ci_upper",
            "groups",
        ),
        filename,
    )
    bard = table.loc[table["gene"].astype(str).eq("BARD1")].copy()
    vhl = table.loc[table["gene"].astype(str).eq("VHL")].copy()
    if len(bard) != 9 or len(vhl) != 4:
        raise ValueError(f"expected 9 BARD1 and 4 VHL rows; found {len(bard)} and {len(vhl)}")

    fig = lc.newfig(width=7.2, height=3.9)
    panels = [
        (lc.axes(fig, [0.17, 0.22, 0.34, 0.59], letter="A"), bard, "BARD1"),
        (lc.axes(fig, [0.63, 0.22, 0.34, 0.59], letter="B"), vhl, "VHL"),
    ]
    for ax, subset, gene in panels:
        lc.style(ax)
        n_groups = int(_num(subset["groups"].iloc[0], f"{gene} group count"))
        ax.set_title(f"{gene} ({n_groups:,} groups)", loc="center", pad=9)
        y_positions = list(range(len(subset) - 1, -1, -1))
        labels: list[str] = []
        for i, (_, row) in enumerate(subset.reset_index(drop=True).iterrows()):
            label = str(row["direction"])
            relation = str(row["replicate_relation"])
            if relation not in {"within", "cross"}:
                raise ValueError(f"unexpected {gene} replicate relation: {relation}")
            color = lc.GREY if relation == "within" else lc.BLUE
            estimate = _num(row["mean_group_covariance_log2_units_squared"], f"{gene} covariance")
            low = _num(row["ci_lower"], f"{gene} CI lower bound")
            high = _num(row["ci_upper"], f"{gene} CI upper bound")
            if low > estimate or estimate > high:
                raise ValueError(f"{gene} estimate falls outside its original CI for {label}")
            ax.errorbar(
                estimate,
                y_positions[i],
                xerr=[[estimate - low], [high - estimate]],
                fmt="o",
                markersize=5.0,
                color=color,
                capsize=2.5,
                elinewidth=1.0,
                zorder=3,
            )
            labels.append(label)

        ax.axvline(0, linewidth=0.9, linestyle=(0, (3, 2)), color=lc.GREY, zorder=0)
        ax.set_yticks(y_positions, labels)
        ax.set_ylim(-0.6, len(subset) - 0.3)
        ax.set_xlim(-0.10, 0.115)
        ax.set_xticks([-0.10, -0.05, 0, 0.05, 0.10])
        ax.grid(axis="x", color=lc.LIGHT, linewidth=0.6)

    relation_legend = [
        Line2D([], [], marker="o", color=lc.GREY, linestyle="none", markersize=5, label="Same repeat"),
        Line2D([], [], marker="o", color=lc.BLUE, linestyle="none", markersize=5, label="Different repeats"),
    ]
    fig.legend(
        handles=relation_legend,
        frameon=False,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.98),
        ncol=2,
        handletextpad=0.45,
        columnspacing=1.5,
        prop={"family": "Arial", "weight": "bold", "size": 9.5},
    )
    fig.text(
        0.5,
        0.065,
        r"Mean group covariance ($\log_{2}$ units$^{2}$)",
        ha="center",
        va="center",
        fontsize=10,
        fontname="Arial",
        fontweight="bold",
        color=lc.INK,
    )
    return fig


def draw_s3(inputs: Path):
    """Return the stacked S3 histogram and all 13 collection-state categories."""
    hist_filename = "figure1_group_size_distribution.csv"
    states_filename = "figure1_full_group_states.csv"
    hist = _read(inputs, hist_filename)
    states = _read(inputs, states_filename)
    _require_columns(hist, ("keys_per_group", "groups"), hist_filename)
    _require_columns(states, ("tier", "state", "state_label", "groups"), states_filename)
    if len(hist) != 5:
        raise ValueError(f"expected 5 historical group-size bins; found {len(hist)}")

    ordered = [states.loc[states["tier"].astype(str).eq(tier)] for tier in _TIER_ORDER]
    if [len(part) for part in ordered] != [4, 4, 5]:
        raise ValueError(f"expected 4, 4, and 5 state rows for the three tiers; found {[len(part) for part in ordered]}")
    plot_df = pd.concat(ordered, ignore_index=True)
    if len(plot_df) != 13:
        raise ValueError(f"expected 13 collection-state rows; found {len(plot_df)}")
    for tier, group in plot_df.groupby("tier", sort=False):
        total = sum(int(_num(value, f"{tier} group count")) for value in group["groups"])
        if total != 13_041:
            raise ValueError(f"{tier} state counts sum to {total}, expected 13,041")

    fig = lc.newfig(width=7.2, height=6.6)
    ax0 = lc.axes(fig, [0.16, 0.76, 0.80, 0.17], letter="A")
    ax1 = lc.axes(fig, [0.49, 0.105, 0.46, 0.515])
    label = fig.text(0.025, 0.638, 'B', fontsize=18, fontname='Arial', fontweight='bold', va='bottom')
    label.set_gid('panel-label-B')
    lc.style(ax0)
    lc.style(ax1)

    hist_keys = [int(_num(value, "historical DNA keys per group")) for value in hist["keys_per_group"]]
    hist_counts = [int(_num(value, "historical group count")) for value in hist["groups"]]
    hist_x = list(range(len(hist)))
    ax0.bar(hist_x, hist_counts, color=lc.BLUE, width=0.68, edgecolor="white", linewidth=0.6)
    ax0.set_xticks(hist_x, [str(key) for key in hist_keys])
    ax0.set_xlabel("Historical DNA keys per amino-acid group")
    ax0.set_ylabel("Groups")
    ax0.set_ylim(0, max(hist_counts) * 1.27)
    ax0.set_yticks([0, 5_000, 10_000])
    ax0.grid(axis="y", color=lc.LIGHT, linewidth=0.6)
    for x, count in zip(hist_x, hist_counts):
        ax0.text(
            x,
            count + max(hist_counts) * 0.025,
            f"{count:,}",
            ha="center",
            va="bottom",
            fontsize=9.5,
            fontname="Arial",
            fontweight="bold",
            color=lc.INK,
            bbox={"facecolor": "white", "edgecolor": "none", "pad": 1.2},
        )

    state_counts = [int(_num(value, "collection-state group count")) for value in plot_df["groups"]]
    state_y = list(range(len(plot_df) - 1, -1, -1))
    bar_colors = [_STATE_COLORS.get(str(state), lc.GREY) for state in plot_df["state"]]
    ax1.barh(state_y, state_counts, color=bar_colors, height=0.72, edgecolor="white", linewidth=0.6)
    ax1.set_yticks(state_y, [str(label) for label in plot_df["state_label"]])
    ax1.set_ylim(-0.65, len(plot_df) - 0.35)
    ax1.set_xlim(0, 16_000)
    ax1.set_xticks([0, 5_000, 10_000, 15_000])
    ax1.set_xlabel("Groups (denominator 13,041 for each tier)")
    ax1.grid(axis="x", color=lc.LIGHT, linewidth=0.6)
    ax1.tick_params(axis="y", pad=5)
    for i, count in enumerate(state_counts):
        x = max(count + 120, 140)
        ax1.text(
            x,
            state_y[i],
            f"{count:,}",
            ha="left",
            va="center",
            fontsize=9.5,
            fontname="Arial",
            fontweight="bold",
            color=lc.INK,
            bbox={"facecolor": "white", "edgecolor": "none", "pad": 1.2},
        )

    offset = 0
    for tier, part in zip(_TIER_ORDER, ordered):
        row_count = len(part)
        tier_positions = state_y[offset : offset + row_count]
        center = (max(tier_positions) + min(tier_positions)) / 2
        ax1.text(
            -0.55,
            center,
            tier,
            transform=ax1.get_yaxis_transform(),
            ha="right",
            va="center",
            fontsize=9.5,
            fontname="Arial",
            fontweight="bold",
            color=lc.INK,
            clip_on=False,
        )
        offset += row_count

    ax1.axhline(8.5, color=lc.LIGHT, linewidth=0.8, zorder=0)
    ax1.axhline(4.5, color=lc.LIGHT, linewidth=0.8, zorder=0)
    fig.text(0.49, 0.015, '* Assumes c.1 is target base 1', fontsize=9.5, fontweight='bold')
    return fig
