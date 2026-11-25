"""tau=0 selection-bias null experiment.

Compare three group-selection pipelines under a known null (tau=0):
A. validation-set-dependent: threshold from each group's noise_var (uses both experiments)
B. fixed external threshold: threshold from known sigma^2
C. cross-fitting: threshold estimated on a training half, applied to the test half

The quantity is R among selected groups (sign concordance), which should be
~50% under tau=0 if selection is independent of the test experiment.
"""
from __future__ import annotations

import json
import math
import argparse
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUT = REPO_ROOT / "results/selection_bias_null"
RNG = np.random.default_rng(20261115)

N2, N3 = 341, 47
N_GROUPS = N2 + N3
SIGMA2 = 0.346
TAU2_BASE = SIGMA2
N_REP = 1000


def sim(ratio: float):
    tau2 = (ratio ** 2) * SIGMA2
    sizes = np.array([2] * N2 + [3] * N3)
    var1 = np.empty(N_GROUPS); var2 = np.empty(N_GROUPS)
    cov = np.empty(N_GROUPS); d1 = np.empty(N_GROUPS); d2 = np.empty(N_GROUPS)
    for i, ni in enumerate(sizes):
        delta = RNG.normal(0, math.sqrt(tau2), size=ni)
        e1 = RNG.normal(0, math.sqrt(SIGMA2), size=ni)
        e2 = RNG.normal(0, math.sqrt(SIGMA2), size=ni)
        y1 = delta + e1; y2 = delta + e2
        x = y1 - y1.mean(); y = y2 - y2.mean()
        var1[i] = (x * x).sum() / (ni - 1)
        var2[i] = (y * y).sum() / (ni - 1)
        cov[i] = (x * y).sum() / (ni - 1)
        if ni == 2:
            d1[i], d2[i] = y1[0] - y1[1], y2[0] - y2[1]
        else:
            o = np.argsort(y1)
            d1[i], d2[i] = y1[o[-1]] - y1[o[0]], y2[o[-1]] - y2[o[0]]
    noise = np.maximum(0, 0.5 * (var1 + var2) - cov)
    conc = (d1 * d2 > 0).astype(float)

    def R(sel):
        return float(conc[sel].mean()) if sel.sum() > 0 else np.nan

    selA = np.abs(d1) > 1.96 * np.sqrt(2 * np.maximum(noise, 1e-12))
    thrB = 1.96 * math.sqrt(2 * SIGMA2)
    selB = np.abs(d1) > thrB
    # cross-fitting: 50/50 split, threshold from train noise only
    idx = RNG.permutation(N_GROUPS)
    tr, te = idx[: N_GROUPS // 2], idx[N_GROUPS // 2:]
    thrC = 1.96 * math.sqrt(2 * max(noise[tr].mean(), 1e-12))
    selC = np.zeros(N_GROUPS, dtype=bool)
    selC[te] = np.abs(d1[te]) > thrC
    return {
        "R_all": R(np.ones(N_GROUPS, dtype=bool)),
        "R_A": R(selA), "n_A": int(selA.sum()),
        "R_B": R(selB), "n_B": int(selB.sum()),
        "R_C": R(selC), "n_C": int(selC.sum()),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUT,
                        help="directory for computed simulation outputs")
    out = parser.parse_args().output_dir.expanduser().resolve()
    out.mkdir(parents=True, exist_ok=True)
    rows = []
    for ratio in [0.0, 1.0]:
        reps = [sim(ratio) for _ in range(N_REP)]
        for key in ["R_all", "R_A", "R_B", "R_C"]:
            vals = np.array([r[key] for r in reps], dtype=float)
            vals = vals[~np.isnan(vals)]
            rows.append({
                "tau_over_sigma": ratio,
                "pipeline": key,
                "mean_R": float(vals.mean()),
                "ci_low": float(np.quantile(vals, 0.025)),
                "ci_high": float(np.quantile(vals, 0.975)),
                "P_R_ge_60": float((vals >= 0.60).mean()),
                "mean_selected": float(np.mean([r["n_" + key[-1]] for r in reps])) if key != "R_all" else N_GROUPS,
            })
    df = pd.DataFrame(rows)
    df.to_csv(out / "selection_bias_null.csv", index=False)
    (out / "selection_bias_null.json").write_text(json.dumps(rows, indent=2))
    print(df.to_string(index=False))


if __name__ == "__main__":
    main()

