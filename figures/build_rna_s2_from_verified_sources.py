"""Show all prespecified RNA/function covariance directions from both studies."""
import argparse
import hashlib
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def main():
    p = argparse.ArgumentParser()
    repo_root = Path(__file__).resolve().parents[1]
    p.add_argument('--data-root', type=Path, default=repo_root / 'data/derived/revision_20261122',
                   help='frozen 20261122 derived data root')
    p.add_argument('--output-dir', type=Path, default=repo_root / 'results/figures',
                   help='directory for generated figure assets')
    args = p.parse_args()
    root = args.data_root.expanduser().resolve()
    out = args.output_dir.expanduser().resolve()
    sources = [root / 'external_data/bard1_analysis_20261122/direction_summary_and_position_cluster_CI.csv',
               root / 'secondary_sources/vhl_rna_crossrep/qc/same_region_sensitivity.csv']
    missing = [str(path) for path in sources if not path.is_file()]
    if missing:
        raise SystemExit('Missing required frozen RNA/function input(s):\n  - ' + '\n  - '.join(missing))
    bard, vhl = [pd.read_csv(path) for path in sources]
    assert len(bard) == 9 and len(vhl) == 4
    out.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({'font.family': 'sans-serif',
                         'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
                         'font.size': 7, 'axes.titlesize': 7.5, 'xtick.labelsize': 6.5,
                         'ytick.labelsize': 6.5, 'svg.fonttype': 'none', 'pdf.fonttype': 42,
                         'axes.spines.top': False, 'axes.spines.right': False,
                         'axes.linewidth': .7, 'legend.frameon': False})
    fig, axes = plt.subplots(1, 2, figsize=(7.2047244, 3.3))
    fig.subplots_adjust(left=.18, right=.98, top=.81, bottom=.20, wspace=.68)
    exported = []
    for panel, ax, table, gene in zip('ab', axes, [bard, vhl], ['BARD1', 'VHL']):
        y = np.arange(len(table))[::-1]
        labels = []
        for i, row in table.iterrows():
            if gene == 'BARD1':
                label = row['direction'].replace(':', ' / ')
                relation = row['replicate_relation']
                effect, low, high = row['mean_group_sample_covariance'], row['cov_ci_lower'], row['cov_ci_upper']
                groups, positions, variants = 580, 445, 1240
            else:
                label = row['rna_measure'] + ' / ' + row['function_measure']
                relation = 'within' if row['pairing'].startswith('within') else 'cross'
                effect, low, high = row['mean_group_covariance'], row['mean_group_covariance_ci95_low'], row['mean_group_covariance_ci95_high']
                groups, positions, variants = 85, 82, 179
            color = '#777777' if relation == 'within' else '#31688e'
            labels.append(label)
            ax.errorbar(effect, y[i], xerr=np.array([[effect-low], [high-effect]]),
                        fmt='o', color=color, markersize=3.6, capsize=2, elinewidth=.9)
            exported.append({'gene': gene, 'direction': label, 'relation': relation,
                             'covariance': effect, 'ci_lower': low, 'ci_upper': high,
                             'target_AA_groups': groups, 'protein_positions': positions,
                             'DNA_target_observations': variants})
        ax.set_yticks(y)
        ax.set_yticklabels(labels)
        ax.axvline(0, lw=.7, linestyle='--', color='.65', zorder=0)
        ax.axhline(5.5 if gene == 'BARD1' else 1.5, color='.85', lw=.6)
        ax.set_xlim(-.10, .115)
        ax.set_ylim(-.6, len(table)-.3)
        ax.set_xticks([-.10, -.05, 0, .05, .10])
        ax.set_xlabel('Mean group covariance (log2 units²)')
        ax.set_title(f'{gene}: {groups} region × AA groups', loc='left', pad=13)
        ax.text(-.28, 1.13, panel, transform=ax.transAxes, fontsize=9, fontweight='bold')
        ax.tick_params(length=2.5, width=.6)
    fig.text(.18, .91, 'Grey: same repeat     Blue: different repeats', fontsize=7)
    fig.savefig(out / 'rna_cross_repeat_checks.svg', facecolor='white')
    fig.savefig(out / 'rna_cross_repeat_checks.pdf', facecolor='white')
    fig.savefig(out / 'rna_cross_repeat_checks.png', dpi=600, facecolor='white')
    fig.savefig(out / 'rna_cross_repeat_checks.tiff', dpi=600, facecolor='white',
                pil_kwargs={'compression': 'tiff_lzw'})
    plt.close(fig)
    pd.DataFrame(exported).to_csv(out / 'rna_cross_repeat_source_data.csv', index=False)
    (out / 'rna_cross_repeat_source_hashes.json').write_text(json.dumps(
        {str(path.relative_to(root)): hashlib.sha256(path.read_bytes()).hexdigest() for path in sources}, indent=2))
    (out / 'rna_cross_repeat_legend.md').write_text('''# Supplementary review figure | Same-repeat and cross-repeat RNA/function associations

All three within-repeat and six cross-repeat BARD1 directions (a), and all four
VHL directions after matching regions within amino-acid groups (b), are shown.
R denotes RNA/DNA at the baseline timepoint; F denotes later/baseline DNA. Each
score is a log2 count ratio with a +1 pseudocount. Points are mean sample
covariances within region-by-amino-acid groups, weighted equally. Bars are
pointwise 95% percentile intervals from 2,000 protein-position cluster bootstrap
draws. BARD1 includes 580 groups, 1,240 target-specific observations (997 unique
DNA variants) and 445 positions; VHL includes 85 groups, 179 DNA variants and
82 positions. BARD1 target overlaps remain together within a position cluster.
Same-repeat ratios share the observed baseline DNA count; cross-repeat ratios
use different baseline count vectors. Shared biological/technical effects can
remain across repeats. Neither the contrast nor its difference identifies a
causal mechanism or a proportion of covariance attributable to count noise.
Group counts and intervals condition on these studies and their supplied tables,
not independent populations of biological repeats. Different genes are shown as
separate study-specific checks, not as a test of differences between genes.
No direction was selected by its association. Full correlations, sign summaries,
count/annotation attrition and the exploratory paired covariance contrast are
provided in the accompanying result tables. These are additional review assets.
''')
    print(json.dumps({'outputs': str(out), 'source_rows': len(exported)}), flush=True)


if __name__ == '__main__':
    main()
