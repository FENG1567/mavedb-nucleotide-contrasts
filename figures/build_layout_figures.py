"""Rebuild the seven publication figures from frozen plotted summaries (MIT)."""
import argparse
from pathlib import Path
from matplotlib import font_manager
import layout_common as lc
import main_layout as main
import supp_layout as supp

def run():
    repo=Path(__file__).resolve().parents[1]
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-dir',type=Path,default=repo/'data/derived/figure_layout_inputs')
    parser.add_argument('--output-dir',type=Path,default=repo/'results/figures')
    parser.add_argument('--font-dir',type=Path,help='Optional directory containing locally licensed Arial TTF files')
    args=parser.parse_args()
    if args.font_dir:
        for f in args.font_dir.glob('*.ttf'):font_manager.fontManager.addfont(str(f))
    try:font_manager.findfont('Arial',fallback_to_default=False)
    except ValueError as exc:raise SystemExit('Arial is required for this figure style. Install a licensed Arial font or use --font-dir.') from exc
    lc.OUTPUT_DIR=args.output_dir
    specs=[(main.draw_f1,'figure1_collection_testability'),(main.draw_f2,'figure2_brca2_paired_contrasts'),
           (main.draw_f3,'figure3_msh2_controlled_endpoints'),(main.draw_f4,'figure4_structural_selection_counterexample'),
           (supp.draw_s1,'Supplementary_Figure_S1_MSH2_sensitivity'),(supp.draw_s2,'Supplementary_Figure_S2_RNA_cross_repeat_checks'),
           (supp.draw_s3,'Supplementary_Figure_S3_collection_detail')]
    for draw,stem in specs:
        lc.finish(draw(args.input_dir),stem)
        print('Rendered',stem,flush=True)

if __name__=='__main__':run()
