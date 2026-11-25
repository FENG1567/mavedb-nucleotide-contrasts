# Publication figure layout

The seven PDFs in `reference_figures/` use the 2026-11-25 layout. Rebuild that
layout from the frozen plotted summaries with:

```bash
python figures/build_layout_figures.py
```

The code uses Python, matplotlib, numpy and pandas from the existing project
requirements. It reads `data/derived/figure_layout_inputs/` and writes PDF, SVG
and PNG to `results/figures/`. These tables retain all previously displayed
observations, categories, estimates and interval bounds; the figure builder
does not refit the scientific analyses. Use the existing analysis pipeline
to reconstruct those summaries from their upstream inputs.

Panel letters are uppercase Arial Bold at 18 pt in the standalone 7.2-inch
PDFs. Other labels are bold Arial, at least 9.5 pt. Word embeds scale the entire
figure proportionately. PDF and SVG retain selectable text.

A licensed Arial font must be installed. Alternatively, pass
`--font-dir /path/to/locally/licensed/Arial/fonts`. Proprietary font files are
not redistributed in this repository.

Original code is licensed under MIT. Author-generated summary tables and
figures retain CC BY 4.0. Third-party data keep the terms listed in
`THIRD_PARTY.md`.
