# Nucleotide variants encoding the same amino acid substitution

Reproducibility materials for *Different DNA Variants, Same Amino Acid Substitution: A MaveDB Analysis Including the Pancreatic Cancer Susceptibility Gene BRCA2*.

Public repository: https://github.com/FENG1567/mavedb-nucleotide-contrasts. Version: [v1.0.0](https://github.com/FENG1567/mavedb-nucleotide-contrasts/releases/tag/v1.0.0). See [中文说明](README_zh.md).

This repository separates original analysis code, author-generated results and
reused source data. **Original code: MIT. Author-generated summaries: CC BY
4.0.** Third-party inputs keep their original terms; see
[THIRD_PARTY.md](THIRD_PARTY.md). Cite the manuscript and original datasets using
[CITATION.cff](CITATION.cff) and the source manifest.

## Contents

| Directory | Purpose |
|---|---|
| `analysis/` | Current analysis implementations and documented optional audits |
| `figures/` | Figure-building code |
| `scripts/` | Execution/acquisition helpers |
| `tests/` | Frozen-result and package-integrity checks |
| `data/derived/` | Frozen author-generated scientific summaries and diagnostics |
| `data/external/` | Two CC0 BRCA2 source exports, pinned inputs and acquisition manifest |
| `tables/` | Four main tables in Markdown and CSV |
| `reference_figures/` | The seven approved figure PDFs |
| `docs/` | Input/output mapping, license boundaries and verification scope |

Generated figures and their source exports go to `results/figures/`, which is ignored by Git. Manuscript DOCX,
internal review correspondence, local environment credentials, raw large
sequencing files and operational logs are not part of this repository.

## Start here

Use Python 3.12 and a virtual environment. See the observed package versions
in `docs/verification_environment.json` and the install instructions in
[docs/REPRODUCING.md](docs/REPRODUCING.md).

```bash
python -m pip install -r requirements.txt
python tests/verify_manifest.py
python tests/verify_frozen_results.py
python scripts/reproduce.py --help
```

The first scientific check recalculates manuscript-facing aggregates from
frozen tables. It does not constitute raw-source reconstruction. The guide
distinguishes figure rebuilding, BRCA2 score-level recomputation, MSH2 raw-count
reconstruction and optional archival analyses. Consult
[docs/VALIDATION.md](docs/VALIDATION.md) for the checks actually completed.

## Result identity

The main BRCA2 Phi is **87.51%**, not the class-control sensitivity of 87.59%.
The MSH2 primary endpoint is **6-TG/mock**, with held-replicate RMSE gains of
4.18%, 4.48% and 4.20%. These are within-study related-replicate transfers.
Full definitions and scope are recorded in
[docs/SCIENTIFIC_CONTRACT.md](docs/SCIENTIFIC_CONTRACT.md).

Input hashes are checked before a frozen source can be substituted. Historical
online records can change; a changed hash is not accepted as equivalent data.
PDF binary hashes may differ across rendering environments; compare source
values, panel content and text before treating a rendering difference as a
scientific discrepancy.

## Current figure layout

The seven reference PDFs use the 2026-11-25 layout. Rebuild them with
`python figures/build_layout_figures.py`. See [figure layout instructions](figures/README_LAYOUT.md)
for the frozen summary inputs, font requirement and export formats.
