# Reproducing the manuscript results

Run commands from the repository root. Use Python 3.12. Install in a virtual
environment with `python -m pip install -r requirements.txt`. The exact Linux
verification environment is recorded separately; it uses distribution-packaged
libraries and is not presented as a complete PyPI lockfile.

For bounded CPU use on Linux/macOS:

```bash
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export MPLBACKEND=Agg PYTHONNOUSERSITE=1
```

In PowerShell set the equivalent `$env:OPENBLAS_NUM_THREADS='1'` variables.
Long computations and downloads for the original project were executed on its
analysis server. External users can run the same commands on their own system.

## 1. Package integrity and frozen scientific aggregates

```bash
python tests/verify_manifest.py
python tests/verify_frozen_results.py
python scripts/reproduce.py collection
```

These recalculate collection tiers, BRCA2 direction counts and Phi from group
moments, and MSH2 gains from group errors. They do not refetch the 92-set archive
or rerun MSH2 count reconstruction. `SHA256SUMS.txt` covers the release files,
excluding itself and ignored generated outputs.

## 2. Four main figures and three supplementary figures

```bash
python scripts/reproduce.py figures
```

The entrypoint runs `figures/build_display_figures.py` once for all four main
figures and Supplementary Figures S1–S3. Inputs are the frozen files under
`data/derived/`; PDFs, PNG previews, SVGs, source tables, the SHA-256 input
manifest and the operational preservation record are written under
`results/figures/` (source exports are in `results/figures/source_data/`).
Figure S1 uses the primary 6-TG/mock sensitivity rows; Figure S2 uses BARD1/VHL
summary tables; Figure S3 preserves the historical group-size distribution
and full identifier-state counts. The seven approved PDF files are kept under
`reference_figures/`. PDF creation dates, embedded fonts and binary hashes can
differ without a change in values. Compare source values and panel content.

The reference PDFs are 170 mm wide and use 18-point Arial Bold panel letters.
For the same typography when Arial is outside Matplotlib's font search path,
set `FIGURE_ARIAL_BOLD` to the path of the installed Arial Bold TrueType font.
The script uses DejaVu Sans Bold for glyphs unavailable in Arial or on hosts
without Arial. The reference PDFs embed the fonts used for each glyph.

The four main tables are supplied as `tables/table_1.csv` through
`table_4.csv` and `tables/main_tables.md`, preserving captions and notes.
These are presentation extracts of the frozen results, not independent
statistical analyses. Their underlying numeric sources are mapped in
`DATA_DICTIONARY.md`.

## 3. BRCA2 score-level analyses

```bash
python scripts/reproduce.py brca2-map
python scripts/reproduce.py brca2
python scripts/reproduce.py brca2-direction
python tests/verify_recomputed_core.py
```

The two source-specific mapped score tables and original CC0 exports are
bundled. The core command reconstructs the same-AA shared-key groups, sample
moments, primary Phi and ordinary group-bootstrap interval. The direction
command reconstructs primary/fixed/reverse selection, exact ties and
position-cluster direction intervals. It retains the historical Phi summary
as separately labeled provenance. The core recomputation is the independent
check of that Phi.

The `brca2-map` command reconstructs mappings from the two bundled raw score
exports using the pinned CDS and checks every mapping field against the archived
mapped tables. Results are written under `results/brca2_mapping/`. The main
analysis can be pointed to those regenerated tables with `--hap1` and `--mes`.
This is not a rerun of the source laboratories'
original count processing or the AVENGERS pipeline. Source-specific processed
score interpretation remains as stated in the manuscript.

## 4. MSH2 raw-count reconstruction

```bash
python scripts/fetch_inputs.py --dataset MSH2_GSE162130
python scripts/reproduce.py msh2
python tests/compare_recomputed.py --reference-dir data/derived/revision_20261123_review/msh2/output --recomputed-dir results/msh2_v4 --output-json results/verification/msh2_comparison.json
```

The acquisition helper enforces the frozen SHA-256 before accepting the GEO
file. It does not put downloaded source counts under this project's CC BY
license. The default analysis writes to `results/msh2_v4/`, reconstructs
6-TG/mock, mock/P0 and DB6/P0 on matched source-defined cohorts, and performs
5,000 primary position-bootstrap draws per held replicate with beta refit.
The benchmark cohort and outcomes must be compared with
`data/derived/revision_20261123_review/msh2/output/` using the comparator CLI.

This reconstruction uses the deposited count-table universe, tile-normalized
WT correction, source-only eligibility and five deterministic position folds.
It is not the original Jia et al. amino-acid LoF/sign-filter pipeline. Estimated
gains describe internal transfer among related culture replicates sharing a
starting pool and WT reference, not independent external validation.

## 5. Selection null

```bash
python scripts/reproduce.py selection
```

This reruns the prespecified null simulation with seed 20261115. The real
workflow source table remains separately frozen; workflow-specific denominators
and selection dependence must not be treated as matched predictive performance.

## 6. Optional archival audits

`collection-raw` and `splicing` expose current source-audit implementations;
run their `--help` to specify an external source project layout. Their complete
third-party raw inputs are not all bundled. See `ARCHIVAL_SCOPE.md` and the
92-record `collection_source_exports.csv`. Successful frozen-input checks
must not be described as a fresh 27-source raw-to-mapped audit.

BARD1 pinned inputs can be fetched with `--dataset BARD1_RNA` and
`--dataset BARD1_annotation`. The upstream repository lacks a license, so
those original tables and external code are not reissued here. VHL original
Table ST1 is obtained from its article. Summary-level Figure S2 reproduction
does not require either original count table.

The splicing historical master request returned HTTP403. The resulting error
page is not data and is not included. Event/transcript/strand/condition identity
remains unresolved; no biological validation is reconstructed from that audit.

## Failure handling

Stop on an input-hash mismatch, a failed frozen-count assertion or a missing
source field. Do not replace missing observations with simulated data, skip
the check, silently use a newer source version or relabel an old endpoint.
Regenerated numerical results belong in `results/`; never overwrite frozen
release inputs to make a comparison pass.

The file listing in `package_inventory.json` excludes itself and `SHA256SUMS.txt` to avoid self-referential hashes. `SHA256SUMS.txt` covers the package inventory and all other distributed files, excluding itself.
