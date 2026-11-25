# Original-count secondary analyses for Figure S2

The scripts under `analysis/secondary/` rerun the BARD1 and VHL analyses from
source counts. They complement the Figure S2 summary renderer, which reads the
frozen summary tables under `data/derived/`. These are source-level, within-study
replicate checks; neither analysis establishes causality or independent-study
validation.

The analysis and QC implementations are project-authored code. The portable
release copies preserve the prespecified calculation and bootstrap logic while
accepting local input and output paths. They do not import or copy code from the
upstream BARD1 repository.

## BARD1

### Inputs and terms

The defaults are the following files under `data/external/bard1/`:

| File | Use | SHA-256 |
|---|---|---|
| `20250922_BARD1RNAscores.tsv` | Sole RNA, D05 and D13 count source | `f2ce427157fe00103835bc38e26f6693e285e68f807d6b6aadc8c3630d3dc2a5` |
| `20260807_BARD1.allscores.tsv` | Annotation fields only; score values are unused | `3f07fdede5e5d4b05329bfe328c9ccca89adc7f3363b0340cd460d4b475d7441` |
| `20260807_BARD1.snvcounts.tsv` | D05 version comparison only; counts do not enter ratios | `057657fd2c61be9e66db79e78b782005a945687fdd63282cc86c4aab98c31571` |

All three files come from `ivanw314/BARD1_SGE_analysis` at commit
`72ba50c8480a8dd88528c811ccc906d41fb3d971`. The upstream repository has no
identified license. All three are listed as fetch-only in
`data/external/sources.csv`. Do not bundle or redistribute these source tables.
The script writes an
`input_manifest.csv` containing each supplied path, byte size, SHA-256, pinned
commit, and use.

Fetch the manifest-listed inputs and run the analysis with the files present at
their default paths:

```bash
python scripts/fetch_inputs.py --dataset BARD1_RNA --dataset BARD1_annotation --dataset BARD1_snvcounts
python analysis/secondary/analyze_bard1_replicate_ratio.py
```

The script accepts `--legacy-counts`, `--annotation`, `--current-counts`, and
`--output-dir` to supply alternate locations. Its defaults are the three paths
above and `results/secondary/bard1/` for generated files.

### Frozen calculation

The legacy count table supplies all RNA, D05 and D13 values. Annotation is
joined by exact genomic SNV identity and is used to identify missense changes.
Eligibility requires all three legacy D05 replicate counts to be present and at
least 10, followed by at least two distinct DNA variants within each measured
`(target, amino-acid change)` group. There is no RNA-depth, D13-depth, RNA-score
or functional-score filter. Zero RNA and D13 counts are retained with the fixed
`+1` pseudocount; missing values remain missing.

For replicate `i`, `R_i = log2(rna_count_Ri + 1) - log2(D05_Ri + 1)` and
`F_i = log2(D13_Ri + 1) - log2(D05_Ri + 1)`. All three within-replicate and all
six ordered cross-replicate directions are reported. Group-specific sample
covariances and variances are cached; the primary summaries weight eligible
target/amino-acid groups equally. Percentile intervals use 2,000 protein-position
cluster-bootstrap draws with seed `20261122`. The newer count table is read only
for the exact-key D05 version comparison.

The output directory includes the direction summary and confidence intervals,
eligible variants and ratios, group moments, first-two-variant signs, bootstrap
draws, source/version checks, route status, report and input manifest.

## VHL

### Input and terms

Supply the supplementary ZIP for article `PMC11250436`, containing member
`41588_2024_1800_MOESM4_ESM.xlsx` and worksheet `ST1`. The ZIP is retrieved from
`https://www.ebi.ac.uk/europepmc/webservices/rest/PMC11250436/supplementaryFiles`.
Its SHA-256 is
`f9e7f13c5722c7dbaea77d6ec2e8e7ff7900d956c4e67d887882d17c04641087`; the ST1
member SHA-256 is
`6a432406170509777fc70d84dcf2288cf8c54cb434da18f570e3a07c4e087a85`. The source article and
supplement are identified as CC BY 4.0 in `data/external/sources.csv`; retain
the article attribution and check any individual-credit notices that apply to
the supplied material. This ZIP is not downloaded by the analysis script.
`openpyxl` is required to read the workbook and is declared in
`requirements.txt`.

Run the main analysis, then the independent count/region QC using the same ZIP
and output directory:

```bash
python analysis/secondary/analyze_vhl_crossrep.py \
  --input-zip /path/to/PMC11250436_supplementaryFiles.zip
python analysis/secondary/qc_vhl_regions.py \
  --input-zip /path/to/PMC11250436_supplementaryFiles.zip
```

Both scripts accept `--output-dir`; the default is `results/secondary/vhl/`.
The main analysis requires explicit `--input-zip`. QC also accepts `--main-script`
to identify the main-analysis script used. By default it hashes the neighboring
`analyze_vhl_crossrep.py`. The resulting `qc/qc_summary.json` records that
script's actual path and SHA-256. The release-portable script has changed path
plumbing relative to the previously executed server script, so this record is
an identity for the current run; it does not claim preservation of the old
script hash.

### Frozen calculation and QC

The main analysis retains missense records with exact genomic and HGVS
identifiers, reconciles duplicate exact-DNA rows only when annotations and all
six replicate counts agree, and applies day-6 DNA count of at least 10 in both
replicates. A retained protein-HGVS group must have at least two distinct DNA
variants. RNA and day-20 zeros remain and ratios use `log2(count + 1)`. The four
within- and cross-replicate pairings are all reported, with mean group sample
covariance, variant-weighted group-centered correlation, and first-two-variant
sign concordance. The percentile intervals use 2,000 protein-position cluster
bootstrap draws with seed `20261122`.

The QC script checks integer, finite, nonnegative counts for the retained rows;
verifies the unique mapping from each protein-HGVS group to protein position;
audits SGE target-region membership; and, when groups span regions, computes a
same-region sensitivity analysis with 2,000 protein-position bootstrap draws
and seed `20261123`. It records hashes for the source workbook member, main
analysis outputs, and the actual main-analysis script used in the run.

The main output directory contains `variant_level.csv`, `group_audit.csv`,
`crossrep_summary.csv`, and `summary.json`. QC writes `region_by_group.csv`,
`same_region_sensitivity.csv` when the sensitivity is applicable, and
`qc_summary.json` under `qc/`.
