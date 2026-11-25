# BRCA2 score-source and direction provenance

## Independent functional datasets

| MaveDB score set | Study unit and assay | Score payload and source processing | Publication |
|---|---|---|---|
| urn:mavedb:00001225-a-1 | Huang et al.; endogenous BRCA2 saturation genome editing in human HAP1 cells, with the MaveDB experiment method specifying LIG4-KO HAP1 cells, Cas9/gRNA and an all-SNV plasmid library. Cells were sampled at day 4 and day 14 for direct sequencing. | The exported CSV has only accession, hgvs_nt, score. MaveDB method metadata describes the continuous functional score as a log2 day-14/day-0 variant-read-frequency ratio, with Loess positional correction from D0–D5 comparisons, filtering below 10 reads, within-exon scaling against silent/nonsense benchmarks, replicate averaging, and cross-exon linear normalization. This is a processed assay score, not a raw read count or an unadjusted read-frequency table. | Huang H et al. Nature 2025;638:528–537. PMID 39779857. DOI 10.1038/s41586-024-08388-8. |
| urn:mavedb:00001242-a-1 | Sahu et al.; CRISPR saturation genome editing in a mouse ES cell line expressing one human BRCA2 copy, using sgRNA-Cas9 and an ssODN library. Cells were sampled at day 3 and day 14 for direct sequencing. | The exported CSV has only accession, hgvs_nt, score. MaveDB method metadata identifies the cell-survival assay and time points but does not specify the full score-normalization formula. The score values have not been individually reconciled to a supplementary-table/model-adjusted layer. MaveDB lists three separate calibration records; these are not additional columns in the exported CSV and were not applied in this direction analysis. | Sahu S et al. Nature 2025;638:538–545. PMID 39779848. DOI 10.1038/s41586-024-08349-1. |

These are two distinct study/assay units, not technical replicates of one cell system: HAP1 LIG4-KO versus a humanized mouse ES line. Each exported row is a variant score; replicate counts are not represented as separate score rows in the three-column files.

## Data sources and analysis

The source data are MaveDB score sets urn:mavedb:00001225-a-1 and urn:mavedb:00001242-a-1. Their source records provide assay and score-method descriptions. The corresponding direction summaries and group-level results are supplied as brca2_direction_summary.json and brca2_direction_group_sensitivities.csv.

Mapped records retain the source score and add transcript, codon, amino-acid, protein-position and class mapping. Analysis converts score and protein position to numeric, excludes unmapped or nonfinite rows, retains missense rows, then computes a separate z score for each assay across that assay's full mapped missense distribution: z = (score − mean) / population SD (ddof=0). Variants are grouped by protein position, reference amino acid and alternate amino acid, and intersected by exact hgvs_nt.

The shared universe is 388 amino-acid groups, 823 exact HGVS keys and 367 protein positions (341 two-variant groups, 47 groups with three variants). The principal statistic R is the share of groups whose selected HAP1 and mES numeric score differences have the same sign. It is a within-group rank-direction statistic, not a clinical classification metric. The Sahu mES biological score direction is summarized below. The HAP1 export was analyzed in its numeric source orientation without recoding it to clinical labels.

## Direction and tie definitions

- Primary estimator: for n=2, use the first two lexically sorted HGVS keys and calculate key 1 minus key 2 in each assay. For n≥3, select the lowest and highest HAP1 z-scored variants and calculate high minus low in HAP1 and the same two variants in mES. The product of the two differences determines concordance.
- Reverse-source sensitivity: for n≥3, select mES-score extremes and evaluate HAP1 on those same variants; the two-variant ordering remains fixed.
- Score-independent fixed-pair sensitivity: use the first two exact HGVS keys in every group, with the same key order in both assays.
- Ties are exact equality of exported score values within a group; per-assay z-scoring preserves equality. There were no within-group score ties, no selected-axis zero differences and no score-extreme tie-break ambiguities in the 388-group universe.
- With no ties, the tie-conditioned chance expectation is 0.5. The JSON reports the conditional exchangeability benchmark as an expectation, not a p-value.

## Direction results

- Primary BRCA2 Phi = 0.8750879 (87.51%), with group-bootstrap 95% interval [0.7680, 0.9866].
- HAP1-extreme direction: R = 193/388 = 49.74%; position-cluster bootstrap 95% interval [0.4479, 0.5476].
- Reverse-source sensitivity: R = 193/388 = 49.74%; position-cluster interval [0.4487, 0.5464].
- Fixed-first-two sensitivity: R = 196/388 = 50.52%; position-cluster interval [0.4552, 0.5544].
- In the n=2 stratum, primary and fixed-pair R = 171/341 = 50.15%. In n≥3, primary R = 22/47 = 46.81%; fixed-pair R = 25/47 = 53.19%.

The source-dependent primary/reverse estimators are near the 0.5 chance benchmark. The fixed-pair result is also near 0.5. These results do not support a general above-chance direction-transfer claim for the two score sets.

## Score provenance limits

The MaveDB CSV exposes a single continuous score field per score set; it does not contain sequencing counts, replicate-level values, VarCall categories, or the three calibration outputs as separate columns. For Huang, the MaveDB method text documents normalization steps, but exported row values have not been individually matched to the paper's supplemental score layer. For Sahu, the article and AVENGERS author code distinguish the continuous Function Score (FS), averaged across two replicates and normalized against synonymous/nonsense controls, from downstream GMM Probability/final.proba and categorical calls. The biological direction is higher score = more retained cell fitness/BRCA2 function; lower score = dropout/loss of function. The generic MaveDB field score is not yet mapped exactly to the code column Function.score versus final.FS, and the source record lists 6,545 variants versus 6,551 reported function scores. The MaveDB values should not be relabeled as PIF, a posterior, a VarCall output or clinical evidence. The two assay scores are source-specific processed functional readouts, not harmonized clinical labels or independent replicate measurements.

## Sahu score layer and biological direction

Sahu et al. describe day 3 as the initial variant representation and day 14 as the later abundance measurement. The Function Score (FS) is the relative dropout between these time points, averaged across two independent replicates. The author AVENGERS pipeline normalizes against synonymous and nonsense controls and keeps continuous FS (Function.score/final.FS) separate from GMM probability (Probability/final.proba) and classification. Higher score corresponds to greater retained cell fitness/BRCA2 function; lower score corresponds to dropout/loss of function. The current MaveDB score field is generic and cannot yet be equated uniquely with Function.score versus final.FS; its metadata reports 6,545 variants while the article reports 6,551 function scores.
