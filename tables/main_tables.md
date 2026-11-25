# Main tables

Four compact tables accompany the four main figures. Full collection details, all MSH2 endpoints, and selection-rule details are reported in Supplementary Tables S4–S6.

## Table 1. Collection-level testability tiers

| Quantity | Count | Unit / scope |
| --- | --- | --- |
| Manifest-defined MaveDB collection | 92 | score sets |
| Analyzable score sets | 89 | score sets |
| Numeric missense rows | 413,090 | numeric rows |
| Historical group-key entries | 31,723 | within-score-set group-key entries |
| Same-amino-acid groups with ≥2 keys | 13,041 | within-score-set groups |
| Score sets with ≥1 multi-key group | 27 | score sets |
| Reference-compatible AND source-scored groups | 3,836 | 19 sets; same groups as reference tier |
| Versioned-request-supported groups | 2,679 | 13 sets |
| Conditional target-sequence-compatible groups | 1,013 | 7 sets; excludes TPK1 |

Counts use different units and scopes. The 3,836 groups are the same set under reference compatibility and source scoring. The versioned-request tier is narrower; the conditional target-sequence tier is reported separately. Source scores do not establish independent physical measurements.

## Table 2. BRCA2 Hap1–mES direction and covariance summaries

| Summary | Groups | Concordant | Estimate | 95% interval |
| --- | --- | --- | --- | --- |
| Primary R (HAP1 extremes) | 388 | 193/388 | 49.74% | 44.79–54.76% |
| R, size 2 | 341 | 171/341 | 50.15% | — |
| R, size 3 | 47 | 22/47 | 46.81% | — |
| R, fixed first two HGVS keys | 388 | 196/388 | 50.52% | — |
| R, reverse mES extremes | 388 | 193/388 | 49.74% | — |
| Pooled covariance ratio Φ | 388 | — | 87.51% | 76.80–98.66% |
| Fractional-rank Φ | 388 | — | 95.985% | — |
| Class-control Φ | 388 | — | 87.59% | — |

Primary R: position-cluster bootstrap 95% CI, 44.79–54.76%; Wilson 95% CI, 44.79–54.70%. No within-group score ties, zero selected contrasts, or endpoint ties; conditional tie-adjusted benchmark 50.00%. Φ is a scale-dependent covariance ratio, not a measurement-noise fraction. The class-control estimate of 87.59% comes from a separate pipeline. The main Φ interval is the ordinary group-bootstrap 95% interval (76.80–98.66%).

## Table 3. MSH2 primary 6-TG/mock held-replicate transfer

| Held replicate | Complete groups | DNA keys | Positions | AA-only RMSE | DNA-informed RMSE | Relative RMSE gain (95% interval) |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | 8,547 | 24,619 | 887 | 2.310 | 2.213 | 4.18% (3.79–4.58%) |
| 2 | 8,416 | 23,997 | 900 | 2.109 | 2.015 | 4.48% (4.03–4.94%) |
| 3 | 8,474 | 24,156 | 892 | 2.226 | 2.133 | 4.20% (3.75–4.65%) |

Primary endpoint: 6-TG/mock. Intervals use 5,000 protein-position-bootstrap draws with β refit in each fold and draw. This is a shared-pool, within-study held-culture transfer; source scores at evaluated positions remain available, and only β is cross-fitted by position. All endpoint RMSE and paired ΔMSE values are in Table S5a; MAE values are in Table S5b.

## Table 4. Rule-specific selection summaries

| Dataset / rule | R (%) | Selected groups | Empirical 2.5–97.5% range |
| --- | --- | --- | --- |
| Gaussian null, A: evaluation-dependent | 99.94 | 96.063 per dataset | 98.97–100.00 |
| Gaussian null, B: fixed threshold | 50.02 | 22.935 per dataset | 30.43–70.83 |
| Gaussian null, C: cross-fitted threshold | 49.61 | 11.631 per dataset | 20.00–80.00 |
| BRCA2, A: evaluation-dependent | 100.00 | 77 | — |
| BRCA2, B: fixed laboratory variance | 56.06 | 66 | — |
| BRCA2, B: fixed technical variance | 50.46 | 327 | — |
| BRCA2, C: position cross-fitting | 61.03 | 11.696 per nonempty split | 37.50–81.82 |

Null ranges describe dataset-to-dataset variability; the BRCA2 C range describes split variation. These empirical ranges are not sampling confidence intervals. Real-data rules select different groups and do not provide a matched-population performance comparison.
