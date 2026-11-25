# Validation record

This record describes what was actually checked for local candidate
`SRR-20261124-C1-reproducibility`. It is separate from the manuscript and
contains no reviewer correspondence.

## Completed checks

| Check | Environment | Result |
|---|---|---|
| MSH2 raw-count reconstruction | `ubuntu-server`, Python 3.12.3, one thread | PASS: 9 primary rows and 45 sensitivity rows matched frozen v4 outputs at absolute tolerance 1e-10. The three primary 6-TG/mock bootstrap intervals used 5,000/5,000 draws with seeds 20261125, 20261126 and 20261127. Input SHA-256 matched `9e315efaf715db4f906aaf5b82b4175847af819d5f65659f03e5070d710e2299`. |
| BRCA2 raw-to-mapped check | `ubuntu-server`, Python 3.12.3 | PASS: both source exports were mapped offline with the pinned CDS and every archived mapped field matched. |
| BRCA2 core result | `ubuntu-server`, Python 3.12.3 | PASS: 388 groups, 823 shared DNA keys, Phi = 0.8750878921486126, and the ordinary group-bootstrap interval matched the frozen record. |
| BRCA2 direction/tie result | `ubuntu-server`, Python 3.12.3 | PASS: primary R = 193/388 = 49.74%; published group fields and estimator summaries matched. Two-member tie-bound diagnostics were corrected to carry observed pair concordance. |
| Collection reconciliation | `ubuntu-server`, Python 3.12.3 | PASS: 92 manifest sets, 89 analyzable sets, 13,041 historical multi-key groups, 31,723 keys, and all 27 raw-to-mapped numeric multiset checks passed in the original source tree. |
| Frozen aggregate checks | `ubuntu-server`, Python 3.12.3 | PASS: collection tiers, BRCA2 aggregates, MSH2 cohorts/gains and 22 splicing groups. |
| Figure generation | `ubuntu-server`, Matplotlib Agg, Python 3.12, one thread | PASS: four main PDFs and Supplementary Figures S1–S3 were generated with PDF signatures and source tables by the single `figures` entrypoint. |
| Git byte policy | local candidate | PASS after `.gitattributes` was changed to `* -text`; CRLF bytes of the two bundled MaveDB exports are preserved. The release manifest must also pass in a clean checkout. |

## Not claimed as complete

The optional BARD1/VHL source-count rerun is not required for submitted Figure
S2 because verified author-generated summary tables are bundled. The BARD1
portable source analysis started successfully; the VHL attempt was blocked when
the temporary supplementary ZIP disappeared from the server workspace. Its
official Europe PMC URL, bundle hash and member hash remain in the source
manifest, and the optional script accepts a user-supplied ZIP.

The full 92-set raw collection audit is reproducible only when the original
score exports and mapped tables are supplied; those third-party files are not
redistributed under a blanket project license. Frozen qualification tables and
the 27-set reconciliation were checked in the original source tree. The
OpenSplice master download was an HTTP 403 response, so no biological splicing
validation is claimed. ClinVar summaries are frozen coverage records, not a
timeless online requery.

The public repository is https://github.com/FENG1567/mavedb-nucleotide-contrasts; the manuscript release is v1.0.0.
No DOI or Zenodo record has been assigned. The historical checks above retain
their stated scope. In a clean checkout, run `tests/verify_manifest.py` and
`tests/verify_frozen_results.py`; the manifest excludes separately acquired
raw-source caches. Source-specific license boundaries are in `THIRD_PARTY.md`.
