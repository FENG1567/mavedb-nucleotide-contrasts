# Third-party sources and license boundaries

The MIT license covers this team's original code. CC BY 4.0 covers this team's
new statistical summaries and annotations. Neither grant replaces upstream
rights. No external repository code has been copied into this package.

| Source | Version / attribution | Included material and terms |
|---|---|---|
| MaveDB BRCA2 HAP1 | Huang et al.; `urn:mavedb:00001225-a-1` | Exact original score export in `data/external/brca2_cc0/`, CC0-1.0 as displayed on the source record. The mapped copy adds this team's sequence annotations; underlying scores remain CC0. |
| MaveDB BRCA2 mES | Sahu et al.; `urn:mavedb:00001242-a-1` | Exact original score export and mapped copy, same source-specific CC0 boundary. |
| Ensembl BRCA2 CDS | ENST00000380152; versioned request ENST00000380152.8 | `data/external/brca2_cc0/ENST00000380152.txt` is a pinned reference sequence **excepted from the directory's CC0 label**. Reuse follows Ensembl data terms, not this team's MIT/CC BY grant. File SHA-256: `c01ae0690f9267aabaaa39bb29e1f5c5cf0dd05354b9ddf236b72a8cdc810b9a`. The returned FASTA identifier/version was not independently verified. |
| Jia et al., GEO GSE162130 | MSH2 deposited raw codon count table | Fetch-only, exact URL and hash in `data/external/sources.csv`. Source submitters retain their rights. The original raw table is not reissued under this team's data license. |
| BARD1 SGE study | `ivanw314/BARD1_SGE_analysis`, commit `72ba50c8480a8dd88528c811ccc906d41fb3d971` | No upstream repository license identified. Source tables and external code are not bundled. Fetch-only pinned URLs/hashes provided. Figure S2 uses this team's newly calculated covariance summaries with explicit source attribution. |
| VHL study | [PMC11250436](https://pmc.ncbi.nlm.nih.gov/articles/PMC11250436/), Table ST1 | Article states CC BY 4.0 except separately credited material. Only this team's new aggregate estimates are bundled. Source table retrieval and credit checking are separately described. |
| Other collection score sets | 92-set frozen MaveDB collection; 27 multi-key source sets | Only this team's identity/qualification/reconciliation outputs are bundled. No blanket CC0 assumption is made for all MaveDB records; raw and mapped source score exports are not bundled. Source IDs, hashes and download endpoints are recorded. |
| OpenSplice | [lehner-lab/OpenSplice](https://github.com/lehner-lab/OpenSplice), ENA PRJEB111846, Figshare DOI `10.6084/m9.figshare.32337414` | Upstream code MIT does not determine data terms. External code and original dataset are not bundled. Historical OpenSplice master request returned a 146-byte HTTP403 page; this is not supplied as a dataset. The relation of historical file ID 64886694 to the Figshare article was not verified. |
| COMPASS | GEO GSE307247 | Fetch-only source data. Bundled fixed-pair outputs describe technical coverage only; they do not establish biological event equivalence. |
| ClinVar | NCBI ClinVar, original 53-record list and 2026-11-23 coverage query | Bundled author-generated crosswalks/coverage summaries; no full ClinVar archive. Database state can change and current online responses need not equal the frozen snapshot. |
| RAD51D | Frozen 200,000-pair read-prefix feasibility summaries | Author-generated descriptive counts only. Raw FASTQ and prefix reconstruction are not redistributed or newly rerun in this package. |
| Hu / AVENGERS and other published pipelines | Cited for source-score interpretation | Referenced, not copied. Do not infer a software license from the paper's article license. |

## Attribution and transformations

Reused scores were mapped to protein consequences using the pinned reference
sequence; the team's computations apply source-specific standardization,
within-AA grouping, covariance, direction, resampling and prediction procedures.
Original score tables are not called raw sequencing data. New derived summary
files are distinct from those original tables. Refer to the manuscript's source
study citations in addition to citing this software release.

## Primary terms and records

- [MaveDB HAP1 record](https://www.mavedb.org/score-sets/urn:mavedb:00001225-a-1)
- [MaveDB mES record](https://www.mavedb.org/score-sets/urn:mavedb:00001242-a-1)
- [CC0 legal code](https://creativecommons.org/publicdomain/zero/1.0/legalcode)
- [Ensembl content and data terms](https://www.ensembl.org/about/articles/ensembl-content-disclaimer)
- [NCBI GEO disclaimer](https://www.ncbi.nlm.nih.gov/geo/info/disclaimer.html)
- [NCBI molecular data policies](https://www.ncbi.nlm.nih.gov/home/about/policies/)
- [ClinVar public data](https://www.ncbi.nlm.nih.gov/clinvar/docs/submission_get_started/)
- [Pinned BARD1 source commit](https://github.com/ivanw314/BARD1_SGE_analysis/commit/72ba50c8480a8dd88528c811ccc906d41fb3d971)
- [VHL article](https://pmc.ncbi.nlm.nih.gov/articles/PMC11250436/)
- [CC BY 4.0 legal code](https://creativecommons.org/licenses/by/4.0/legalcode)

The splicing variant coverage file has source PSI, delta-PSI, delta-logit and
pooled numeric values removed; only author-computed coverage/event diagnostics
and source identifiers remain. Group-level differences are new computed
contrasts. ClinVar classification, title, accession, review-status and phenotype
strings in the crosswalk and term-count table are explicitly attributed NCBI
source excerpts, excluded from the authors' CC BY grant; see LICENSE-DATA.md.

License evidence was checked on 2026-11-24. `data/external/sources.csv` describes
source acquisition, not a grant of additional permissions.
