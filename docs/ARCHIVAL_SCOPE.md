# Archival versus raw-source reconstruction

The main inferential BRCA2 and MSH2 computations and the figure inputs are
explicitly separated from ancillary historical feasibility work.

## Collection

The frozen qualification tables under
`data/derived/revision_20261122/dna_audit/` contain the team's computed mappings,
identity flags and target checks. They support traceable reconstruction of
the published tiers. They are not new independent molecular measurements.

`collection_source_exports.csv` lists archived raw score inputs, original
hashes and MaveDB API endpoints. Only the two main BRCA2 exports have been
individually cleared for bundling. The 27-set source-to-mapped numeric audit
also requires original mapped tables and pinned reference evidence; those
third-party score copies are not bundled under a blanket project license.
The full-source audit must be given a reconstructed source project layout.
Frozen-result checks do not repeat its historical 27-set source audit.

## Supplementary checks

- Figure S2 can be rebuilt from the two supplied author-generated summary
  tables. BARD1 original inputs are fetched from a pinned commit. VHL Table ST1
  is obtained from the article. Original count-level reconstruction is a
  separate path and is not implied by successful summary-level plotting.
- The splicing outputs retain a fixed universe of 22 pairs. The historical
  master download was an HTTP403 page, not usable data. We preserve the stated
  coverage limitations; no fabricated master file replaces it.
- ClinVar coverage is time-dependent. The original 53-record snapshot and the
  34 records returned by a later ESummary query reflect different fields and
  dates. Bundled crosswalks document those observations.
- RAD51D outputs are descriptive 200,000-read-pair-prefix feasibility summaries,
  not unique molecule counts or new validation experiments. Raw FASTQ processing
  was not newly executed for this release.

These boundaries prevent a frozen plot or an archival output from being
misrepresented as a newly validated end-to-end data pipeline.
