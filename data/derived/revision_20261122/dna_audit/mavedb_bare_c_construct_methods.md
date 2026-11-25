# MaveDB bare-c target-construct validation

## Scope and coordinate rule

This check reads the previously retrieved 11 MaveDB score-set metadata JSON files and the existing `dna_variant_details.csv.gz`; it performs no new metadata download. It evaluates simple exonic `c.NREF>ALT` substitutions against the DNA string supplied in each record's `targetGenes[].targetSequence.sequence`. For this construct-level check only, it assumes `c.1` corresponds to the first supplied target-sequence base and translates from that base in the forward frame. It records the external identifiers and offsets verbatim but does not apply them to c. coordinates.

The official [MaveDB API model](https://github.com/VariantEffect/mavedb-api/blob/bfc04148cc5841c9ff63d92be752e430809c36b5/src/mavedb/models/target_gene_external_identifier.py) describes the offset as the target's offset within an external identifier. The [view model](https://github.com/VariantEffect/mavedb-api/blob/bfc04148cc5841c9ff63d92be752e430809c36b5/src/mavedb/view_models/external_gene_identifier_offset.py) requires a nonnegative integer, while the [target sequence model](https://github.com/VariantEffect/mavedb-api/blob/bfc04148cc5841c9ff63d92be752e430809c36b5/src/mavedb/models/target_sequence.py) stores sequence type and sequence. Those schema sources do not specify offset indexing or guarantee that c. numbering begins at target-sequence base 1. Therefore the results below test an explicit construct assumption and do not create an accession-qualified transcript/genome mapping.

## Validation

For each simple c. SNV, the script checks range, reference base, reference and alternate codons, DNA-derived amino-acid position/substitution, then compares that translation with the supplied p.HGVS. Complex, indel, intronic, and otherwise unsupported notations are counted as unassessed by this SNV check. No filtering on measured scores or functional outcomes was used.

See `mavedb_bare_c_construct_by_set.csv` for the per-score-set summary and `mavedb_bare_c_construct_validation.csv` for row-level details. The JSON report records target sequence hashes, metadata JSON hashes, schema references, and the interpretation boundary.
