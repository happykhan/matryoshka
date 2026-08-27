# Result-directory contract

`matryoshka run INPUT --out DIRECTORY` is the stable external-alpha workflow.

## Files

| Path | Contract |
|---|---|
| `annotation.json` | Full v1 annotation document validated by `schema/annotation-v1.schema.json` |
| `annotation.gff3` | One GFF3 document with sequence regions, unique stable IDs and `Parent` links |
| `annotation.cell` | Nested CellGen/Wolvercote cell-format representation of the same hierarchy |
| `run.json` | Compact counts and relative output locations |
| `hierarchy/*.svg` | One whole-record, scale-accurate rendering of the original nested hierarchy |
| `locus-map/*.svg` | One readable locus map per supported target, including a symbol key |
| `locus-table/*.svg` | Matching annotation table per locus, including evidence notes and a compact key |
| `proof/proof.json` | Machine-checkable Tn1/Tn2/Tn3 component, grammar, match and output verdicts |
| `proof/components.tsv` | One row per independently sequence-detected Tn1/Tn2/Tn3 component |
| `proof/matches.tsv` | One row per assembled Tn1/Tn2/Tn3 locus, including the rule-based type call, optional reference context, definition ID/version and sequence-difference counts |
| `proof/report.html` | Portable human-readable proof report linking all representations |
| `detectors/` | Present only when Matryoshka ran optional external detectors |

`run.json` contains a `locus_outputs` array with one entry per generated
locus map. Each entry records the sequence record, call, family, coordinates,
viewport and direct relative paths to its locus map, locus table and hierarchy
view. Consumers should use this index rather than guessing filenames.

The proof status is `PASS` only when every required Tn1/Tn2/Tn3 component was
detected from sequence, the component order/orientation grammar is complete, the
expert component profiles support a type, and zero internal components were merely
projected. A known whole-element reference is secondary evidence and is not required
for `PASS`. Inputs without a Tn1/Tn2/Tn3 locus are reported as `NO_TN123_LOCI`, not
as a false pass.

An expected incomplete locus with sequence-detected components and canonical-family
homology receives a `PARTIAL` locus verdict. A run containing only such loci reports
`PARTIAL_TN123_EVIDENCE`; a run containing complete passes as well reports
`PASS_WITH_PARTIAL_TN123_EVIDENCE`. These are conservative biological outcomes, not
pipeline failures.

Reviewed subtype context, such as the ISEcp1-associated insertion in Tn2.1, is
counted separately from projected required components. It may enrich a supported
figure but cannot satisfy the six-part Tn1/Tn2/Tn3 grammar.

## Coordinates and identity

Coordinates are 1-based and inclusive. Every JSON feature has a deterministic `id`;
children name that value in `parent_id`. The same biological call on the same record,
with the same hierarchy and order, receives the same ID on repeated runs. Record names
are part of the identity, preventing collisions across multi-FASTA inputs.

The JSON envelope records:

- output-schema version;
- Matryoshka version;
- reference-database version and selected profile;
- input path and SHA-256 checksum;
- command, Python and platform provenance;
- paths to any supplied or generated detector outputs.

Changing a documented schema incompatibly requires a new schema version. New optional
fields may be added within v1. Consumers should ignore unknown fields and must not infer
biological confidence solely from the feature name.

## Profiles

`validated` is the default and contains the expert-example-backed alpha subset. The
`tn123-components` profile runs the Tn1/Tn2/Tn3 expert discovery path without any
complete Tn1/Tn2/Tn3 reference lookup. `all` adds experimental and legacy references.
The selected profile is embedded in both `annotation.json` and `run.json`.

The diagrams are renderings of the annotation hierarchy, not independent evidence.
For Tn1/Tn2/Tn3, solid component arrows represent independent sequence matches used
by the assembly grammar; outlined dashed arrows represent components projected from a
curated map after a supported parent call. Use the matching table and machine-readable
annotation for identities, coverages, boundary evidence and grammar status.
