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
| `mara/*.svg` | One readable MARA-style locus map per supported target, including a symbol key |
| `mara-table/*.svg` | Matching annotation table per locus, including evidence notes and a compact key |
| `proof/proof.json` | Machine-checkable Tn1/Tn2/Tn3 component, grammar, match and output verdicts |
| `proof/components.tsv` | One row per independently sequence-detected Tn1/Tn2/Tn3 component |
| `proof/matches.tsv` | One row per assembled and known-element-matched Tn1/Tn2/Tn3 locus |
| `proof/report.html` | Portable human-readable proof report linking all representations |
| `detectors/` | Present only when Matryoshka ran optional external detectors |

The proof status is `PASS` only when every required Tn1/Tn2/Tn3 component was
detected from sequence, the component order/orientation grammar is complete, a known
whole-element reference supports the name, and zero internal components were merely
projected. Inputs without a Tn1/Tn2/Tn3 locus are reported as `NO_TN123_LOCI`, not as
a false pass.

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

`validated` is the default and contains the expert-example-backed alpha subset. `all`
adds experimental and legacy references. The selected profile is embedded in both
`annotation.json` and `run.json`.

The diagrams are renderings of the annotation hierarchy, not independent evidence.
For Tn1/Tn2/Tn3, solid component arrows represent independent sequence matches used
by the assembly grammar; outlined dashed arrows represent components projected from a
curated map after a supported parent call. Use the matching table and machine-readable
annotation for identities, coverages, boundary evidence and grammar status.
