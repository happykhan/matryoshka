# Result-directory contract

`matryoshka run INPUT --out DIRECTORY` is the stable external-alpha workflow.

## Files

| Path | Contract |
|---|---|
| `annotation.json` | Full v1 annotation document validated by `schema/annotation-v1.schema.json` |
| `annotation.gff3` | One GFF3 document with sequence regions, unique stable IDs and `Parent` links |
| `run.json` | Compact counts and relative output locations |
| `mara/*.svg` | One readable MARA-style locus map per supported target, including a symbol key |
| `mara-table/*.svg` | Matching annotation table per locus, including evidence notes and a compact key |
| `detectors/` | Present only when Matryoshka ran optional external detectors |

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
Use the matching table and machine-readable annotation to distinguish directly
detected parent loci and boundary evidence from components projected from a curated
reference map.
