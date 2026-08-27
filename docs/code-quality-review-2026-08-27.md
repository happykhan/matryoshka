# Code-quality review — 2026-08-27

Scope: the complete Python package, command-line workflows, output writers,
definition loaders, tests, packaging and the Tn1/Tn2/Tn3/Tn21/Tn1721/Tn1722
component-rule path.

## Outcome

No release-blocking correctness issue remains in the reviewed path. Static checks
and the complete regression suite pass. The review found two behavioural defects
and three maintainability problems that were fixed in this change.

| Severity | Finding | Resolution |
| --- | --- | --- |
| P1 | A fixed 2.5 kb `tnpA` neighbourhood rejected complete loci with large internal insertions, including Tn1Mer. | Candidate assembly now evaluates terminal-IR pairs within a declared 60 kb maximum span. A regression test covers Tn1Mer. |
| P1 | Repeated homologous IRs could produce both a tight complete locus and a redundant larger parent. | Complete candidates are boundary-ranked and containing duplicates are suppressed. The accession panel exercises this on KX156773.1. |
| P1 | Tn21/Tn1721/Tn1722 maps were produced from whole-reference projection rather than independent component grammar. | Component profiles, counts, order, orientation, span and the Tn1721 duplication junction are now external YAML rules; `component-rules` runs without whole-element lookup. |
| P2 | `run` and `annotate` had separate GFF3 implementations with different escaping and identifier rules. | Both paths now use the stable multi-record writer. |
| P2 | GenBank was available only through `annotate`, and new gene/IR component types degraded to `misc_feature`. | Every run writes `annotation.gbk`; gene components are CDS and IR/direct-repeat components are repeat features. |
| P2 | Dense locus labels shared a fixed baseline and overlapped. | A four-lane collision allocator places feature labels deterministically; a regression test protects the layout. |
| P2 | Reviewed unit maps were embedded as a large Python constant. | Rules and feature maps now live in `unit_definitions.yaml`; the Python adapter is generic. |

## Structural audit

The largest modules remain `reference_scan.py` (about 1,660 lines),
`locus_map.py` (about 840 lines) and `__main__.py` (about 790 lines). They are
large but their current responsibilities are covered by tests and no new
family-specific biological conditions were added to them. A later internal
refactor should separate BLAST execution/HSP assembly from reference-family
adapters, and split map primitives from document composition. That refactor is
not required to add or review another biological definition: the extension
surface is YAML plus canonical component sequences.

## Strict second pass: one-command workflow and priority roadmap

The detector-orchestration change received a second maintainability audit with a
McCabe threshold of 12. The implementation was restructured before acceptance:

- optional detector execution, state transitions and best-effort/strict failure
  policy live in `detector_workflow.py`, not in the already large CLI module;
- detector runtimes and versions are typed records rather than loosely shaped
  dictionaries; AMRFinderPlus software/database versions and portable output paths
  are part of the result provenance;
- overlapping reference/rule/curated calls moved out of `__main__.py` into a focused
  deduplication policy, reducing the CLI module from roughly 850 to 750 lines and
  removing its only McCabe violation;
- expert-definition validation was split into named type/unit/component validators,
  so extending YAML does not add another condition chain to a monolithic loader;
- IR/TSD boundary confirmation and locus-table evidence notes were decomposed into
  small evidence-specific helpers without changing their tested behaviour.

Four pre-existing functions remain above the deliberately strict complexity threshold:
the IntegronFinder parser, locus-map document composer, GenBank feature converter,
and Tn1/Tn2/Tn3 insertion extractor. The map composer is the largest concern at 841
lines/module and complexity 32; it should be separated into a typed feature partition
and track renderers before substantially expanding the drawing vocabulary again. These
are maintainability debt, not hidden biological coverage, and are not blockers for the
isolated detector/roadmap change.

## Verification

- `ruff check matryoshka tests scripts`: pass
- `mypy matryoshka`: pass
- `pytest -q`: pass (188 tests)
- Python byte-code compilation: pass
- `git diff --check`: pass
- Component-only exact controls: Tn21-like, Tn1721-like and Tn1722-like
- Related negative controls: Tn4401, Tn5393 and Tn5403 are not promoted to those families
- Variation control: mutated Tn1722 plus a 120 bp insertion remains Tn1722-like
- Real-accession panel: six complete Tn1-like calls; incomplete records retain
  component evidence without a forced name
