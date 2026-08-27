# Changelog

All notable changes to Matryoshka are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- Rich staged console reporting for the one-command workflow, including detector
  status, reference-scan progress, per-contig results and a final summary; `--quiet`
  suppresses presentation output for automated use.
- `matryoshka expert-rules` Markdown/HTML export and matching files in every result
  bundle. The report translates the executable YAML into plain-language component,
  grammar, classification, exact-name, variant and fragment decisions for review.
- A versioned 17-target priority-transposon roadmap with machine-readable YAML,
  Markdown/JSON export, explicit component-rule evidence levels and a shared
  validation contract.
- BactScout-style detector orchestration: `matryoshka run` now executes supported
  AMRFinderPlus, ISEScan and IntegronFinder environments by default, retains raw
  outputs and records completed/provided/unavailable/failed/disabled provenance.
- `matryoshka preflight` and Pixi `preflight`/`analyse` tasks for a one-FASTA,
  one-command workflow.
- Pixi-managed AMRFinderPlus on Intel macOS and Apple silicon through the osx-64
  Bioconda build/Rosetta, in addition to the Linux detector environments.
- First-run AMRFinderPlus database initialisation when Bioconda has installed the
  executable but no valid default database; existing databases are left untouched.
- AMRFinderPlus software and database versions in the detector provenance stored in
  `annotation.json` and `run.json`.
- Reviewable component grammars for Tn21, Tn1721 and Tn1722, with independent
  component scans and `*-like` assembly under the whole-locus-free
  `component-rules` profile.
- GenBank (`annotation.gbk`) in every standard result bundle, alongside GFF3.
- A reproducible 16-record NCBI Tn1 accession panel and conservative result
  ledger covering complete loci, large insertions, fragments and mosaics.
- A self-contained Word/PDF Tn1/Tn2/Tn3 evidence report with embedded locus
  maps, locus tables, hierarchy views and component ledgers for every reviewed
  accession, plus natural pEK499 fragments and the arbitrary-contig test.
- A reproducible pEK499 output bundle retaining expert-reviewed two natural
  Tn2-related regions as partial fragments.
- Reviewable, versioned YAML expert definitions for Tn1/Tn2/Tn3 type,
  subtype, variant and fragment rules, with JSON/Markdown export, a public
  JSON Schema and a generic definition loader.
- Reviewed whole-locus definitions for Tn2c (HM749967), Tn2.1 (CP028717),
  Tn1Mer (GQ160960) and the legacy 9 bp-duplication Tn3 sequence (V00613).
- A reproducible 14-case validation ledger and review report covering exact
  real accessions, pEK499 fragments, minor variation, truncation and related
  elements that must not be misnamed as Tn1/Tn2/Tn3.
- A complete reviewed-definition demonstration bundle containing every locus
  locus map and table, plus a `run.json` index with direct paths for all maps.
- A source-backed locus component ontology with 30 raw component classes,
  18 compound-element assembly grammars, named validation examples and
  `matryoshka catalog` JSON/TSV export.
- Locus-map symbols for cassette remnants/arrays, attI/Pc, ISCR ori/ter sites,
  captured segments, direct-repeat features, replicons, ncRNA, broad resistance
  regions and explicit unknown fragments/gaps.
- BLAST reference scan against bundled MGE references (Tn4401a/b, Tn1546, Tn21,
  Tn1331, Tn5393, Tn6022/Tn6019/Tn6172, ISEcp1, Tn7, Tn552, In2, GIsul2,
  PlasmidFinder Enterobacteriales replicons, Tn3 *res* consensus, rolling-circle
  ter-site motifs).
- Declarative rule tables for flanked-cargo transposons, one-ended captures,
  rolling-circle captures, and signature-based detection.
- Tn4401 variant discrimination (a vs b) via subject-coverage + HSP-count scoring.
- Plasmid replicon typing via PlasmidFinder DB.
- Per-feature confidence scoring (0.0–1.0 + high/medium/low/speculative label).
- Single-IS26 translocatable unit detection (low-confidence).
- Tn3-family *res*-site positional annotation.
- Acinetobacter resistance island compound detection (AbGRI1, AbaR3).
- Multi-contig FASTA support with per-contig output.
- GenBank flat-file output via Biopython.
- Scale-accurate linear SVG visualisation with auto-scaling canvas.
- A one-command `matryoshka run` workflow producing versioned JSON, valid
  multi-record GFF3, run metadata and paired locus map/table directories.
- A Tn1/Tn2/Tn3 acceptance-proof bundle linking every independent component
  call to the assembly grammar, known-element match, locus outputs, original
  nested hierarchy and CellGen/Wolvercote cell-format representation. The
  bundle includes JSON verdicts, TSV ledgers and a portable HTML report.
- Deterministic feature IDs, a public v1 JSON Schema and embedded tool,
  reference-database, input-checksum and detector provenance.
- A checksummed machine-readable reference manifest with validated and
  experimental profiles.
- Concurrent reference scanning with user-controlled workers and progress.
- Public NCBI acceptance plasmids, compact detector parser fixtures, a locked
  cross-platform Pixi environment and a core Docker image.
- locus single-line SVG output with directional genes, IR flags, exact
  `res` sites, long-gap notation, hollow expected/filled confirmed TSD lollipops
  and fragment edges. Every map now includes a self-contained symbol key.
- locus tabular SVG output with position, hierarchical feature name and
  strand glyph, accession/FID, type, evidence notes and a compact symbol key.
- Canonical Tn1 (NC_008357), Tn2 (AY123253) and Tn3 (HM749966) detection with
  curated internal maps and conservative `Tn1/2/3` fragment reporting.
- Locus-based Tn1/Tn2/Tn3 detection in arbitrary-length sequences, including
  collinear multi-HSP assembly, nearest-reference `Tn*-like` calls, reverse
  orientation, internal insertion/deletion status and separate views for
  multiple loci on the same contig.
- Independent Tn1/Tn2/Tn3 component scanning for both terminal IRs, `blaTEM`,
  `tnpR`, `res` and `tnpA`, followed by orientation-aware grammar assembly.
  Exact member names now require the component grammar as well as whole-locus
  reference evidence; an interrupted `tnpA` can be assembled from collinear HSPs.
- Reference-only annotation, so the bundled BLAST scan can run without dummy
  detection-tool output files.
- Seven curated complete ISEcp1-blaCMY transposition-unit references from
  expert-reviewed supplied corpus, with fragment-safe reporting and locus views.
- Exact ISKpn6, ISKpn7 and IS1999 references for structural rules that cannot
  safely rely on an IS-family label.
- A paper- and expert-feedback-backed locus annotation specification covering the
  feature database, annotation grammar, symbols and acceptance fixtures.
- Curated component-aware locus maps for expert-reviewed Tn21, Tn1721, Tn1722, Tn4401,
  Tn5393 and Tn5403 references, with stable FIDs and exact IS names.
- Collinear assembly for minor/indel variants of those curated unit references;
  complete parents override competing terminal fragments.
- Whole-cassette reconstruction from IntegronFinder ORF plus oriented attC
  evidence, with parser children retained through the CLI hierarchy rebuild.
- Component-based class-1-integron and complete Tn402 inference using exact
  5'-CS, 3'-CS, IRi, IRt and tni references.
- Boundary-adjacent TSD detection with exact repeat coordinates, evidence
  strength and negative controls for random/non-adjacent matches. A small
  coordinate-refinement window opens only when independent terminal-repeat
  evidence supports it.
- Feature-specific TSD lengths and explicit boundary evidence states for
  sequence-confirmed, searched-but-not-found and unavailable flanking sequence.
- Hierarchy enforces biologically-sensible parent/child pairings.
- CHANGELOG, LICENSE, CITATION.cff, CONTRIBUTING.md, references manifest.

### Fixed
- Locus-map labels now use collision-aware lanes instead of overlapping at
  dense sites; gene arrows use the grey evidence palette and unresolved or
  inserted sequence remains visibly distinct.
- Tn1-family assembly now evaluates compatible terminal-IR pairs across a
  declared maximum span, recovering the large Tn1Mer insertion while rejecting
  redundant larger boundaries around a tighter complete locus.
- GFF3 output now uses one stable writer for both `run` and `annotate`, avoiding
  divergent escaping and duplicate feature IDs; GenBank gene/IR types are
  emitted as CDS/repeat features rather than generic miscellaneous features.
- Non-canonical subtype insertions can no longer create stray Tn1/Tn2/Tn3
  fragment calls. Partial-family evidence is assessed against canonical
  references, while subtype references require substantially complete loci.
- The default validated profile now detects the complete the expert reviewer-supplied Tn21
  parent, rather than reporting only its internal In2 integron.
- Fresh installations no longer depend on a missing sibling `cell-format`
  checkout. Unsupported legacy circular SVG/PNG choices are no longer
  advertised; the supported linear, locus-map and locus-table SVG renderers remain bundled.
- Oppositely oriented or merely family-level IS6 calls are no longer merged
  into a generic "IS26 island". Only explicitly named, complete, directly
  oriented IS26 pairs emit individual pseudo-compound candidates, with outer
  direct-repeat evidence recorded separately.
- Hierarchy algorithm rewritten from O(n³) to O(n²) with clean semantics.
- ISEcp1 cargo pattern no longer matches blaTEM via "LAT" substring.
- ISEcp1 cargo inference now follows the IRR-facing side of the IS orientation
  and remains boundary-incomplete until IRalt/ter and the outer DR are confirmed.
- Variant clustering keeps multi-FASTA records separate even when every record
  begins at coordinate 1.
- The bundled ISEcp1 reference is now the exact, expert-reviewed 1,656 bp ISfinder sequence;
  the previous 166 kb accession-scale sequence could mislabel an entire region
  as one insertion sequence.
- The In2 reference is now the exact 11 kb IRi-to-IRt region from Tn21 rather
  than the previous 53 kb accession-scale sequence, preventing a tni fragment
  from being promoted to a complete integron.
- Generic `blaOXA` calls are explicitly reported as unresolved family-level
  evidence and down-weighted until an exact allele is available.
- Locus maps and tables distinguish independently sequence-detected components
  (solid arrows) from reference-projected components (outlined dashed arrows),
  and describe the component grammar and boundary evidence in their keys.
- Wolvercote parser no longer crashes on gene names containing parentheses.
- Wheel builds no longer add the bundled reference directory twice.

## [0.1.0a1] — Planned external alpha
Initial public release.
