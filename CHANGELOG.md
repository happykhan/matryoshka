# Changelog

All notable changes to Matryoshka are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
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
- Scale-accurate linear SVG/PNG visualisation with auto-scaling canvas.
- TSD detection with ±3bp boundary offset tolerance.
- Hierarchy enforces biologically-sensible parent/child pairings.
- CHANGELOG, LICENSE, CITATION.cff, CONTRIBUTING.md, references manifest.

### Fixed
- Hierarchy algorithm rewritten from O(n³) to O(n²) with clean semantics.
- ISEcp1 cargo pattern no longer matches blaTEM via "LAT" substring.
- Wolvercote parser no longer crashes on gene names containing parentheses.

## [0.1.0] — Planned first tag
Initial public release.
