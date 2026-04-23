# Matryoshka

**Nested mobile genetic element annotation for bacterial genomes.**

Matryoshka takes a resolved bacterial assembly plus the outputs of standard MGE detection
tools (ISEScan, AMRFinder+, IntegronFinder) and produces a *nested* annotation — showing
what is contained inside what — rather than a flat coordinate list.

It runs rule-based inference for known composite transposons, BLAST against a bundled
reference library for structurally conserved elements, and confirms boundaries with
target-site-duplication (TSD) and inverted-repeat (IR) detection. Every feature is
tagged with a confidence score.

---

## What it identifies

| Category | Elements |
|---|---|
| Insertion sequences | IS6 (IS26, IS257), IS1380 (ISEcp1), IS30 (ISApl1), IS91/ISCR (rolling-circle), plus anything ISEScan calls |
| Composite transposons | Tn4401 (blaKPC, variant a/b discriminated), Tn1999 (blaOXA-48), Tn_mcr1 (ISApl1+mcr-1), IS26 islands (merged pairs), IS26 TUs (single-IS, low-confidence) |
| Tn3-family unit transposons | Tn1546 (vanA), Tn21, Tn1331, Tn5393 — via BLAST or cargo signature |
| One-ended / rolling-circle capture | ISEcp1 → blaCTX-M / CMY / ACC / LAT / DHA; ISCR1-mediated + IS91 capture units |
| Genomic islands | AbGRI1 (Tn6022 + Tn6172 flanking comM), AbaR3-like (Tn6019), GI*sul2* |
| Integrons | Class 1 / 2 / 3 with cassettes (via IntegronFinder) |
| Plasmid features | Replicon typing via PlasmidFinder (Enterobacteriales, 159 replicons) |
| Structural sub-features | Tn3-family *res* sites (positional), ter-site motifs for rolling-circle elements |
| AMR cargo | Everything AMRFinder+ calls |

## How the hierarchy works

Each feature becomes a child of the *smallest* feature that strictly contains it. IS
elements can't parent AMR / cassette / integron features (a spatial overlap doesn't
imply biological containment), but transposons and integrons can. The result is a
tree of `MGEFeature` nodes you can walk, render, or flatten.

## Example — pKpQIL (K. pneumoniae, blaKPC plasmid)

```
Tn4401a              8882–18788   (BLAST, 99.84% ID, 100bp Pout deletion, conf=1.0)
├── res_Tn4401a      9082–9202    (positional, conf=0.3)
├── IS21_203         13998–15928
└── blaKPC-3         15984–16862

IS26_island_1        25540–32163
├── IS6_292          25540–26359
├── blaOXA           26527–27348
├── blaTEM-1         28051–28908
└── IS6_292          29089–32163

IncFIB(pQil)         53120–53859  (PlasmidFinder, 100% ID, conf=1.0)
```

---

## Install

Requires: [pixi](https://pixi.sh). All other dependencies (Python 3.11+, BLAST+,
matplotlib, biopython) are managed by pixi.

```bash
git clone https://github.com/happykhan/matryoshka
cd matryoshka
pixi install
```

Detection tools (ISEScan, AMRFinder+, IntegronFinder) live in isolated pixi
environments and are invoked as subprocesses. They are optional — you can also pass
existing outputs on the command line.

## Quickstart

```bash
# Annotate a plasmid using pre-computed detection-tool outputs
pixi run python -m matryoshka annotate plasmid.fasta \
    --isescan     isescan_output.tsv \
    --amrfinder   amrfinder_output.tsv \
    --integrons   integronfinder_output.integrons \
    --format      json \
    --out         plasmid.annotation.json

# Render a scale-accurate linear SVG
pixi run python -m matryoshka annotate plasmid.fasta \
    --isescan isescan_output.tsv --amrfinder amrfinder_output.tsv \
    --format linear --out plasmid.svg
```

Multi-FASTA input (e.g. chromosome + plasmids) is handled automatically — each
contig gets its own hierarchy. For per-contig SVG/PNG output pass `--out some_dir/`.

## Output formats

| `--format` | Contents |
|---|---|
| `json` | Full hierarchy with attributes, confidence, TSD/IR evidence |
| `gff3` | One row per feature with Parent= relationships, IR/TSD/source attributes |
| `genbank` | Biopython SeqFeature round-tripped record — loadable in Artemis / SnapGene |
| `wolvercote` | Compact nested cell-format notation (see [cell-format](https://github.com/happykhan/cell-format)) |
| `svg` | Circular schematic via wolvercote renderer |
| `png` | Rasterised circular diagram |
| `linear` | Scale-accurate linear map (matryoshka's flagship output) |

## Reference library

Bundled in `matryoshka/references/`. Every reference is traceable to a published
accession — see the header line of each FASTA for origin.

Refresh or extend via the fetcher scripts:

```bash
NCBI_EMAIL=you@example.org pixi run python scripts/fetch_mge_references.py --force
```

References include Tn4401 variants (a/b), Tn1546, Tn21, Tn1331, Tn5393,
Acinetobacter islands (Tn6022, Tn6019, Tn6172, Tn6022Δ), GI*sul2*, ISEcp1, Tn7,
Tn552, In2, mcr-1 exemplars, Tn3 *res* consensus, rolling-circle ter-site motifs
and PlasmidFinder's Enterobacteriales replicon DB (159 records).

## Validation plasmids

| Plasmid | Accession | Tests |
|---|---|---|
| pKpQIL | NC_014016 | Tn4401a detection, variant discrimination, IncFIB(pQil) replicon |
| pOXA-48a | JN626286 | Tn1999 inference, TSD=TGCTG confirmation, IncL replicon |
| pEK499 | EU935739 | IS26 island merging (16 children), nested ISEcp1_capture, IncFIA+IncFII |

Run the test suite:

```bash
pixi run pytest tests/ -q
```

## Confidence scoring

Every feature is tagged with a confidence score in `attributes.confidence` (0.0–1.0)
and a human label in `attributes.confidence_label`:

| Score | Label | Typical evidence |
|---|---|---|
| 1.00 | high | BLAST hit ≥98% identity and ≥95% subject coverage |
| 0.85 | high | Rule-inferred composite transposon with TSD and IR |
| 0.70 | medium | Rule-inferred transposon, no boundary evidence |
| 0.60 | low | IS element with IR |
| 0.40 | speculative | Single-IS26 translocatable unit |
| 0.30 | speculative | Positional res site (not sequence-confirmed) |

See [`matryoshka/confidence.py`](matryoshka/confidence.py) for the full rubric.

## Status

Alpha. Core annotation pipeline is stable and validated on three reference
plasmids. Known coverage gaps vs the Partridge 2018 review are tracked in
[GAPS.md](GAPS.md).

## Citing

If you use Matryoshka in published work, please cite the Partridge 2018 review
whose framework it implements:

> Partridge SR, Kwong SM, Firth N, Jensen SO. *Mobile Genetic Elements Associated
> with Antimicrobial Resistance*. Clin Microbiol Rev. 2018;31(4):e00088-17.
> doi:10.1128/CMR.00088-17. PMCID: [PMC6148190](https://pmc.ncbi.nlm.nih.gov/articles/PMC6148190/).

A Matryoshka-specific citation will follow once the preprint is available.

## License

MIT. See [LICENSE](LICENSE).

## Contributing

Contributions welcome — particularly new transposon inference rules or
reference-sequence additions. See [GAPS.md](GAPS.md) for open work.

New transposon rules are declarative (see `FLANKED_RULES`, `ONE_ENDED_RULES`,
`SIGNATURE_RULES` in `matryoshka/transposon.py`); adding most Tn families is a
one-line addition to a dict.
