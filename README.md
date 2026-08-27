# Matryoshka

Nested mobile genetic element annotation and MARA-style visualisation for bacterial
assemblies.

Matryoshka detects known mobile-element structures, reconstructs their biological
containment hierarchy, and writes machine-readable annotations plus one readable
diagram per locus. The external-alpha default is deliberately conservative: it scans
the reference sets validated against Sally Partridge's supplied worked examples.

## Try it

Matryoshka needs Python 3.11–3.13 and NCBI BLAST+ (`blastn` and `makeblastdb`).

```bash
git clone https://github.com/happykhan/matryoshka.git
cd matryoshka

# macOS
brew install blast
python3 -m venv .venv
source .venv/bin/activate
python -m pip install .

# Ubuntu/Debian: install BLAST first with:
# sudo apt-get install ncbi-blast+
```

Run the bundled Tn1/Tn2/Tn3 acceptance example:

```bash
matryoshka run tests/test-data/partridge-examples/Tn1-Tn2-Tn3.fasta \
  --out results
```

The result directory contains:

- `annotation.json` — versioned hierarchy with stable feature IDs and provenance;
- `annotation.gff3` — one valid multi-record GFF3 document;
- `run.json` — short run summary;
- `mara/*.svg` — locus-based MARA-style maps;
- `mara-table/*.svg` — MARA-style annotation tables.

Each Tn1/Tn2/Tn3 target receives its own viewport, so the element stays legible in a
long plasmid or chromosome. Near matches are labelled `Tn1-like`, `Tn2-like` or
`Tn3-like`; supported insertions, deletions, fragments and reverse-complement matches
retain their evidence in the JSON and table.

For these three elements, the diagram is assembled from sequence evidence rather than
drawn from a name alone. Matryoshka independently scans for both terminal IRs,
`blaTEM`, `tnpR`, `res` and `tnpA`; checks their order and orientation; and then uses
the closest whole-locus reference to assign Tn1, Tn2, Tn3 or a qualified `*-like`
name. Missing required components prevent an exact named call. Solid arrows in the
figure are sequence-detected components; outlined dashed arrows are reference
projections and are used only where a supported parent lacks a component call.

### Pixi alternative

[Pixi](https://pixi.sh) provides a locked core environment for Linux and Intel/Apple
Silicon macOS:

```bash
pixi install
pixi run smoke
```

ISEScan, AMRFinder+ and IntegronFinder environments are currently Linux-only. The core
validated-reference workflow works on all listed platforms.

## Run your own sequence

```bash
matryoshka run assembly.fasta --out results --threads 4
```

The default `--profile validated` includes the Sally-backed Tn1/Tn2/Tn3, curated unit
transposon, integron and ISEcp1-TPU references. `--profile all` enables broader legacy
and experimental references and should be treated as exploratory.

BLAST-based detection is always run. Extra detector evidence can be supplied without
rerunning tools:

```bash
matryoshka run assembly.fasta --out results \
  --isescan isescan.tsv \
  --amrfinder amrfinder.tsv \
  --integrons sample.integrons
```

On Linux, installed detectors can be invoked automatically:

```bash
# Run any detector available on PATH (or through this checkout's Pixi environments).
matryoshka run assembly.fasta --out results --detectors available

# Require all three; fail clearly if one is unavailable.
matryoshka run assembly.fasta --out results --detectors all
```

The lower-level `matryoshka annotate` command emits one selected format and remains
useful for scripts. `matryoshka run` is the supported reproducible workflow for new
users.

## What it currently identifies

| Category | Implemented support |
|---|---|
| Tn1/Tn2/Tn3 | Independent IR/blaTEM/tnpR/res/tnpA detection, orientation-aware component assembly, canonical naming, closest-reference minor variants, indels, reverse orientation and conservative partial/ambiguous calls |
| Curated MARA units | Tn21, Tn1721, Tn1722, Tn4401, Tn5393 and Tn5403 with component-aware maps |
| ISEcp1 transposition units | Seven complete supplied ISEcp1–blaCMY references and orientation-aware incomplete candidates |
| Insertion sequences | Exact curated IS calls plus ISEScan calls, including IS26/IS257, ISEcp1, ISApl1 and IS91/ISCR families |
| Composite/signature elements | Tn4401, Tn1999, Tn6330, Tn2006, Tn125, Tn10, Tn5, Tn4001, Tn1546 and Tn1331 rules |
| Integrons | IntegronFinder parsing, cassette reconstruction, class-1-integron and complete Tn402 component inference |
| Other references | Tn4401 variants, Acinetobacter islands, GIsul2, Tn7, Tn552, Tn1546, Tn1331, replicons and documented motifs in the `all` profile |

The source-backed [MARA component inventory](docs/mara-component-inventory.md) defines
30 raw component classes and 18 compound-element grammars. It distinguishes what is
implemented, partially implemented and still missing; the catalogue is also available
to software:

```bash
matryoshka catalog --format tsv
matryoshka catalog --format json --out mara-catalog.json
```

Biological and representational targets are defined in the
[MARA parity specification](docs/mara-parity-spec.md). Unfinished capabilities remain
explicit in [GAPS.md](GAPS.md); this is not yet a general replacement for expert MARA
annotation.

## Supported outputs

The result-directory contract is documented in
[docs/output-contract.md](docs/output-contract.md), and its JSON Schema is
[annotation-v1.schema.json](docs/schema/annotation-v1.schema.json).

The lower-level `annotate --format` choices are:

| Format | Contents |
|---|---|
| `json` | Legacy single-output hierarchy |
| `gff3` | Feature rows with `Parent` relationships |
| `genbank` | Biopython GenBank feature annotations |
| `wolvercote` | Compact nested cell-format text |
| `linear` | Scale-accurate whole-record SVG |
| `mara` | One feature-specific locus SVG per validated target, with an on-figure symbol key |
| `mara-table` | One MARA-style annotation-table SVG per locus, with evidence notes and a compact key |

Legacy circular CellGen SVG/PNG rendering is not distributed in this alpha because it
required an unpublished sibling package. Use the bundled linear and MARA SVG outputs.

## Reference and validation data

The machine-readable [reference manifest](matryoshka/references/manifest.yaml) pins the
database version, file hashes, record IDs, accessions, provenance, profile and scan
thresholds for every bundled FASTA. The human-readable companion is
[matryoshka/references/MANIFEST.md](matryoshka/references/MANIFEST.md).

Three complete public NCBI plasmids are committed as acceptance fixtures:

| Plasmid | Accession | Acceptance purpose |
|---|---|---|
| pKpQIL | NC_014016.1 | Tn4401 reference detection |
| pOXA-48a | JN626286.1 | Tn1999/TSD and negative Tn4401 control |
| pEK499 | EU935739.1 | ISEcp1 and false-positive controls |

Their SHA-256 hashes and refresh script are in `tests/test-data/acceptance/`. Compact
detector parser fixtures are clearly labelled as hand-authored. Sally's private working
documents are not redistributed; see [docs/redistribution.md](docs/redistribution.md).

## Development and verification

```bash
python -m pip install -e ".[dev,viz]"
pytest -q
ruff check matryoshka tests scripts
python -m build
```

CI tests Python 3.11, 3.12 and 3.13 with BLAST installed, validates the checksummed
reference manifest, builds the wheel, installs it in a clean environment and runs the
complete result-directory workflow.

## Version and status

`0.1.0a1` is an external alpha. The validated profile is suitable for demonstration
and repeatable evaluation of the documented subset, not unattended clinical or
surveillance decisions. Calls retain coordinates, reference evidence and confidence so
an expert can review them.

Changes are recorded in [CHANGELOG.md](CHANGELOG.md). Please report reproducible issues
through [GitHub Issues](https://github.com/happykhan/matryoshka/issues).

## Licence and citation

Code is MIT licensed. Reference-data provenance and individual licences are recorded in
the manifest. Citation metadata is in [CITATION.cff](CITATION.cff).
