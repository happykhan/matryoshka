# Matryoshka

Component-first mobile genetic element annotation and locus visualisation for
bacterial assemblies.

Matryoshka scans an arbitrary FASTA sequence for constituent mobile-element features,
assembles compatible features into loci, applies reviewable expert rules, and draws
the detected evidence. For the validated Tn1/Tn2/Tn3 workflow, a whole-element
reference hit is not required to discover a locus or assign a qualified type-like
call. Complete reviewed sequences are used separately to confirm exact names and
subtypes.

## How Tn1/Tn2/Tn3 classification works

The primary classification path is component-first:

| Stage | What the program does | Role of reviewed sequences |
|---|---|---|
| Detect components | Search independently for both terminal IRs, `blaTEM`, `tnpR`, `res` and `tnpA` | Reviewed examples supply local component profiles, not a pre-drawn locus |
| Assemble a locus | Test component count, order, orientation, spacing and terminal boundaries on either strand | No complete-element match is required |
| Apply expert rules | Compare the detected component profiles with the YAML-defined Tn1/Tn2/Tn3 signatures and require the configured score and margin | The component evidence assigns `Tn1-like`, `Tn2-like` or `Tn3-like` |
| Describe variation | Retain substitutions, insertions, deletions, interruptions, missing ends and ambiguous evidence | Differences remain visible rather than being replaced by a canonical structure |
| Confirm a reviewed name | Optionally compare the assembled complete locus with reviewed whole-element sequences | This can upgrade a compatible call to an exact type/subtype or add closest-reference context; it cannot create or override the component-rule call |
| Render outputs | Draw the hierarchy, locus maps and locus table from the detected/assembled feature record | Projected context is visually distinguished from independently detected components |

The rules and thresholds are maintained in human-readable definition files, not
embedded as family-specific conditionals in the classifier. Tn1/Tn2/Tn3 live in
[`matryoshka/tn123_definitions.yaml`](matryoshka/tn123_definitions.yaml); Tn21,
Tn1721 and Tn1722 live in
[`matryoshka/unit_definitions.yaml`](matryoshka/unit_definitions.yaml).

This is still homology-based component detection. The current profiles are intended
to recover close and moderately varied versions of recognised Tn1/Tn2/Tn3 components,
including split alignments and indels. Very remote homologues or entirely new component
classes will require additional nucleotide, protein or profile-HMM detectors; the tool
does not claim unrestricted de novo mobile-element discovery.

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
matryoshka run tests/test-data/reviewed-examples/Tn1-Tn2-Tn3.fasta \
  --out results
```

To prove that the type-like call does not depend on complete-element lookup, run the
44,647 bp arbitrary-contig demonstration with the component-only profile:

```bash
matryoshka run demo-output/arbitrary-tn123/arbitrary-demo.fasta \
  --out results/arbitrary-component-only \
  --profile tn123-components --detectors none
```

That checked-in input produces three independently assembled loci:

| Result | Evidence recovered from the input | Complete-element lookup | Verdict |
|---|---|---|---|
| `Tn1-like` | Six required components; 25 substitutions across the Tn1-derived locus | Not run | `PASS` |
| `Tn2-like` | Six required components; split `tnpA` alignment retaining and drawing a 799 bp insertion | Not run | `PASS` |
| `Tn3-like` | Six required components in reverse orientation | Not run | `PASS` |

The generated proof page reports the component coordinates, identities, coverage,
grammar result, component-profile scores, structural differences and links to every
corresponding output. The committed example is
[`demo-output/arbitrary-tn123/proof-bundle/proof/report.html`](demo-output/arbitrary-tn123/proof-bundle/proof/report.html).

The result directory contains:

- `annotation.json` — versioned hierarchy with stable feature IDs and provenance;
- `annotation.gff3` — one valid multi-record GFF3 document;
- `annotation.gbk` — the same hierarchy as GenBank feature annotations;
- `annotation.cell` — nested CellGen/Wolvercote cell-format text;
- `run.json` — short run summary;
- `hierarchy/*.svg` — the original scale-accurate nested hierarchy view;
- `locus-map/*.svg` — one readable locus map per supported target;
- `locus-table/*.svg` — the corresponding locus annotation tables;
- `proof/report.html` — a human-readable component-to-call proof report;
- `proof/proof.json`, `components.tsv` and `matches.tsv` — machine-checkable
  component, grammar and known-element evidence.

Each Tn1/Tn2/Tn3 target receives its own viewport, so the element stays legible in a
long plasmid or chromosome. Near matches are labelled `Tn1-like`, `Tn2-like` or
`Tn3-like`; supported insertions, deletions, fragments and reverse-complement matches
retain their evidence in the JSON and table. Solid arrows in the figure are
sequence-detected components; outlined dashed arrows are reference projections and
are used only where a supported parent lacks a component call.

For a Tn1/Tn2/Tn3 locus to receive `PASS` in the proof report, all required
components must be independently sequence-detected, their grammar and orientation
must be valid, the component-rule classifier must assign a type, and no internal
component may be supplied only by reference projection. A whole-element match is not
required for `PASS`. Missing required components prevent a complete named call.
Expected fragments receive a `PARTIAL` locus verdict and a
`PARTIAL_TN123_EVIDENCE` run status rather than being misrepresented as either a
complete pass or a pipeline failure.

The biological naming rules are data, not embedded conditionals. They are maintained
in [`matryoshka/tn123_definitions.yaml`](matryoshka/tn123_definitions.yaml), validated
by the definition loader, documented by a [JSON Schema](docs/schema/tn123-definitions-v1.schema.json),
and exportable in review-friendly or machine-readable form:

```bash
matryoshka definitions --format markdown --out tn123-definitions-review.md
matryoshka definitions --format yaml --out tn123-definitions.yaml
matryoshka definitions --format json --out tn123-definitions.json

# Review the Tn21/Tn1721/Tn1722 component grammars and feature profiles
matryoshka definitions --set reviewed-units --format markdown \
  --out reviewed-unit-definitions.md
```

The current definitions include three expert-reviewed canonical elements plus
Tn2c (HM749967), Tn2.1 (CP028717), Tn1Mer (GQ160960) and the legacy V00613 Tn3
sequence. See the [definition and validation report](docs/tn123-definition-and-validation-report.md)
and its [real-accession ledger](docs/validation/tn123-real-accession-results.tsv).
The [reviewed-definition demonstration](demo-output/tn123-reviewed-definitions/README.md)
includes the complete result bundle and every generated locus map/table.
The self-contained Word/PDF evidence report in `reports/` embeds those figures,
the original hierarchy views, component ledgers, the natural pEK499 partial
fragments and the arbitrary-contig demonstration. The complete pEK499 bundle is
in `demo-output/tn123-pEK499/proof-bundle/`.

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

The default `--profile validated` includes the expert-reviewed Tn1/Tn2/Tn3, curated unit
transposon, integron and ISEcp1-TPU references. `--profile all` enables broader legacy
and experimental references and should be treated as exploratory.

To exercise and audit discovery without any complete-element lookup, use the
component-rule profile:

```bash
matryoshka run assembly.fasta --out results --profile component-rules \
  --detectors none
```

In that mode, a compatible candidate can still receive a `Tn1-like`, `Tn2-like`,
`Tn3-like`, `Tn21-like`, `Tn1721-like` or `Tn1722-like` call from its detected
component composition, order, orientation, spacing and component-profile scores.
Exact reviewed subtype names require the optional secondary whole-locus confirmation
available in the `validated` profile.

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
| Tn1/Tn2/Tn3 | YAML-defined IR/blaTEM/tnpR/res/tnpA grammar, orientation-aware component-profile classification, optional exact reference confirmation, indels, reverse orientation and conservative partial/ambiguous calls |
| Reviewed unit transposons | Component-rule discovery for Tn21, Tn1721 and Tn1722; reference-supported maps for Tn4401, Tn5393 and Tn5403 |
| ISEcp1 transposition units | Seven complete supplied ISEcp1–blaCMY references and orientation-aware incomplete candidates |
| Insertion sequences | Exact curated IS calls plus ISEScan calls, including IS26/IS257, ISEcp1, ISApl1 and IS91/ISCR families |
| Composite/signature elements | Tn4401, Tn1999, Tn6330, Tn2006, Tn125, Tn10, Tn5, Tn4001, Tn1546 and Tn1331 rules |
| Integrons | IntegronFinder parsing, cassette reconstruction, class-1-integron and complete Tn402 component inference |
| Other references | Tn4401 variants, Acinetobacter islands, GIsul2, Tn7, Tn552, Tn1546, Tn1331, replicons and documented motifs in the `all` profile |

The source-backed [locus component inventory](docs/locus-component-inventory.md) defines
30 raw component classes and 18 compound-element grammars. It distinguishes what is
implemented, partially implemented and still missing; the catalogue is also available
to software:

```bash
matryoshka catalog --format tsv
matryoshka catalog --format json --out locus-catalog.json
```

Biological and representational targets are defined in the
[locus annotation specification](docs/locus-annotation-spec.md). Unfinished capabilities remain
explicit in [GAPS.md](GAPS.md); this is not yet a general replacement for expert
mobile-element annotation.

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
| `locus-map` | One feature-specific locus SVG per validated target, with an on-figure symbol key |
| `locus-table` | One locus annotation-table SVG per locus, with evidence notes and a compact key |

Legacy circular CellGen SVG/PNG rendering is not distributed in this alpha because it
required an unpublished sibling package. Use the bundled linear, locus-map and
locus-table SVG outputs.

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
detector parser fixtures are clearly labelled as hand-authored. Private expert-review
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
