# Contributing to Matryoshka

Thanks for your interest — contributions are very welcome, particularly:

- New transposon / integron / genomic-island inference rules
- New reference sequences for the BLAST library
- Improved boundary detection (TSD, IR, res sites)
- Additional output formats
- Documentation and examples

## Development setup

```bash
git clone https://github.com/happykhan/matryoshka
cd matryoshka
pixi install
pixi run pytest            # run tests
pixi run -e dev ruff check . # lint (once dev env is set up)
```

For pip-only development without pixi:

```bash
pip install -e ".[dev,viz]"
pytest
```

## Adding a new transposon rule

Most transposon additions are a **one-line change** to a rule table in
[`matryoshka/transposon.py`](matryoshka/transposon.py). No new code needed.

### Flanked-cargo composite (two flanking IS + cargo)

```python
FLANKED_RULES.append(FlankedRule(
    family="TnXXXX", name="TnXXXX",
    cargo_match="KEY_GENE",             # substring match in AMR name
    left_family="ISxx", right_family="ISyy",
    upstream_max=500, downstream_max=1000,
))
```

### One-ended capture (single IS + cargo, like ISEcp1)

```python
ONE_ENDED_RULES.append(OneEndedRule(
    family="ISxxx_capture", name="ISxxx_capture",
    is_family="ISfamily",
    cargo_matches=("gene1", "gene2"),   # tuple of substrings
    max_gap=1500,
))
```

### Signature-only (cargo-gene cluster with no flanking IS)

```python
SIGNATURE_RULES.append(SignatureRule(
    family="TnXXXX", name="TnXXXX",
    required_genes=("geneA", "geneB", "geneC"),
    max_span=10_000,
))
```

Add a short test in `tests/test_transposon_rules.py` that constructs a
synthetic `MGEFeature` list and asserts the rule fires.

## Adding a BLAST reference

1. Add an entry in [`scripts/fetch_mge_references.py`](scripts/fetch_mge_references.py)
   (or drop a FASTA directly into `matryoshka/references/`).
2. Ensure the FASTA header encodes metadata: `>name accession element_type=X family=Y …`
3. Register appropriate thresholds in `REFERENCE_PARAMS` in
   [`matryoshka/reference_scan.py`](matryoshka/reference_scan.py).
4. Update `matryoshka/references/MANIFEST.md` with the source and SHA-256.

## Code style

- Python 3.11+
- `from __future__ import annotations` at module top
- Type hints on all public functions
- Docstrings explain *why*, not *what*
- `ruff check .` must be clean before opening a PR
- `pytest` must be green

## Tests

Run the full suite:

```bash
pixi run pytest tests/ -q
```

BLAST-dependent tests skip automatically when `blastn` is not on PATH. Tests
that need the bundled reference FASTAs skip automatically if those aren't
present.

## Opening a pull request

1. Fork + branch from `main`.
2. Include tests for new behaviour.
3. Update `CHANGELOG.md` under `[Unreleased]`.
4. Describe the biological rationale (paper reference, accession, structural
   feature) in the PR body — this is the hardest part for reviewers to verify
   from code alone.

## Code of Conduct

Please be kind. Follow the
[Contributor Covenant 2.1](https://www.contributor-covenant.org/version/2/1/code_of_conduct/).
