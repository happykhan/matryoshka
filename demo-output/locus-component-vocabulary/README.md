# Locus component-vocabulary demo

This acceptance artifact demonstrates the raw-feature symbols independently of
sequence-detection coverage. It intentionally combines an integron, an ISEcp1
TPU, ISCR sites, plasmid context and incomplete evidence so every symbol can be
reviewed.

- `locus-component-vocabulary.svg`: locus-map diagram.
- `locus-component-table.svg`: companion Position/Name/FID/Type/Notes locus table.
- `locus-component-catalog.tsv`: raw component list.
- `locus-component-catalog.json`: complete component and assembly grammar.

Regenerate from the repository root:

```bash
pixi run python scripts/make_locus_component_demo.py demo-output/locus-component-vocabulary
```
