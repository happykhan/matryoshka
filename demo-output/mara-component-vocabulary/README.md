# MARA component-vocabulary demo

This acceptance artifact demonstrates the raw-feature symbols independently of
sequence-detection coverage. It intentionally combines an integron, an ISEcp1
TPU, ISCR sites, plasmid context and incomplete evidence so every symbol can be
reviewed.

- `mara-component-vocabulary.svg`: MARA-compatible diagram.
- `mara-component-table.svg`: companion Position/Name/FID/Type/Notes table.
- `mara-component-catalog.tsv`: raw component list.
- `mara-component-catalog.json`: complete component and assembly grammar.

Regenerate from the repository root:

```bash
PYTHONPATH=. python scripts/make_mara_component_demo.py demo-output/mara-component-vocabulary
```
