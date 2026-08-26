# Tn1/Tn2/Tn3 MARA demonstration

This demonstration uses Sally Partridge's `Tn1-Tn2-Tn3.fasta` without ISEScan,
AMRFinder+ or IntegronFinder placeholder files. Matryoshka's bundled reference scan
detects the three complete canonical elements and projects their curated internal
maps into the new MARA-style SVG output.

| Input | Call | Reference | Identity | Covered reference |
| --- | --- | --- | ---: | ---: |
| `Tn1_NC_008357` | Tn1 | NC_008357 | 100% | 100% |
| `Tn2_AY123253` | Tn2 | AY123253 | 100% | 100% |
| `Tn3_HM749966` | Tn3 | HM749966 | 100% | 100% |

The diagrams now follow MARA's compact Figure 1/2 visual grammar: solid,
black-outlined transposon blocks, black directional gene arrows, terminal inverted
repeats as small flags, bounded `res` sites, jagged fragment edges, red dashed
unannotated sequence and paired boundary lollipops. Hollow lollipops indicate the
expected 5-bp TSD for these elements when flanking sequence is unavailable; filled,
labelled pairs are reserved for sequence-confirmed TSDs.

Short homologous segments below 98% reference coverage are reported as
`Tn1/2/3` fragments, with a jagged outline, rather than being forced into an exact
Tn1, Tn2 or Tn3 call.

Reproduce the output with:

```bash
python -m matryoshka annotate tests/test-data/partridge-examples/Tn1-Tn2-Tn3.fasta \
  --format mara --out demo-output/tn123-mara
```

The pre-existing nested visualisation is unchanged and remains available with
`--format linear`; comparison files are in `../tn123-linear/`.

MARA-style annotation tables can be generated from the same hierarchy with
`--format mara-table`; examples are in `../tn123-mara-table/`.
