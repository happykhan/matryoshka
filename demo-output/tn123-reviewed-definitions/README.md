# Reviewed Tn1/Tn2/Tn3 definition output bundle

This directory is the complete `matryoshka run` output for the seven reviewed
whole-element definitions in `matryoshka/references/tn1_tn2_tn3.fasta`:

- canonical Tn1 (NC_008357);
- canonical Tn2 (AY123253);
- canonical Tn3 (HM749966);
- Tn2c (HM749967);
- Tn2.1 (CP028717);
- Tn1Mer (GQ160960); and
- the declared V00613 legacy Tn3 subtype.

Open `proof-bundle/proof/report.html` for the linked evidence report. The
`proof-bundle/mara/` directory contains a MARA locus SVG for each definition;
`proof-bundle/mara-table/` contains the corresponding tables. Tn2.1 also has a
separate locus entry for its nested ISEcp1-associated insertion, so the bundle
contains eight locus maps in total.

`proof-bundle/run.json` has a `mara_locus_outputs` index giving the record,
call, coordinates and exact relative path of every MARA map and table. All seven
Tn1/Tn2/Tn3 whole-element calls receive a proof `PASS`.

Reproduce the bundle with:

```bash
matryoshka run matryoshka/references/tn1_tn2_tn3.fasta \
  --out demo-output/tn123-reviewed-definitions/proof-bundle \
  --detectors none
```
