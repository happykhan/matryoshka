# Tn1 real-accession component-rule evaluation

Run date: 2026-08-27

This evaluation tests whether the component grammar can recover Tn1-like loci
without searching the input against a complete Tn1 sequence. Sixteen NCBI
records returned by a Tn1/transposon/plasmid metadata search were downloaded as
unannotated FASTA and processed with `--profile component-rules --detectors none`.
The NCBI annotations were used to select the panel, not as input to Matryoshka.

```bash
python scripts/evaluate_tn1_accessions.py --out results/tn1-accessions
```

Six records produced a complete component-rule call. PP591960.1, PP591959.1,
AM261760.1, KX156773.1 and HM804085.1 yielded compact Tn1-like loci. GQ160960.2
yielded a 10,491 bp Tn1-like locus containing the known large internal mer-region
insertion. The latter exposed and fixed a former fixed-window limitation: candidate
boundaries are now chosen from compatible terminal-IR pairs and the declared
maximum span, rather than an arbitrary 2.5 kb neighbourhood around `tnpA`.

The remaining ten records did not receive a complete Tn1/Tn2/Tn3 name. Nine retain
one or more independently detected components in the annotations and figures, while
AY150843.2 has no supported group component profile. This is the intended conservative
result for partial, interrupted, dispersed or misleading metadata hits.

The complete ledger is in
[`tn1-accession-panel-results.tsv`](tn1-accession-panel-results.tsv). Each normal run
also emits JSON, GFF3, GenBank, cell format, hierarchy SVG, locus maps and locus tables.
