# Arbitrary-sequence Tn1/Tn2/Tn3 demonstration

`arbitrary-demo.fasta` is a deterministic 44,647 bp synthetic contig generated
by `scripts/make_tn123_demo.py`. It contains:

- a Tn1 reference with 25 substitutions (called `Tn1-like`, 99.5% identity);
- a Tn2 reference with an 800 bp internal insertion (called `Tn2-like` and
  assembled as one interrupted locus); and
- an exact reverse-complement Tn3 (called `Tn3`, reverse strand).

`calls.json` is the machine-readable hierarchy. `mara/` contains one readable
diagram per detected locus and `mara-table/` contains the corresponding
annotation tables using absolute positions on the original contig. Every SVG has
an on-figure key. In all three loci, both terminal IRs, `blaTEM`, `tnpR`, `res`
and `tnpA` are independently detected from sequence and then assembled in the
required orientation. Solid arrows show those detected components; an outlined
dashed arrow would identify a curated reference projection.

`proof-bundle/` is the one-command acceptance result. Open
`proof-bundle/proof/report.html` to review the full chain for each locus. The
same directory contains the raw component ledger, known-element match ledger,
machine-checkable PASS/PARTIAL/FAIL proof, original hierarchical SVG, MARA map and table,
and `annotation.cell` CellGen/Wolvercote nested cell-format output.

The match ledger also names the versioned YAML definition used for each call and
records the type, subtype, substitutions, inserted bases and deleted bases. The
current regenerated bundle gives three `PASS` verdicts: the two qualified variants
retain their differences, while only the unchanged Tn3 receives an exact name.

Reproduce it with:

```bash
matryoshka run demo-output/arbitrary-tn123/arbitrary-demo.fasta \
  --out demo-output/arbitrary-tn123/proof-bundle
```
