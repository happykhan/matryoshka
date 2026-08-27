# Arbitrary-sequence Tn1/Tn2/Tn3 demonstration

`arbitrary-demo.fasta` is a deterministic 44,647 bp synthetic contig generated
by `scripts/make_tn123_demo.py`. It contains:

- a Tn1 reference with 25 substitutions (called `Tn1-like`, 99.5% identity);
- a Tn2 reference with an 800 bp internal insertion (called `Tn2-like` and
  assembled as one interrupted locus); and
- a reverse-complement Tn3-derived sequence (called `Tn3-like`, reverse strand).

`calls.json` is the machine-readable hierarchy. `locus-map/` contains one readable
diagram per detected locus and `locus-table/` contains the corresponding
annotation tables using absolute positions on the original contig. Every SVG has
an on-figure key. In all three loci, both terminal IRs, `blaTEM`, `tnpR`, `res`
and `tnpA` are independently detected from sequence and then assembled in the
required orientation. Solid arrows show those detected components; an outlined
dashed arrow would identify a curated reference projection.

`proof-bundle/` is the one-command acceptance result. Open
`proof-bundle/proof/report.html` to review the full chain for each locus. The
same directory contains the raw component ledger, known-element match ledger,
machine-checkable PASS/PARTIAL/FAIL proof, original hierarchical SVG, locus map and table,
and `annotation.cell` CellGen/Wolvercote nested cell-format output.

The proof ledger records the expert-rule type, categorical `tnpR`/`tnpA` backbone
haplotype, component structure and optional reference context. This bundle deliberately uses the
`tn123-components` profile, so it contains no complete Tn1/Tn2/Tn3 lookup. All three
loci nevertheless receive `PASS`; the interrupted Tn2-like locus retains and draws
its approximately 799 bp insertion from split `tnpA` component evidence.

Reproduce it with:

```bash
matryoshka run demo-output/arbitrary-tn123/arbitrary-demo.fasta \
  --out demo-output/arbitrary-tn123/proof-bundle \
  --profile tn123-components --detectors none
```
