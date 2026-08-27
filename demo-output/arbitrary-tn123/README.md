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
an on-figure key; the table Notes distinguish detected parent evidence from
curated reference-projected internal components.
