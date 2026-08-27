# pEK499 Tn1/Tn2/Tn3 partial-fragment output

This is the complete `matryoshka run` output for the natural pEK499 sequence
(accession EU935739.1). The validated built-in reference and component scans
retain two separate Tn2-related regions at 38,747–40,671 and 60,316–62,561.
Both are reported as `Tn1/2/3 fragment`; neither is promoted to complete Tn2.

The two fragment calls have locus maps and locus tables in the proof
bundle. Other validated-profile mobile-element loci found in the same record
are also drawn, demonstrating that the output is locus based rather than
restricted to a single expected element.

Reproduce the bundle with:

```bash
matryoshka run tests/test-data/acceptance/pEK499.fasta \
  --out demo-output/tn123-pEK499/proof-bundle \
  --detectors none
```
