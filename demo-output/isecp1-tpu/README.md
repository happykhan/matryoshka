# ISEcp1 TPU demonstration

These outputs were generated from the seven complete CMY-like transposition-unit
sequences in `tests/test-data/reviewed-examples/TPU_CMY-2-like.fasta`.

- `locus-map/` contains one locus diagram for each complete reference match.
- `locus-table/` contains the matching position/FID/type/evidence table.
- The white arrow is the exact 1,656 bp ISEcp1 reference and points toward IRR.
- Hollow lollipops show the expected 5 bp outer direct repeat where the supplied
  element-only sequence does not include enough flanking DNA to confirm it.
- Every map and table carries a key explaining the arrows, flags and filled or
  hollow boundary lollipops.

These files validate complete-unit identity and ISEcp1 orientation. They do not
yet claim full internal gene/IRalt component reconstruction; that remaining work
is tracked in `docs/locus-annotation-spec.md`.
