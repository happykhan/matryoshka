# Priority transposon implementation roadmap

The priority scope is 17 named elements. The governing machine-readable plan is
[`matryoshka/priority_transposons.yaml`](../matryoshka/priority_transposons.yaml);
`matryoshka roadmap` exports it for review.

## Current position

| Evidence level | Targets | Meaning |
| --- | --- | --- |
| Validated component rule | Tn1, Tn2, Tn3, Tn21, Tn1721, Tn1722 | Independent components assemble a complete family-like locus without whole-element lookup. |
| Implemented, broader validation required | Tn402, Tn4401, Tn10, Tn1331, Tn1999 | Some component or signature logic exists, but boundaries, specificity or variant panels are incomplete. |
| Reference-supported | Tn5393, Tn5403, Tn7, Tn1696 | The locus can be recognised from curated reference evidence; this is not yet general component-driven discovery. |
| Definition required | Tn6029, Tn2670 | The reviewed component grammar and validation data must be encoded before a named rule is enabled. |

## Delivery sequence

1. Complete Tn402, Tn4401, Tn5393 and Tn5403 definitions from the supplied component
   sequences. Replace broad cargo shortcuts with explicit ends, transposition genes,
   resolution sites, cargo and permitted interruptions.
2. Encode Tn1696 as a Tn21-subgroup rule with distinguishing integron/mer structure;
   add Tn7 protein/end/target-site profiles; tighten Tn10, Tn1331 and Tn1999 rules.
3. Reconstruct Tn6029 from reviewed evidence without using IS26 proximity as a name.
4. Encode Tn2670 last as a nested validation case. Its purpose is to prove that the
   same independently detected Tn21, integron, Tn402, cassette and outer composite
   components can be assembled into overlapping/nested parents and rendered without
   flattening the evidence.

## Definition of done for each target

- The YAML describes required/optional components, counts, order, orientation, span,
  ends, TSD expectation, permitted variation, exclusions and naming thresholds.
- Atomic profiles are independently searchable; a complete-element reference is not
  required for a family-like call.
- Exact naming requires a complete compatible grammar plus secondary reviewed-reference
  evidence.
- Canonical, reverse, mutated, inserted/deleted, partial, long-contig and near-neighbour
  controls pass.
- JSON, GFF3, GenBank, CellGen, hierarchy SVG, locus map, locus table and evidence ledger
  all describe the same detected components and classification.
