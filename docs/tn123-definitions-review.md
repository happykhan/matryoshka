# Tn1/Tn2/Tn3 expert definitions

Definition version: `2026-08-27`

Reviewable biological rules for detecting the Tn1/Tn2/Tn3 group, assigning a type or named subtype, describing variants, and retaining incomplete fragments without over-naming them.

## Family rule

A unit transposon bounded by two 38 bp terminal inverted repeats and containing blaTEM, tnpR, res and tnpA in the conserved element-relative order. This grammar establishes a Tn1/Tn2/Tn3-group unit. It does not by itself distinguish Tn1, Tn2 and Tn3 because their structures and most of their nucleotide sequences are shared.

Required component order:

`terminal_IR -> blaTEM -> tnpR -> res -> tnpA -> terminal_IR`

A complete group-level call requires all six component occurrences in the stated order. Reverse-complement copies use the reversed order. Missing components, conflicting order or a missing end produce a fragment or an unresolved Tn3-family unit rather than an exact Tn1, Tn2 or Tn3 call.

## Classification rules

Type assignment: First establish a complete Tn1/Tn2/Tn3-group locus from independently detected components. Score Tn1, Tn2 and Tn3 using the locally detected blaTEM, tnpR, res and tnpA component profiles and the declared weights. Assign a rule-based type only when the best component profile clears the score and margin thresholds. The score threshold permits a locally interrupted or deleted component while the required-component and order rules remain satisfied. Otherwise retain a complete but unresolved Tn1/Tn2/Tn3-group element. A whole-locus reference match is not required to discover or type a candidate.

Secondary reference comparison: Compare an assembled candidate with reviewed whole-element references only after component-rule classification. Record the closest reference, identity, coverage and differences. An exact match may confirm a known definition, but a closest match must not create the family or type call.

Exact definition match: Exact means identical to a declared reference definition. A sequence containing even one substitution or indel is a variant unless that sequence has its own reviewed subtype definition.

Close variant: A complete candidate that is clearly closest to one type but is not identical to a declared definition is reported as type-like, with its substitutions, insertions and deletions retained in the evidence.

Fragment: Partial homology is retained with coordinates, detected components and closest reference. A fragment is not promoted to a complete named type.

## Type: Tn1

Canonical Tn1 uses the NC_008357 sequence selected by Sally. It has the complete Tn1/Tn2/Tn3 grammar and carries blaTEM-2. Its local blaTEM, tnpR, res and tnpA profiles define the Tn1 component signature. Exact Tn1 additionally requires a gap-free, mismatch-free match to this declared sequence.

| Order | Role | Name | Coordinates | Strand |
| ---: | --- | --- | --- | :---: |
| 1 | terminal_IR | IRL | 1..38 | . |
| 2 | blaTEM | blaTEM-2 | 148..1008 | - |
| 3 | tnpR | tnpR | 1191..1748 | - |
| 4 | res | res | 1754..1867 | . |
| 5 | tnpA | tnpA | 1911..-34 | + |
| 6 | terminal_IR | IRR | -38..-1 | . |

## Type: Tn2

Canonical Tn2 uses Sally's AY123253 sequence, described as Tn2* in the 2005 paper and selected as Tn2 in her working document. It has the complete group grammar and carries blaTEM-1b. Its local component profiles, particularly the diagnostic sequence around res, distinguish a Tn2-like candidate from Tn1 and Tn3 without requiring a whole-element match.

| Order | Role | Name | Coordinates | Strand |
| ---: | --- | --- | --- | :---: |
| 1 | terminal_IR | IRL | 1..38 | . |
| 2 | blaTEM | blaTEM-1b | 148..1008 | - |
| 3 | tnpR | tnpR | 1191..1748 | - |
| 4 | res | res | 1754..1867 | . |
| 5 | tnpA | tnpA | 1911..-34 | + |
| 6 | terminal_IR | IRR | -38..-1 | . |

## Type: Tn3

Canonical Tn3 uses HM749966, the sequence Sally selected in preference to the older V00613 record. It has the complete group grammar and carries blaTEM-1a. The combined local component profiles, not blaTEM protein identity or gene order alone, distinguish a Tn3-like candidate from Tn1 and Tn2.

| Order | Role | Name | Coordinates | Strand |
| ---: | --- | --- | --- | :---: |
| 1 | terminal_IR | IRL | 1..38 | . |
| 2 | blaTEM | blaTEM-1a | 148..1008 | - |
| 3 | tnpR | tnpR | 1191..1748 | - |
| 4 | res | res | 1754..1867 | . |
| 5 | tnpA | tnpA | 1911..-34 | + |
| 6 | terminal_IR | IRR | -38..-1 | . |

## Named definitions and subtypes

### Tn1 (`Tn1_NC_008357`)

Type: `Tn1`. Subtype: `canonical`. Source accession: `NC_008357`.

Exact, complete, gap-free match to NC_008357 plus complete component grammar.

Review status: Sally-selected canonical sequence.

### Tn2 (`Tn2_AY123253`)

Type: `Tn2`. Subtype: `canonical`. Source accession: `AY123253`.

Exact, complete, gap-free match to AY123253 plus complete component grammar.

Review status: Sally-selected canonical sequence.

### Tn3 (`Tn3_HM749966`)

Type: `Tn3`. Subtype: `canonical`. Source accession: `HM749966`.

Exact, complete, gap-free match to HM749966 plus complete component grammar.

Review status: Sally-selected canonical sequence.

### Tn2c (`Tn2c_HM749967`)

Type: `Tn2`. Subtype: `Tn2c`. Source accession: `HM749967`.

Exact, complete, gap-free match to the 4950 bp Tn2 interval in HM749967. It retains the Tn2 grammar but carries the blaTEM-1c sequence and nine substitutions relative to AY123253.

Review status: Sally called this Tn2c in MARA and marked it for later addition.

Known differences from the parent definition:

- 9 nucleotide substitutions relative to canonical Tn2 AY123253
- blaTEM-1c allele

### Tn2.1 (`Tn2_1_CP028717`)

Type: `Tn2`. Subtype: `Tn2.1`. Source accession: `CP028717`.

Exact match to the 8979 bp interrupted Tn2 structure extracted from CP028717. The Tn2 backbone is complete, including both terminal IRs, and a 4024 bp ISEcp1-associated region interrupts the sequence between canonical reference positions 124 and 129.

Review status: Sally described this as Tn2 with an inserted ISEcp1 TPU.

Known differences from the parent definition:

- 4024 bp ISEcp1-associated insertion near the blaTEM end
- 5 bp duplicated Tn2 sequence flanks the insertion junction

### Tn1Mer (`Tn1Mer_GQ160960`)

Type: `Tn1`. Subtype: `Tn1Mer`. Source accession: `GQ160960`.

Exact match to the 10491 bp interrupted Tn1 structure in GQ160960. The outer Tn1 grammar is complete, while a defective 5537 bp Tn5036-like transposon interrupts the Tn1 backbone. This is a named interrupted subtype, not the canonical Tn1 reference.

Review status: Sally said to ignore this as the canonical Tn1 because it contains Tn5036-like.

Known differences from the parent definition:

- defective Tn5036-like transposon inserted into the Tn1 backbone
- blaTEM-1 rather than canonical Tn1 blaTEM-2

### Tn3 (`Tn3_V00613`)

Type: `Tn3`. Subtype: `V00613_legacy_9bp_duplication`. Source accession: `V00613`.

Exact match to V00613 is Tn3 with the declared 9 bp duplication. It is a recognised Tn3 sequence variant, not an exact match to canonical HM749966.

Review status: Sally noted the 9 bp duplication and preferred HM749966 as canonical.

Known differences from the parent definition:

- 9 bp duplication relative to canonical Tn3 HM749966

## Related but different elements

Shared components may be detected in another Tn3-family transposon. If the complete Tn1/Tn2/Tn3 grammar is absent and no declared whole-element definition matches, draw the supported components and report an unresolved Tn3-family unit or fragment. Do not force a Tn1, Tn2 or Tn3 name.

- Tn1696
- Tn1721
- Tn5403
