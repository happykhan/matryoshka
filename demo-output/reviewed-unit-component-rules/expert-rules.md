# Expert rule review

Rule versions: Tn1/Tn2/Tn3 `2026-08-27`; reviewed units `2026-08-27`.

This report is generated directly from the same YAML definitions used by the program. It translates those executable rules into review language; it does not maintain a separate manual description.

## Scope of this review

This edition covers the six families with executable component-driven rules: Tn1, Tn2, Tn3, Tn21, Tn1721 and Tn1722. Other roadmap targets are not included here until their component grammars are implemented; reference-supported names must not be mistaken for reviewed expert rules.

## How a sequence becomes a call

1. **Detect parts independently.** Search for terminal repeats, genes, resolution sites, resistance genes and diagnostic junctions.
2. **Assemble a candidate.** Check copy number, order, orientation and maximum span.
3. **Apply the family rule.** A complete component grammar can create a family or family-like call without a whole-element match.
4. **Assign a type where rules permit it.** Compare the locally detected components with the profiles declared for each type.
5. **Compare with canonical sequences second.** This confirms an exact known definition and records differences; it does not drive discovery.
6. **Keep incomplete evidence.** Components and fragments are drawn even when there is not enough evidence to name a complete element.

## Meaning of result labels

| Label | Meaning |
| --- | --- |
| Exact name, for example **Tn1** | Complete rule plus a 100% complete, mismatch-free match to a declared exact definition |
| **-like** | Complete component rule, but not identical to a declared exact definition |
| Unresolved group member | Complete family structure but insufficient evidence to choose a type |
| Fragment | Some supported parts are present, but required grammar is incomplete |

## Tn1/Tn2/Tn3 group

### Plain-language decision

The software first looks for the parts of the transposon independently. It calls a complete Tn1/Tn2/Tn3-group locus only when all required parts occur in the correct order between two ends. It then compares the local versions of the genes and resolution site to decide whether the locus is most consistent with Tn1, Tn2 or Tn3. A whole-element match is checked afterwards and is not allowed to create the initial family call.

Expected element boundaries: two 38 bp terminal inverted repeats. A 5 bp target-site duplication is supporting boundary evidence, not a substitute for the required internal parts.

### Evidence that must be found

| Component | Minimum number | What it represents |
| --- | ---: | --- |
| `terminal_IR` | 2 | Short terminal inverted repeat marking an element end |
| `blaTEM` | 1 | Beta-lactam resistance gene |
| `tnpR` | 1 | Resolvase gene used during transposition |
| `res` | 1 | Resolution site acted on by the resolvase |
| `tnpA` | 1 | Transposase gene used to move the element |

Required order on the forward strand:

`terminal_IR → blaTEM → tnpR → res → tnpA → terminal_IR`

The reverse-complement order is also accepted.

### How Tn1, Tn2 and Tn3 are separated

The resistance gene is required cargo for this group, but its allele does **not** determine the transposon type. Classification uses the transposition backbone sequences.

| Reported type | Canonical example | Primary backbone evidence | Cargo reported separately |
| --- | --- | --- | --- |
| Tn1 | `NC_008357` | tnpR and tnpA; res shown as context | blaTEM-2 |
| Tn2 | `AY123253` | tnpR and tnpA; res shown as context | blaTEM-1b |
| Tn3 | `HM749966` | tnpR and tnpA; res shown as context | blaTEM-1a |

There is no cross-gene weighting and no averaged type score. Each detected backbone gene is first assigned independently to a declared sequence-profile group. The resulting pair is then looked up in the expert haplotype rule:

| tnpR profile group | tnpA profile group | Rule-based result |
| --- | --- | --- |
| Tn1/Tn3-like tnpR | Tn1/Tn2-like tnpA | Tn1-like |
| Tn2-like tnpR | Tn1/Tn2-like tnpA | Tn2-like |
| Tn1/Tn3-like tnpR | Tn3-like tnpA | Tn3-like |
| Any other confident combination | Any other confident combination | Tn1/Tn2/Tn3-group mosaic |
| Missing or ambiguous | Missing or ambiguous | Unresolved group element |

A gene profile group is called only when its best local match reaches 90% and exceeds the other declared group by at least 0.5 percentage points. These thresholds decide whether an individual gene can be placed in a group; they do not determine the type by averaging genes. res is structural context and blaTEM is reported cargo.

The canonical component sequences support that choice:

| Component | Tn1 vs Tn2 | Tn1 vs Tn3 | Tn2 vs Tn3 | Interpretation |
| --- | ---: | ---: | ---: | --- |
| tnpR | 94.444% | 99.642% | 94.624% | Primarily separates Tn2 |
| tnpA | 99.102% | 97.572% | 97.772% | Primarily separates Tn3 |
| res | 79.130% | 96.522% | 79.825% | Context hotspot; not scored |
| blaTEM | 99.303% | 99.419% | 99.652% | Cargo only; not scored |

### Reviewed exact names and subtypes

| Name | Parent type | Accession | What makes this exact definition distinct |
| --- | --- | --- | --- |
| Tn1 | Tn1 | `NC_008357` | Exact, complete, gap-free match to NC_008357 plus complete component grammar. |
| Tn2 | Tn2 | `AY123253` | Exact, complete, gap-free match to AY123253 plus complete component grammar. |
| Tn3 | Tn3 | `HM749966` | Exact, complete, gap-free match to HM749966 plus complete component grammar. |
| Tn2c | Tn2 | `HM749967` | 9 nucleotide substitutions relative to canonical Tn2 AY123253; blaTEM-1c allele |
| Tn2.1 | Tn2 | `CP028717` | 4024 bp ISEcp1-associated insertion near the blaTEM end; 5 bp duplicated Tn2 sequence flanks the insertion junction |
| Tn1Mer | Tn1 | `GQ160960` | defective Tn5036-like transposon inserted into the Tn1 backbone; blaTEM-1 rather than canonical Tn1 blaTEM-2 |
| Tn3 | Tn3 | `V00613` | 9 bp duplication relative to canonical Tn3 HM749966 |

### What the program will not do

Shared components may be detected in another Tn3-family transposon. If the complete Tn1/Tn2/Tn3 grammar is absent and no declared whole-element definition matches, draw the supported components and report an unresolved Tn3-family unit or fragment. Do not force a Tn1, Tn2 or Tn3 name.

Examples deliberately not forced into a Tn1/Tn2/Tn3 name: Tn1696, Tn1721, Tn5403.

## Tn21

### Plain-language decision

Call a complete Tn21-like unit when two Tn21 terminal IRs delimit a Tn21-specific tnpA-tnpR-res-tnpM transposition region in the declared order. The class 1 integron, partial or complete tni region, internal insertion sequences and mer module are described when present but are not required to discover the Tn21-subgroup backbone. Exact Tn21 requires secondary confirmation against AF071413.

The candidate must span approximately 15,000–35,000 bp and contain 38 bp terminal inverted repeats. A 5 bp target-site duplication may support the inferred boundaries.

### Evidence that must be found

| Component | Minimum number | What it represents |
| --- | ---: | --- |
| `terminal_IR` | 2 | Short terminal inverted repeat marking an element end |
| `tnpA` | 1 | Transposase gene used to move the element |
| `tnpR` | 1 | Resolvase gene used during transposition |
| `res` | 1 | Resolution site acted on by the resolvase |
| `tnpM` | 1 | Tn21-family transposition-region gene |

Required order on the forward strand:

`terminal_IR → tnpA → tnpR → res → tnpM → terminal_IR`

Optional features are drawn when detected but do not create the family call: `integron_5cs`, `cassette_region`, `integron_3cs`, `tni`, `mer_region`.

### Naming decision

Complete component grammar produces **Tn21-like**. The unqualified name **Tn21** requires an exact, complete secondary match to the reviewed example `AF071413`. Missing required components remain visible but are not promoted to a complete family call.

## Tn1722

### Plain-language decision

Call a complete Tn1722-like unit when two Tn1721-family terminal IRs delimit the mcp-tnpR-tnpA backbone in the declared order. Exact Tn1722 requires secondary whole-element confirmation against the 5640 bp reviewed interval from AB366441.

The candidate must span approximately 5,000–8,500 bp and contain 38 bp terminal inverted repeats. A 5 bp target-site duplication may support the inferred boundaries.

### Evidence that must be found

| Component | Minimum number | What it represents |
| --- | ---: | --- |
| `terminal_IR` | 2 | Short terminal inverted repeat marking an element end |
| `mcp` | 1 | Methyl-accepting chemotaxis-like protein gene |
| `tnpR` | 1 | Resolvase gene used during transposition |
| `tnpA` | 1 | Transposase gene used to move the element |

Required order on the forward strand:

`terminal_IR → mcp → tnpR → tnpA → terminal_IR`

No optional components are declared for this rule.

### Naming decision

Complete component grammar produces **Tn1722-like**. The unqualified name **Tn1722** requires an exact, complete secondary match to the reviewed example `AB366441`. Missing required components remain visible but are not promoted to a complete family call.

## Tn1721

### Plain-language decision

Tn1721 consists of the complete Tn1722-like backbone followed by an internal copy of IRR, the tetR-tetA-pecM region and a partial duplication of the tnp region ending at the terminal IRR. Require all three IR occurrences and the unique tet-to-duplicated-tnp junction; tet genes or a Tn1722-like backbone alone are insufficient. Exact Tn1721 requires secondary whole-element confirmation against the 11128 bp reviewed interval from AB366441.

The candidate must span approximately 10,000–16,000 bp and contain 38 bp terminal inverted repeats. A 5 bp target-site duplication may support the inferred boundaries.

### Evidence that must be found

| Component | Minimum number | What it represents |
| --- | ---: | --- |
| `terminal_IR` | 3 | Short terminal inverted repeat marking an element end |
| `mcp` | 1 | Methyl-accepting chemotaxis-like protein gene |
| `tnpR` | 1 | Resolvase gene used during transposition |
| `tnpA` | 1 | Transposase gene used to move the element |
| `tetR` | 1 | Regulator of tetracycline resistance |
| `tetA` | 1 | Tetracycline resistance gene |
| `pecM` | 1 | Membrane-protein gene in the Tn1721 variable region |
| `tet_tnp_junction` | 1 | Sequence junction linking the tet region to duplicated transposition DNA |

Required order on the forward strand:

`terminal_IR → mcp → tnpR → tnpA → terminal_IR → tetR → tetA → pecM → tet_tnp_junction → terminal_IR`

No optional components are declared for this rule.

### Naming decision

Complete component grammar produces **Tn1721-like**. The unqualified name **Tn1721** requires an exact, complete secondary match to the reviewed example `AB366441`. Missing required components remain visible but are not promoted to a complete family call.

## Questions for expert review

For each family, please check:

- Are the required components and their minimum copy numbers correct?
- Is the required order correct on the element-relative forward strand?
- Should any currently optional component be required, or vice versa?
- Are the acceptable candidate-size limits biologically reasonable?
- Is the canonical accession appropriate only for exact naming?
- Which known variants need their own exact subtype definition?
- Which related elements should be explicit negative controls?

## Source YAML

The executable sources are `matryoshka/tn123_definitions.yaml` and `matryoshka/unit_definitions.yaml`. Edit and validate those files rather than editing this generated report.
