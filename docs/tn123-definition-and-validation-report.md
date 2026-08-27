# Automated Tn1, Tn2 and Tn3 definition and validation report

**Definition version:** 2026-08-27  
**Software:** Matryoshka external alpha  
**Scope:** Tn1, Tn2, Tn3, reviewed close subtypes, minor sequence variants and partial fragments

## Summary

Matryoshka now separates biological definitions from the software that applies them. The Tn1/Tn2/Tn3 rules are maintained in a human-readable YAML file, checked by the definition loader, documented by a JSON Schema, and available as YAML, JSON or a formatted review document from the command line. Adding or revising a type no longer requires adding a new set of conditions to the Python source.

For an arbitrary FASTA record, the program first detects the component sequences, assembles them into a candidate unit, and only then compares the complete locus with reviewed whole-element definitions. It assigns an exact name only when the component grammar is complete and the locus is identical to a declared definition. A previously undeclared sequence with minor changes is reported as `Tn1-like`, `Tn2-like` or `Tn3-like`. A partial match is retained as `Tn1/2/3 fragment`, with its closest reference recorded, but is not promoted to a complete type.

The current validation ledger contains 14 checks and all pass. It covers seven reviewed real-accession definitions, Sally Partridge's two reviewed Tn2 fragments in pEK499, three related Tn3-family elements that must not be misnamed, one complete sequence containing substitutions, and one truncated sequence. A complete result-bundle run on the seven reviewed definitions gives seven independently assembled loci and seven `PASS` proof verdicts.

This establishes a reproducible Tn1/Tn2/Tn3 workflow. It does not claim that the present reference library is an exhaustive definition of every Tn3-family transposon or a complete replacement for expert MARA annotation.

## What happens to an arbitrary sequence

The supported workflow is:

```text
arbitrary FASTA
    |
    +-- scan independently for terminal IRs, blaTEM, tnpR, res and tnpA
    |
    +-- assemble nearby, correctly oriented components into a candidate locus
    |
    +-- check the required component order and both element ends
    |
    +-- compare the complete locus with reviewed whole-element definitions
    |
    +-- assign exact type/subtype, qualified type-like call, unresolved unit,
    |   or fragment
    |
    +-- detect boundary evidence and retain insertions, deletions and mismatches
    |
    +-- write JSON, GFF3, CellGen format, hierarchy SVG, MARA map,
        MARA table and an evidence proof report
```

The component scan and whole-locus comparison are independent pieces of evidence. A reference match alone is not allowed to manufacture the six required components. Conversely, shared components can be drawn even when the whole locus does not match any declared element.

## The group-level biological rule

A complete Tn1/Tn2/Tn3-group candidate requires:

1. two 38 bp terminal inverted repeats;
2. one `blaTEM` gene;
3. one `tnpR` resolvase gene;
4. one `res` resolution site; and
5. one `tnpA` transposase gene.

In the forward element orientation, the required order is:

```text
terminal IR -> blaTEM -> tnpR -> res -> tnpA -> terminal IR
```

The reverse-complement form must have the corresponding reversed order. This grammar establishes membership of the Tn1/Tn2/Tn3 group; it does not distinguish Tn1, Tn2 and Tn3 because their structure and most of their sequence are shared.

## How type and subtype names are assigned

### Exact declared definition

An exact name requires all of the following:

- complete component grammar;
- both terminal ends covered;
- 100% nucleotide identity;
- 100% reference coverage;
- no mismatches;
- no inserted bases; and
- no deleted bases.

This is deliberately strict. One substitution makes the sequence a variant unless the altered sequence has its own reviewed subtype definition.

### Tn1, Tn2 or Tn3 type

A complete candidate is compared with the canonical Tn1, Tn2 and Tn3 whole-locus sequences. The closest type must have at least 95% identity, at least 80% reference coverage, and a score at least 0.5 percentage points better than the next different type. This whole-locus comparison includes the sequence around `res`, where the informative differences described by Partridge and Hall are concentrated.

If the closest complete sequence is not identical to a declared definition but has at least 98% identity, the visible name is qualified as `Tn1-like`, `Tn2-like` or `Tn3-like`. The evidence records its mismatch, insertion and deletion counts.

### Named subtype

A subtype is another declared whole-locus sequence with its own expert rule. The current reviewed set contains:

| Visible call | Parent type | Source accession | Definition |
| --- | --- | --- | --- |
| Tn2c | Tn2 | HM749967 | Tn2 grammar with `blaTEM-1c` and nine substitutions relative to AY123253 |
| Tn2.1 | Tn2 | CP028717 | complete Tn2 backbone interrupted by a 4,024 bp ISEcp1-associated region |
| Tn1Mer | Tn1 | GQ160960 | complete outer Tn1 structure interrupted by a defective Tn5036-like element |
| Tn3 | Tn3 | V00613 | recognised legacy sequence containing the reported 9 bp duplication; HM749966 remains canonical |

Subtype references are not used as fragment bait. Their inserted material can exist elsewhere independently, so a non-canonical subtype reference is considered for naming only when both ends and at least 80% of that subtype definition are covered. Partial-family evidence is assessed against the canonical Tn1, Tn2 and Tn3 references.

### Partial sequence

A match covering at least 400 reference bases is retained, but a missing end or incomplete component grammar prevents a complete name. The visible call is `Tn1/2/3 fragment`; its coordinates, strand, closest type and detected components remain available for review and drawing.

### Related but different element

A different transposon may contain one or more shared genes, IRs or sequence blocks. Those supported parts can still be annotated and drawn. If the complete Tn1/Tn2/Tn3 grammar is absent and no declared whole-element definition matches, the program does not force a Tn1, Tn2 or Tn3 name. This is the expected outcome for examples such as Tn1696, Tn1721 and Tn5403.

## Worked definition: Tn1

Sally selected the Tn1 sequence from accession NC_008357. The human rule is:

> Call exact Tn1 only when the complete six-part grammar is independently supported and the whole locus is a complete, gap-free, mismatch-free match to the declared NC_008357 Tn1 sequence. The element carries `blaTEM-2`.

The corresponding reviewable definition is:

```yaml
types:
  Tn1:
    canonical_reference: Tn1_NC_008357
    source_accession: NC_008357
    bla_allele: blaTEM-2
    components:
      - {role: terminal_IR, name: IRL, start: 1, end: 38, strand: "."}
      - {role: blaTEM, name: blaTEM-2, start: 148, end: 1008, strand: "-"}
      - {role: tnpR, name: tnpR, start: 1191, end: 1748, strand: "-"}
      - {role: res, name: res, start: 1754, end: 1867, strand: "."}
      - {role: tnpA, name: tnpA, start: 1911, end: -34, strand: "+"}
      - {role: terminal_IR, name: IRR, start: -38, end: -1, strand: "."}

definitions:
  Tn1_NC_008357:
    display_name: Tn1
    type: Tn1
    subtype: canonical
    source_accession: NC_008357
    reference_kind: canonical
    component_reference: true
    review_status: Sally-selected canonical sequence
```

Negative coordinates are counted from the right end of the reviewed reference, allowing the same layout convention to work for definitions of different lengths.

## Where the definitions live and how to extend them

The authoritative biological rules are in `matryoshka/tn123_definitions.yaml`. Their structure is documented by `docs/schema/tn123-definitions-v1.schema.json` and checked by the generic loader. The remaining Python code is a component assembler and classifier rather than a separate set of element-specific naming conditions.

The exact file used by the program can be exported for Sally's review:

```bash
matryoshka definitions --format markdown --out tn123-definitions-review.md
matryoshka definitions --format yaml --out tn123-definitions.yaml
matryoshka definitions --format json --out tn123-definitions.json
```

To add a reviewed subtype:

1. add its public reference sequence and provenance;
2. add one entry under `definitions` with the parent type, visible name, accession, review status and expert rule;
3. declare any reviewed inserted or nested features under `additional_features`;
4. rebuild the reference FASTA and checksum manifest;
5. add a positive accession test and a related-element negative test; and
6. rerun the definition validation ledger and complete result-bundle test.

If a future element needs a different component grammar rather than another Tn1/Tn2/Tn3 subtype, it should receive a separate definition family. It should not be squeezed into the Tn1/Tn2/Tn3 rule merely because some components are shared.

## Direct repeats

The five-base direct repeat is boundary evidence, not the primary basis for assigning Tn1, Tn2 or Tn3. Once a candidate boundary is available, Matryoshka examines the flanking sequence immediately outside both ends for the expected repeat. A matching pair is reported with its sequence and coordinates. If the necessary flank is absent, the result says that the repeat is untestable. If the flank is present but no pair is found, the result says that it was searched for and not found.

A short repeated sequence by itself is weak evidence because five-base matches can occur by chance. It therefore supports a component- and reference-backed call but cannot create an exact element name. Natural fragments can retain an original repeat at a truncation junction; that repeat is reported without pretending the missing element is complete.

## Validation results

The reproducible machine-readable results are in `docs/validation/tn123-real-accession-results.tsv` and `.json`.

### Reviewed accession definitions

| Accession | Input definition | Observed call | Result |
| --- | --- | --- | :---: |
| NC_008357 | canonical Tn1 | exact Tn1 | PASS |
| AY123253 | canonical Tn2 | exact Tn2 | PASS |
| HM749966 | canonical Tn3 | exact Tn3 | PASS |
| HM749967 | Tn2c candidate | exact Tn2c | PASS |
| CP028717 | interrupted Tn2.1 | exact Tn2.1, including subtype context | PASS |
| GQ160960 | interrupted Tn1Mer | exact Tn1Mer, including Tn5036-like context | PASS |
| V00613 | legacy Tn3 with 9 bp duplication | exact declared Tn3 subtype, not canonical HM749966 | PASS |

For all seven, whole-locus identity and reference coverage are 100%, the independently detected component grammar is complete, and the proof report passes. Required components are never supplied solely by projection. Reviewed subtype context can be projected into the figure, but it is counted separately from the six sequence-detected components.

### Natural partial fragments: pEK499 (EU935739)

| Coordinates | Observed call | Closest type | Reference coverage | Result |
| --- | --- | --- | ---: | :---: |
| 38,747–40,671 | Tn1/2/3 fragment | Tn2 | 38.89% | PASS |
| 60,316–62,561 | Tn1/2/3 fragment | Tn2 | 45.37% | PASS |

These are retained at the boundaries described in Sally's comments. They are not combined into an invented complete Tn2 and are not labelled exact Tn2.

The proof bundle marks each locus `PARTIAL` and the pEK499 run `PARTIAL_TN123_EVIDENCE`. That status means the incomplete evidence was retained correctly; it is distinct from both a complete-element `PASS` and a validation `FAIL`.

### Minor variation and truncation

The deterministic complete Tn1 test with 25 substitutions is reported as `Tn1-like`, with 99.5% identity, 100% reference coverage and 25 mismatches. The exact Tn1 name is withheld. A Tn1 sequence missing its first 700 bases is reported as `Tn1/2/3 fragment`, closest to Tn1, with 85.86% reference coverage.

### Related elements

Tn1696, Tn1721 and Tn5403 receive no named Tn1, Tn2 or Tn3 call in the Tn1/Tn2/Tn3 validation. Their own supported components or independently defined element references remain available to the main annotation and drawing workflow. A genuine complete Tn2 nested inside another element is still allowed; the rule prevents false naming, not legitimate nesting.

## Outputs used for review

`matryoshka run sequence.fasta --out results` writes:

- `annotation.json`, containing the full hierarchy, definition identifier, subtype, difference counts and evidence;
- `annotation.gff3`;
- `annotation.cell`, preserving the original nested CellGen/Wolvercote representation;
- `hierarchy/*.svg`, the original scale-accurate hierarchy view;
- `mara/*.svg`, one readable locus map per supported target;
- `mara-table/*.svg`, with definition, subtype, identity, coverage, differences and boundary notes;
- `proof/report.html`, a human-readable component-to-call report; and
- `proof/proof.json`, `components.tsv` and `matches.tsv`, the machine-checkable evidence ledger.

Every MARA map and table contains its own key. The maps distinguish sequence-detected components from expert-definition context so that a reader can see what was observed and what was added from a reviewed subtype description.

The repository includes a ready-to-open example at `demo-output/arbitrary-tn123/proof-bundle/proof/report.html`. Its input is one 44,647 bp arbitrary contig containing a Tn1 with 25 substitutions, a Tn2 with an 800 bp internal insertion, and an exact reverse-complement Tn3. All three loci pass the component-to-call proof. The first two remain qualified `Tn1-like` and `Tn2-like` calls with their differences recorded; only the unchanged Tn3 receives an exact canonical name.

## Review questions for Sally

The machinery is now extensible; the remaining expert decisions are content decisions rather than software constraints:

1. confirm whether `Tn2c` is the preferred visible name for HM749967;
2. confirm whether the CP028717 interrupted structure should remain `Tn2.1` and whether the inserted region should be labelled more specifically than “ISEcp1-associated insertion”;
3. confirm whether `Tn1Mer` should be exposed as a named subtype or retained only as an interrupted Tn1 example;
4. confirm the preferred display wording for V00613; and
5. supply or approve further type/subtype definitions and counterexamples for the next expansion.

The YAML definition, generated Markdown review file and accession validation ledger are intended to make those decisions inspectable without requiring Sally to read the implementation code.
