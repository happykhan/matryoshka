# MARA parity specification

This document defines the target for Matryoshka's automated MARA-compatible
annotation. The aim is biological and representational parity with Sally
Partridge's manual practice, not merely a visual imitation of a MARA figure.

## Evidence base

The specification is derived from:

- Partridge and Tsafnat (2018), *Automated annotation of mobile antibiotic
  resistance in Gram-negative bacteria: the Multiple Antibiotic Resistance
  Annotator (MARA) and database*, JAC 73:883-890,
  [doi:10.1093/jac/dkx513](https://doi.org/10.1093/jac/dkx513).
- Partridge (2011), *Analysis of antibiotic resistance regions in
  Gram-negative bacteria*, FEMS Microbiology Reviews 35:820-855,
  [doi:10.1111/j.1574-6976.2011.00277.x](https://doi.org/10.1111/j.1574-6976.2011.00277.x).
- Partridge et al. (2018), *Mobile Genetic Elements Associated with
  Antimicrobial Resistance*, Clinical Microbiology Reviews 31:e00088-17,
  [doi:10.1128/CMR.00088-17](https://doi.org/10.1128/CMR.00088-17).
- Partridge et al. (2009), *Gene cassettes and cassette arrays in mobile
  resistance integrons*, FEMS Microbiology Reviews 33:757-784,
  [doi:10.1111/j.1574-6976.2009.00175.x](https://doi.org/10.1111/j.1574-6976.2009.00175.x).
- Tsafnat, Copty and Partridge (2011), *RAC: Repository of Antibiotic
  resistance Cassettes*, Database 2011:bar054,
  [doi:10.1093/database/bar054](https://doi.org/10.1093/database/bar054).
- Partridge and Hall (2005), *Evolution of transposons containing blaTEM
  genes*, AAC 49:1267-1268,
  [doi:10.1128/AAC.49.3.1267-1268.2005](https://doi.org/10.1128/AAC.49.3.1267-1268.2005).
- Sally Partridge's supplied Tn1/Tn2/Tn3, ISEcp1 TPU and other-transposon
  FASTAs and worked examples, plus her 1 June 2026 comments on the Matryoshka
  demonstration.

## Feature-database contract

Every reference feature must be versioned data rather than a hard-coded drawing
assumption. A feature record needs:

- unique stable feature ID (FID) and unique name;
- feature type and optional biological class;
- exemplar accession and inclusive exemplar coordinates;
- reference sequence and orientation convention;
- feature-specific identity and coverage criteria;
- optional nucleotide/amino-acid constraints needed to distinguish variants;
- expected IR and DR lengths and allowed variants;
- alternative names, phenotype and nomenclature notes;
- provenance, reference and database version.

Exact allele and IS names take priority over family-level calls. A family-level
name is only a fallback and must be labelled as such.

## Annotation grammar

1. Find exact or near-exact complete features first.
2. Reserve their coordinates before resolving fragments and overlaps. This
   implements Sally's recommendation that a complete feature overrides a
   competing fragment at a shared boundary.
3. Permit biologically real overlaps: overlapping ORFs, genes within
   transposons, and short shared boundary sequence.
4. Annotate partial features only when sufficient diagnostic sequence remains.
   The original MARA thresholds were at least 19 bp from a feature end or 32 bp
   from its middle; Matryoshka may use stricter feature-specific thresholds when
   necessary to control false positives.
5. Record which end is present, which end is truncated, and the cause when it
   can be explained by a contig boundary or insertion of another feature.
6. Report unexplained adjacent fragments as a possible divergent feature or
   assembly problem, not as a confident reconstruction.
7. Find DR at element boundaries using the feature-specific expected length.
   Distinguish confirmed sequence, expected-but-unavailable and absent evidence.
8. Keep an interruption and both surviving pieces of its host feature in one
   nested structure when collinearity supports that interpretation.
9. Preserve every call's absolute sequence coordinates, orientation, reference
   evidence and confidence in the table and machine-readable outputs.

## Biological feature vocabulary

The exhaustive raw-component ontology, compound assembly rules and named
source examples now live in
[`mara-component-inventory.md`](mara-component-inventory.md) and the
machine-readable [`mara_component_catalog.yaml`](../matryoshka/mara_component_catalog.yaml).
The table below remains the concise parity checklist.

| Feature | Required interpretation | MARA diagram symbol | Current state |
| --- | --- | --- | --- |
| Resistance gene | Exact allele/framework where nomenclature permits; generalized family only for unresolved or truncated copies | black directional arrow | Partial: AMRFinder allele calls supported and generic OXA calls are explicitly marked unresolved; allele validation and framework database incomplete |
| Gene cassette | Whole cassette including its attC context, not merely the resistance ORF | pale-blue box | Supported from IntegronFinder ORF plus the nearest oriented attC; ORF-only calls remain explicit fragments and an exact cassette FID library is incomplete |
| attC | Cassette recombination site, including targeted interruptions | boundary/site mark | Partial |
| Integron 5'-CS and 3'-CS | Conserved segments as independent complete or truncated features | orange boxes | Exact reference scanning and partial status supported, including the AF071413 5'-CS GGG variant and partial 3'-CS; broader variant library incomplete |
| IRi and IRt | 25 bp ends of class 1 In/Tn, with integrase and tni-end identity preserved | flags and matching DR lollipops | Exact short-reference scanning and orientation-aware rendering supported; inference remains conservative when one end is absent |
| Tn402 tni region | tniABQ, res and tniR component, complete or truncated | purple unit-transposon box | Complete and partial supplied-reference detection supported; complete Tn402 is emitted only with IRi, IRt and a complete tni region |
| Insertion sequence | Exact IS name/isoform, orientation defined by transposase direction and IRR | white block arrow pointing to IRR | Exact support incomplete; family-level ISEScan is not sufficient for structural inference |
| Group II intron | Exact intron, orientation defined by reverse transcriptase | grey block arrow | Pending |
| Unit transposon | Complete named backbone and IR, with internal cargo/interruptions nested | coloured box with IR flags | Supported with curated internal maps for Tn1/Tn2/Tn3, Tn21, Tn1721, Tn1722, Tn4401, Tn5393 and Tn5403; minor substitutions and internal indels assemble into one locus |
| Composite transposon | Same/allowed IS copies plus intervening region and supporting outer DR | enclosing structure plus constituent white IS arrows | Partial; rules must require exact IS identity and orientation |
| IS26 pseudo-compound transposon | Directly oriented IS26 pair delimiting an individually mobile region; outer DR evidence stated | constituent IS26 arrows and bounded region | Corrected: individual candidates only; no generic merged island |
| IS26 translocatable unit | One IS26 plus the specific adjacent segment that can move; never inferred from opposite-orientation bounding copies | one IS26 plus adjacent cargo region | Low-confidence prototype; boundary reconstruction pending |
| ISEcp1 TPU | ISEcp1 from IRL through IRR and captured segment to IRalt, with outer 5 bp DR when intact | ISEcp1 arrow plus captured region and DR | Partial: seven supplied complete CMY-like TPU exemplars are reference-scannable; general IRalt reconstruction and internal component maps remain pending |
| ISCR capture | oriIS-to-terIS rolling-circle feature and captured adjacent region; no conventional IR/DR | directional rolling-circle element plus captured region | Partial motif/rule support |
| IR fragment | Named left/right end with completeness and interrupted state | flag, complete or truncated | Tn1/2/3 supported; general library pending |
| DR/TSD | Feature-specific length and exact repeated sequence, with evidence state | same-colour lollipops | Supported, but all feature classes must use it consistently |
| Unknown insertion/gap | Unclassified sequence between known collinear pieces | white jagged interruption or dashed red gap if >50 bp | Supported for Tn1/2/3; generalization pending |

## Rendering contract

- Draw approximately to scale within a readable locus viewport.
- Use one line unless the locus cannot remain legible; do not shrink an entire
  chromosome into one unreadable line.
- Show unannotated gaps over 50 bp as dashed red lines with their exact length.
- Use directional gene arrows.
- Draw complete gene cassettes as pale-blue boxes and integron conserved
  segments as orange boxes.
- Draw IS as white block arrows, with the pointed end marking IRR.
- Draw Group II introns as equivalent grey block arrows.
- Draw common unit transposons with stable colours and their IR as flags whose
  flat side marks the outer boundary.
- Use a jagged edge only on the truncated side of a partial feature.
- Draw each confirmed DR pair as same-colour lollipops and print its sequence in
  the table.
- Preserve overlapping and nested features rather than flattening them.
- The companion table must include Position, directional Name, FID, Type and
  Notes, including truncation, identity, IR/DR and evidence limitations.

## Expert corrections that are regression requirements

- Do not infer an IS26 island or translocatable unit from oppositely oriented
  IS26 copies. Do not merge a complex IS26-rich region into one object.
- Do not interpret a generic `IS6` call as IS26. Exact identity is required.
- Distinguish full IS26 from IS26 fragments and retain their orientation.
- Annotate `ISKpn7`, not merely `IS21`, and `ISKpn6` within Tn4401. Tn4401 in
  pKpQIL is itself inserted in ISKpn31.
- Annotate `IS1999`, not merely `IS4`, in Tn1999.
- Require the exact `blaOXA` allele; `blaOXA` alone is not an acceptable final
  annotation because OXA families are not a close single group.
- Treat IRi as the integron boundary at the integrase end. Do not end a class 1
  integron automatically at the end of qacEdelta1; the 3'-CS includes sul1 and
  may be extended by downstream modules.
- Retain truncated genes such as catB3 and explain the interrupting IS26.
- Prefer complete Tn5403 boundary sequence over an overlapping Tn1721 fragment;
  then annotate the remaining Tn1721 fragment.
- Keep the ISEcp1 region adjacent to Tn1721 separate from Tn1721 itself.
- For rearranged ISEcp1-blaCTX-M-15 structures, retain individual recognizable
  ISEcp1 and Tn2 fragments and any original DR evidence rather than inventing a
  complete TPU.

## Acceptance corpus

The following supplied sequences are mandatory regression fixtures:

- the three canonical Tn1/Tn2/Tn3 sequences and minor/interrupted synthetic
  derivatives;
- all ISEcp1 TPU CMY-2-like examples, including the 3078+503 captured segment;
- Tn21, Tn1721/Tn1722, Tn4401 and its component regions, Tn5393 and IRR variant,
  and Tn5403 from `other_Tn.fasta`;
- pKpQIL, pEK499, pOXA-48a and the R100 MARA example when complete source
  sequences and expert coordinate truth are available.

For every fixture, tests must check detected identity, coordinates, orientation,
completeness, nesting, DR/IR evidence, table rows and rendered symbols.
