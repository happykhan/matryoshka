# locus component inventory and assembly grammar

This is the biological contract for Matryoshka's locus-compatible annotation.
It separates three questions that must not be conflated:

1. **What raw feature is present?** Exact genes, sites, repeats, boundaries and
   small units are identified independently.
2. **What larger structure do those features justify?** A named compound
   element is emitted only when its required components, order, orientation and
   boundaries satisfy an assembly rule.
3. **How is the evidence drawn?** Every accepted component has a stable symbol;
   partial, interrupted and uncertain evidence stays visible.

The machine-readable source of truth is
[`matryoshka/locus_component_catalog.yaml`](../matryoshka/locus_component_catalog.yaml).
Run `matryoshka catalog --format tsv` for the component list or
`matryoshka catalog --format json` for the components, grammar and named
validation examples together.

## Raw component classes

| Layer | Components | Defining evidence | locus-compatible symbol |
| --- | --- | --- | --- |
| Genes | exact resistance allele, unresolved/truncated resistance framework, non-resistance gene, transposition gene, integrase | allele/reference or functional model; coding direction retained | black directional arrow |
| Cassettes | complete gene cassette, cassette remnant, attC, cassette array | usually one ORF (occasionally two) plus its cognate attC; tandem order and orientation retained | pale-blue box, `#` for remnant, site mark, array outline |
| Integron platform | attI, Pc/P2, 5'-CS, 3'-CS | class-specific site/promoter motifs and conserved-segment references | orange sites/flags and orange blocks |
| Element boundaries | IRL, IRR, IRi, IRt, IRalt, subterminal IR, DR/TSD | element-specific end sequence and positional support; exact DR sequence and evidence state | IR flag; paired same-colour lollipops |
| Transposition machinery | tnpA, tnpR, tniA/B/Q/R, tns genes, res/resI/resII/resIII | curated component sequence or protein model in the expected order | gene arrows and dashed `res` site |
| Small mobile units | exact IS, group-II intron, ISCR | exact/variant reference and boundaries; ISCR uses rcr/oriIS/terIS rather than conventional IRs | white IS arrow, grey intron arrow, pale-yellow ISCR arrow |
| Capture machinery | oriIS, terIS, captured segment, ISEcp1 IRalt | mechanism-specific sites, orientation and outer boundary evidence | site marks and pale captured-region block |
| Structural evidence | interruption, annotation gap, unknown fragment | split collinearity, coordinate gap or non-diagnostic partial match | jagged insertion, red dashed line, `#` fragment |
| Plasmid outer layer | replicon, oriV, oriT, backbone/function genes, ncRNA | PlasmidFinder/curated sites and functional models | teal replicon/sites, gene arrows, blue ncRNA block |

“Atomic” here means a raw annotation object. A gene cassette and an insertion
sequence have internal biology, but they are atomic inputs to the next level of
compound-element assembly because their identities and boundaries must be established
before an integron or transposon can be claimed.

## Compound-element definitions

| Larger element | Required components and relations | Must not be inferred from |
| --- | --- | --- |
| Insertion sequence | exact IS identity, diagnostic/transposase sequence and element-specific ends | family label alone; coincident cargo |
| Gene cassette | one, occasionally two, ORFs plus cognate attC | an ORF without attC evidence |
| Cassette array | ordered, consistently oriented cassettes attached to an attI/Pc platform | nearby unordered resistance genes |
| Minimal integron | intI + attI + Pc; cassette array optional | intI alone |
| Class-2 integron | intI2 or its characteristic truncated form, attI2 and class-2 cassette context; ybeA-A remnant and Tn7 tns context are supporting features | a cassette array without class-2 platform evidence |
| Class-3/4/5 integron | class-specific intI, cognate attI and class-specific flanking evidence | integrase similarity alone |
| Structural class-1 In/Tn | IRi, class-1 platform, and IRt or a supported downstream tni boundary | ending the structure automatically at qacEdelta1; sul1 belongs to the 3'-CS where present |
| Tn402-type unit | 25 bp IRi/IRt, tniA/B/Q, res and tniR; class-1 integron optional | a partial tni match |
| Tn3-family unit | diagnostic backbone, terminal IR pair, tnpA, res and tnpR; cargo optional | a short Tn1/2/3 fragment promoted to an exact member |
| Tn21-subgroup unit | Tn21-like transposition region, res and named terminal ends; mer/integron/internal IS optional | an adjacent module without boundary evidence |
| Hybrid transposon | two or more reference-backed transposon segments in coherent order with supported outer boundaries | nearby fragments forced into one canonical parent |
| Composite transposon | two permitted complete IS copies, family-specific orientation, intervening region and supporting outer DR | different/family-only IS calls or unsupported outer boundaries |
| IS26 pseudo-composite transposon | two exact, complete IS26 copies in the permitted direct orientation and an independently bounded intervening region | opposite IS26 orientations, an IS26 fragment, generic IS6 calls, or a merged “IS26 island” |
| IS26 translocatable unit | one exact IS26 plus the specific adjacent segment, supported TU boundary and mobility orientation | the whole region between opposing outer IS26 copies |
| ISEcp1 TPU | exact ISEcp1 with IRL/IRR, captured segment and IRalt; outer 5 bp DR when intact | interrupted parts promoted to a complete TPU; adjacent Tn1721 absorbed as cargo |
| ISCR capture unit | exact ISCR/rcr, oriIS/terIS relationship and adjacent captured segment | conventional paired IR/DR requirements or unoriented adjacent cargo |
| Multiresistance region | at least two independently annotated resistance/mobile modules plus explicit context/boundaries | an assertion that the entire region moves as one unit |
| Plasmid | replicon or other strong backbone evidence and sequence boundary; may contain oriV/oriT, transfer/maintenance genes, ncRNA and accessory regions | circularity alone or an accessory region alone |

## States and precedence

Every component and compound element can be `complete`, `partial`, `fragment`,
`interrupted` or `uncertain`. The annotation engine follows these rules:

- complete exact features are placed before fragments and overlapping
  alternatives;
- a complete feature may explain and override a competing terminal fragment,
  but the surviving part of the fragment is retained;
- exact IS/allele identity takes priority over a family label;
- interrupted host pieces remain linked to the interrupting feature;
- DR/TSD is recorded as confirmed, expected but unavailable, or absent evidence;
- unexplained fragments remain fragments rather than being forced into a named
  parent element;
- a broad multiresistance region groups modules but does not confer mobility on
  the whole group.

## Named examples in expert-reviewed material

These are validation instances, not new component types.

**Insertion sequences and ISCR:** IS1, IS1R, IS5, IS10, IS26, ISEcp1,
IS1326, IS1353, IS6100, IS1111 family, IS4321, IS5075, ISPa21e, ISKpn6,
ISKpn7, ISKpn14, ISKpn31, IS1999, ISCfr1, IS1182-like, ISCR3-like and
ISCR27 fragments.

**Unit and related transposons:** Tn1, Tn1f, Tn2, Tn2c, Tn3, Tn7, Tn10, Tn21,
Tn1696, Tn1721, Tn1722, Tn4401, Tn5393, Tn5403, Tn1331, Tn402, Tn4352,
Tn1999, Tn2603, Tn5036, Tn5060, Tn6029 and TnPa38.

**Integron and transposon components supplied as sequences:** 5CS_std,
5CS_GGG AF071413, 3CS_2239bp U12338, IRi AF071413, IRt AF071413,
Tn402_tni U67194; Tn21 IRL/IRR, mer and tnp regions; Tn1721 IRL/IRR and
tet region; Tn4401 IRL/IRR, KPC region, right end and tnp region; Tn5393
IRL/IRR/IRRv; and the identical Tn5403 terminal IRs.

**Cassette, resistance and accessory examples:** blaTEM variants, blaKPC,
exact blaOXA alleles including OXA-9, OXA-30, OXA-48 and OXA-181,
blaCTX-M-15, blaCMY-2/23/31/36, blaGES-5, aacA4 variants including C329,
aac(6')-Ib-cr, aac(3)-IIe, aadA1a/aadA2/aadA5, dfrA12/dfrA17,
catA1/catB3, sul1/sul2, qacEdelta1, tet(A)/tet(B), strA/strB aliases,
chrA, mph(A)-mrx-mphR(A), mer, gcuF, ereA3 fragment, qnrB19, qnrE1,
ACC and rmtC. A generic `blaOXA` is explicitly unresolved and cannot be
presented as an allele.

**Group-II intron examples:** Apu1 and Apu2.

## Current implementation boundary

The catalogue is complete enough to define and draw all raw classes above, but
definition/rendering and sequence detection are deliberately reported
separately.

- **Working detection:** canonical/minor/indel/partial Tn1/Tn2/Tn3, with
  independent terminal-IR, blaTEM, tnpR, res and tnpA calls assembled by an
  orientation-aware grammar before exact naming;
  supplied Tn21, Tn1721/Tn1722, Tn4401, Tn5393 and Tn5403 references;
  exact ISKpn6/ISKpn7/IS1999; supplied complete ISEcp1-CMY TPUs; 5'-CS,
  3'-CS, IRi, IRt and Tn402 tni references; IntegronFinder cassette/attC;
  AMRFinder genes; DR/TSD and interruptions where sequence evidence allows.
- **Partial detection:** broader exact IS/allele libraries, general ISEcp1
  IRalt reconstruction, ISCR capture, complete cassette identities, general
  Tn3/Tn21-family variants and IS26 TU boundaries.
- **Defined and drawable, reference/model work still required:** attI and Pc
  variants, group-II introns, oriIS, oriV, oriT, ncRNA and the broad plasmid
  functional categories.

The renderer therefore shows supplied or detected evidence faithfully, but it
does not fabricate a compound element whose grammar has not been satisfied.

## Source traceability

- `Partridge_MARA_2018.pdf`: feature vocabulary, matching precedence,
  complete/partial notation, table fields and visual symbols.
- `Partridge_RAC_2011.pdf`: cassette, attC, cassette-array and integron
  definitions and orientation.
- `Partridge_Tn1-2-3_2005.pdf` and `Tn1_Tn2_Tn3.docx`: Tn1/Tn2/Tn3
  backbone, IRs, res, variants and supplied sequences.
- `OtherTn__examples.docx` and `other_Tn.fasta`: named integron components,
  Tn21, Tn1721/Tn1722, Tn4401, Tn5393 and Tn5403 parts.
- `ISEcp1_TPU_examples.docx`, `TPU_CMY-2-like.fasta` and `TPU_updated.txt`:
  ISEcp1 TPU boundaries, DR, capture examples and additional cargo.
- `Partridge_2011_Fplasmid_blaCTX-M-15.pdf`,
  `Partridge_2011_Fig_S3.ppt` and `R100_PlasAnn_TnCentral_MARA.pptx`:
  nested modules, breakpoints, rearrangements and resistance-region context.
- Expert-review comments PDF: corrections that constrain exact
  identity, boundary choice, precedence and IS26/ISEcp1 containment.
- `Islam_PlasAnn_2026.pdf`: plasmid backbone/function outer layer, kept
  explicitly separate from the core mobile-resistance grammar.
