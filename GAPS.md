# Coverage Matrix — matryoshka vs Partridge et al. 2018

Mapping of matryoshka's current capability against the concepts and
elements reviewed in Partridge, Kwong, Firth & Jensen 2018 *Mobile Genetic
Elements Associated with Antimicrobial Resistance* (Clin Microbiol Rev,
PMC6148190).

Legend
- ✅ covered (rule or BLAST-based detection works)
- ⚠️ partial (detected but with limitations below)
- ❌ not covered
- ➖ out of scope (biochemistry/phenotype, not a structural call)

## IS families

| Family | TSD length | Covered? | Notes |
|--------|-----------|----------|-------|
| IS6 (IS26, IS257) | 8 bp | ✅ | ISEScan + boundary confirm + IS26_island composite |
| IS1380 (ISEcp1, ISKpn23) | 5 bp | ✅ | ISEScan + ISEcp1_capture one-ended rule + BLAST hit for ISEcp1 reference |
| IS30 (ISApl1) | 2 bp | ✅ | Tn_mcr1 composite rule |
| IS91 / ISCR (rolling-circle) | 0 bp | ✅ | No TSD expected (family-specific). ISCR_capture rule fires on adjacent AMR |
| IS3 | varies | ⚠️ | Detected by ISEScan; no composite rule |
| IS110, IS200, IS1, IS5, IS66 | varies | ⚠️ | Detected as generic IS features; no family-specific rules |

## Transposon families

| Element | Accession | Covered? | Mechanism | Notes |
|---------|-----------|----------|-----------|-------|
| Tn4401 (blaKPC) | EU176011, CP069050, GU386376 | ✅ | BLAST + IS-flanked rule | Variants **a** and **b** shipped (clean Tn4401 sub-regions extracted from GenBank). Best-variant picker discriminates by subject coverage + HSP count. Variants c/d/e/f/g/h flagged `variant=unknown` with deletion_bp |
| Tn1999 (blaOXA-48) | — | ✅ | IS4-flanked rule | Composite transposon — no res site |
| Tn1546 (vanA) | M97297 | ✅ | vanA signature + BLAST | |
| Tn21 (class 1 integron + mer) | AF071413 | ⚠️ | BLAST only | No signature rule (needs mer/Hg gene calls not in AMRFinder) |
| Tn1331 (multi-AMR) | AF479774 | ✅ | Gene-signature (aac(6')-Ib + blaOXA-9 + aadA1) + BLAST |
| Tn5393 (strAB) | AF262622 | ✅ | BLAST reference hit |
| Tn3 itself | V00613 | ❌ | No explicit reference shipped |
| Tn2 (blaTEM) | AY123253 | ❌ | No rule |
| Tn7 | AP002527 | ⚠️ | BLAST only | No attTn7/glmS detection (needs gene prediction) |
| Tn552 | X52734 | ⚠️ | BLAST only | Staph-specific; no inference rule |
| Tn_mcr1 (ISApl1-mcr-1) | CP016184 etc. | ✅ | Composite rule, 2 bp TSD |
| Tn6022 (AbaR backbone) | CP012952 | ✅ | BLAST |
| Tn6019 (AbaR3) | FJ172370 | ✅ | BLAST → AbaR3 parent island emitted |
| Tn6172 | KU744946 | ✅ | BLAST |
| Tn6021 / Tn6164 | CP012005 | ⚠️ | Not currently downloaded (CP012005 is very large, omitted for speed) |
| Tn6022Δ | JN247441 | ✅ | BLAST |
| Acinetobacter islands (AbGRI1, AbaR3) | — | ✅ | Compound detection: Tn6022+Tn6172 emits AbGRI1 parent; Tn6019 emits AbaR3 parent |

## Integrons

| Concept | Covered? | Notes |
|---------|----------|-------|
| Class 1 integrons (intI1 + attI + cassettes + attC) | ✅ | IntegronFinder parser |
| Class 2 integrons (intI2 with internal stop codon) | ✅ | IntegronFinder parser |
| Class 3 integrons | ✅ | IntegronFinder parser |
| In2 archetype (reference) | ⚠️ | U67194 shipped but no subtype classifier |
| In2-like vs In4-like vs complex | ❌ | Requires 3'-CS structure analysis |
| 3'-CS partial duplications (ISCR1-driven) | ❌ | |
| attC disruptions (IS1111, group II introns) | ❌ | |
| Pc / P2 promoter variants | ❌ | |
| SOS-response regulated intI1 | ➖ | Not structural |
| Gene cassettes in sedentary chromosomal integrons (SCI) | ❌ | |

## Insertion-sequence mechanisms

| Mechanism | Covered? | Notes |
|-----------|----------|-------|
| Composite transposon (two-IS flanked) | ✅ | FLANKED_RULES + IS26_island |
| Translocatable unit (single-IS26 TU) | ⚠️ | `infer_is26_translocatable_units` emits low-confidence TUs for cargo not covered by an IS26 island. Confidence=0.40. Partridge's "50×-preferred" targeting bias not modelled |
| One-ended transposition (ISEcp1) | ✅ | ONE_ENDED_RULES |
| Rolling-circle (IS91 / ISCR) | ✅ | ROLLING_CIRCLE_RULES |
| Cointegrate formation & resolution | ❌ | Not inferrable from structural snapshots |
| RecA-dependent IS26 resolution | ❌ | |
| Programmed frameshift transposase control | ❌ | |

## TSDs / IRs

| Concept | Covered? | Notes |
|---------|----------|-------|
| Family-specific TSD lengths (TSD_LENGTHS dict) | ✅ | IS6=8, IS1380=5, IS30=2, IS91/ISCR=None, Tn3/Tn7/Tn4401/Tn1999=5 |
| TSD boundary offset tolerance (±3bp) | ✅ | `find_tsd` offset window |
| IR detection with mismatch tolerance | ✅ | `find_ir` (2 mismatch default) |
| Hybrid promoter creation (IS26/IS257 -35 + adjacent -10) | ➖ | Regulatory, not structural |
| Terminal IR conservation across Tn3 family (38 bp) | ⚠️ | Not explicitly annotated; BLAST picks up Tn3 backbone |
| Res site (Tn3-family TnpR recognition) | ⚠️ | Expected-position annotation only (not sequence-confirmed). Motif BLAST too divergent to call reliably |

## Plasmids

| Feature | Covered? | Notes |
|---------|----------|-------|
| Replicon typing (PlasmidFinder-style) | ✅ | PlasmidFinder Enterobacteriales DB (159 replicons) shipped — BLAST at ≥95% identity, ≥60% coverage |
| Inc group nomenclature | ✅ | Via PlasmidFinder (IncF/IncHI1/IncHI2/IncN/IncP/IncQ/ColE1 etc.) |
| Replication mode (RC / theta / strand-displacement) | ❌ | |
| Conjugation system (T4SS, MPF types) | ❌ | |
| Relaxase / MOB typing | ❌ | Needs mob_recon (skipped — OOM) |
| oriT sequences | ❌ | |
| Maintenance systems (par / TA / res) | ⚠️ | Tn3-family res annotated; par/TA not |
| Multireplicon plasmids | ❌ | |

## Integrative Conjugative Elements (ICE)

| Feature | Covered? |
|---------|----------|
| ICE presence / integration state | ❌ |
| Site-specific integration at attB | ❌ |
| Reversible excision / circular intermediate | ❌ |

## Resistance islands (genomic context)

| Island | Covered? | Notes |
|--------|----------|-------|
| AbGRI1 (Tn6022 + Tn6172 flanking comM) | ✅ | Compound BLAST detection |
| AbaR3-like (Tn6019 backbone) | ✅ | BLAST |
| GIsul2 (sul2 + tet(31), ISCR2-flanked) | ✅ | BLAST (KX709966) |
| Other genomic islands | ❌ | |

## Output / reporting

| Feature | Covered? |
|---------|----------|
| GFF3 with IR / TSD / parent attributes | ✅ |
| JSON hierarchy with attributes | ✅ |
| GenBank flat file output | ✅ |
| Wolvercote cell-format (circular schematic) | ✅ |
| Scale-accurate linear SVG / PNG | ✅ |
| Per-contig output on multi-FASTA input | ✅ |
| Confidence scores per element | ✅ | 0.0–1.0 scalar + label (high/medium/low/speculative) in `attributes.confidence` |
| HTML / Markdown batch report | ❌ | |
| Rolling-circle ter-site motif library | ⚠️ | `rolling_circle_ter_sites.fasta` scaffolded (terIS91, terISCR1, terISCR2, oriIS91) — experimental, needs refined consensus |

## Known false-positive / false-negative modes

- **ISEcp1 BLAST boundaries** come from a whole-plasmid reference (FJ621588) and can run ~500 bp long. This sometimes places flanking AMR genes inside the ISEcp1 coordinate span; the hierarchy uses a type-based filter to prevent IS → AMR parent nesting, but the visible coordinates remain oversized.
- **Tn4401 variant calling** requires the full set of a–h references. Currently only Tn4401b is shipped; non-b hits are flagged `variant=unknown` with deletion size in bp.
- **Cross-family Tn3 backbone matches** (e.g. Tn6019 hit in E. coli) are suppressed by a 10 kb minimum alignment length on the Acinetobacter references. Tighter thresholds elsewhere may be needed as the reference library grows.
- **Single-IS26 translocatable units** are not emitted. Only IS26/IS26 *pairs* (merged into islands) are detected. A single IS26 adjacent to cargo is reported as two separate features.
- **Plasmid features (MOB, oriT, ICE)** are entirely uncovered pending a memory-safe plasmid-backbone caller. mobsuite crashed the VM and is intentionally disabled.
- **Promoter / regulatory calls** (hybrid promoters, P2 class-1-integron variants, SOS-responsive intI1 expression) are out of scope: matryoshka annotates structural elements, not regulatory behaviour.

## Additional transposon coverage (from 2026-04 paper re-read)

| Element | Accession | Status | Rule / reference |
|---|---|---|---|
| Tn3 (archetype) | V00613 | ✅ | BLAST (`tn3_family_extras.fasta`) |
| Tn2 (blaTEM) | AY123253 | ✅ | BLAST |
| Tn1696 (Tn21 subfamily + mer) | U12338 | ✅ | BLAST |
| Tn1721 (tet(A)) | X61367 | ✅ | BLAST |
| Tn6452 (mcr-5) | KY807920 | ✅ | BLAST |
| Tn6330 (ISApl1-mcr-1, canonical name) | CP016184 | ✅ | Flanked rule (renamed from Tn_mcr1) |
| Tn2006 (ISAba1-blaOXA-23) | — | ✅ | Flanked rule (IS4 + blaOXA-23, 9bp TSD) |
| Tn125 (ISAba125-blaNDM) | — | ✅ | Flanked rule (IS30 + NDM, 3bp TSD) |
| Tn402 / Tn5053 (class-1 integron progenitor) | — | ❌ | Needs tniA/B/Q/R profile — medium complexity |
| Tn916 family (tet(M), vanB2, erm(B)) | — | ❌ | Needs backbone BLAST; Gram+ |
| SXT/R391 ICE (attB in prfC) | — | ❌ | Needs integrase HMM — hard |
| SCCmec (ccr + mecA) | — | ❌ | Needs ccr typing — hard |
| SGI1 / SGI2 (Salmonella) | — | ❌ | attB at 3' trmE — medium |

## Additional IS / IS-family coverage

| Family | Status | Notes |
|---|---|---|
| ISCR1 specifically (vs generic ISCR) | ⚠️ | BLAST finds it, but no 3'-CS-adjacent rule for complex class-1 detection |
| ISCR3/4/5/6/14/15/27 | ❌ | Table 2 AMR associations (blaSPM-1, blaAIM-1, rmtB/D, floR, blaOXA-45) |
| ISAba1 / ISAba125 / ISAba14 / ISAba3 | ⚠️ | ISEScan catches ISAba1; ISAba125 added as BLAST reference (ISEScan misses it). Rules for Tn2006 (OXA-23), Tn125 (NDM) |
| IS1326, IS1353 | ❌ | Diagnostic keys for In2-like vs In4-like integron subtyping |
| IS1111 / IS4321 / IS5075 | ❌ | attC / IRt interruption — needed for broken-integron detection |
| IS256 (erm(B), cfr) | ❌ | Gram+ |
| IS257/IS431 (Staph IS6 subfamily) | ⚠️ | Matched by generic IS6 rule |
| TSD_LENGTHS table expanded | ✅ | Now covers Tn2/Tn21/Tn1696/Tn1721/Tn1546/Tn1331/Tn5393/Tn6452/Tn4401/Tn1999/Tn6330/Tn2006/Tn125/Tn402/Tn5053/Tn552 |

## Integron enhancements

- **Canonical cassette-array string output** ✅ — `attributes.cassette_array` emits `|gene1|gene2|gene3|` on every integron feature. `attributes.cassette_count` gives the count.
- **In2-like vs In4-like vs complex class-1 subtyping** ❌ — requires IS1326/IS1353/ISCR1/IS6100 position analysis
- **Class 4 / 5 mobile integrons** ❌ — rare, low priority

## ICE / Genomic islands entirely uncovered

| Family | Status | Complexity |
|---|---|---|
| SXT/R391 (Vibrio-cholerae-origin, floR/sulII) | ❌ | hard — integrase + attB HMM |
| Tn4371 family ICE (IncP-like T4SS) | ❌ | hard |
| pKLC102 / PAPI-1 / PAGI-n (tRNA-Lys integrators) | ❌ | hard |
| Tn916 family (Tn1549 vanB2 ICE) | ❌ | medium — conserved backbone BLAST |
| SCCmec (MRSA) | ❌ | hard — ccr typing |
| SGI1 / SGI2 | ❌ | medium — attB at trmE |
| ICEberg resource integration | ❌ | would unlock all of the above |

## ISEScan limitations

ISEScan produced no output for 8 of 30 TnCentral benchmark sequences (27%).
Investigation of each case (May 2026 benchmark):

| Sequence | Length | Ground truth IS | Classification | Notes |
|---|---|---|---|---|
| In1 (AY046276) | 6418bp | None | (a) genuine no-IS-found | Integron without IS elements; ISEScan correctly has nothing to report |
| In336 (KP873172) | 3989bp | None | (a) genuine no-IS-found | Integron without IS elements; ISEScan correctly has nothing to report |
| Tn1331 (KC354802.1) | 7996bp | None | (a) genuine no-IS-found | Tn3-family unit transposon; no IS elements in structure |
| Tn1546 (M97297.1) | 10851bp | None | (a) genuine no-IS-found | Tn3-family unit transposon (vanA); no IS elements in structure |
| Tn2 (KT002541) | 4950bp | None | (a) genuine no-IS-found | Tn3-family unit transposon; no IS elements in structure |
| Tn3 (V00613) | 4957bp | None | (a) genuine no-IS-found | Archetype Tn3; no IS elements in structure |
| TnPMLUA4 (KC964607.1) | 4473bp | IS26 @1-820, IS26 @3654-4473 | (c) IS near boundary | Both IS26 copies start at position 1 / end at sequence terminus. ISEScan likely requires flanking context beyond the IS element boundaries to detect the IS. Short sequence (4473bp) may compound the issue |
| Tn125 (JN872328) | 10099bp | ISAba125 @1-1087, ISCR21 @7080-8535, ISAba125 @9013-10099 | (d) ISEScan profile gap | ISAba125 (IS30 family) is not in ISEScan's profile database. Both copies are also at sequence boundaries, compounding the issue. **Mitigated**: ISAba125 added as a BLAST reference in `matryoshka/references/isaba125.fasta` |

**Summary**: Of 8 ISEScan no-output cases, 6 are genuine (no IS elements exist in
the ground truth), 1 is an IS-at-boundary issue (TnPMLUA4), and 1 is a profile
gap (ISAba125 in Tn125, now mitigated by BLAST). The 27% failure rate is misleading
because most "failures" are correct behaviour on sequences without IS elements.

**True ISEScan limitations** (2 of 30 = 7%):
- **IS at sequence boundaries**: When IS copies start at position 1 or end at the
  sequence terminus, ISEScan may not detect them. This affects transposon-only
  TnCentral records where the element *is* the entire sequence. For whole-genome
  input, IS elements will have flanking context and this is not a concern.
- **ISAba125 profile gap**: ISEScan lacks a profile for ISAba125 (IS30 family).
  HMMER-based detection would help here but is not needed now that BLAST reference
  detection covers ISAba125.

**Not worth tackling with HMMER**: The boundary issue (1 case) is specific to
isolated transposon sequences, not real genomic data. The ISAba125 gap (1 case)
is already mitigated by BLAST reference. Neither justifies adding HMMER complexity
at this stage.

## Boundary inaccuracy analysis (9 elements, May 2026)

For 9 elements where detection was correct but boundaries deviated >50bp from
TnCentral annotations. Investigated whether BLAST hit coordinates could be
extended to the nearest TSD to improve boundaries.

**Finding: TSD-based extension is not applicable to any of these 9 cases.**

| Element | Cause | Detail |
|---|---|---|
| Tn5393 | BLAST reference shorter | Reference AF262622 covers only the Tn5393 unit; TnCentral includes flanking IS1133+ISEcp1 (3046bp upstream) not in the reference |
| Tn4401b | Variant discrimination | BLAST picks Tn4401a (7019bp) instead of Tn4401b (10006bp). Variant picker issue, not alignment truncation |
| In0 | Partial reference match | In2 archetype reference (U67194) covers only 2682bp of the 8525bp In0 element |
| In1 | IntegronFinder offset | IntegronFinder reports 1683-6216 instead of 1-6418. Boundary convention mismatch |
| In104 | IntegronFinder miss | Detected as IS26_TU (153-2184) instead of the full 5411bp integron |
| In336 | IntegronFinder offset | IntegronFinder reports 392-3787 instead of 1-3989. Boundary convention mismatch |
| Tn2012 | One-ended capture extent | ISEcp1_capture stops at cargo end (2583); element extends 157bp further with no downstream TSD |
| TnEcp1.1 | One-ended capture extent | ISEcp1_capture stops at cargo end (2656); element extends 761bp further with no downstream TSD |
| Tn2.1 | Complex multi-IS element | ISEcp1_capture covers 130-3392 of a complex 8979bp element with ISEcp1::IS1F structure |

**Root causes by category:**
- **3 BLAST reference mismatches**: Reference is a different size than TnCentral's annotation. Fix: improve references, not boundary logic.
- **3 IntegronFinder boundaries**: IntegronFinder uses different boundary conventions than TnCentral. Structural difference, not error.
- **3 ISEcp1 one-ended capture**: One-ended transposition by definition has no downstream IS/TSD boundary marker. The captured region extends beyond the last cargo gene but we have no structural signal to anchor to.

## Future-session priorities

1. **Tn402/Tn5053 tni module detection** — enables class-1 integron subtype discrimination (In2-like vs In4-like vs complex).
2. **ISAba1 / ISAba125 / ISAba14 specific rules** beyond Tn2006/Tn125 — add a-family-specific TSD table + reference library for Acinetobacter IS types.
3. **Tn916-family backbone BLAST** — tet(M) / vanB2 / erm(B) Gram+ ICE coverage.
4. **Shufflon detection for IncI1/I2 plasmids** — invertible pilV region needs a special representation in MGEFeature.
5. **attB catalogue (prfC, glmS, tRNA-Lys, tRNA-Gly, comM, trmE)** — detect *where* an ICE/island integrated.
6. **SCCmec ccr-complex typing** — highest-cited MRSA MGE, fully missing.
7. **Pull CP012005** once storage budget allows — unlocks Tn6021 / Tn6164 AbaR variants.
8. **ICEberg integration** — mirror / ship a subset of ICEberg as a BLAST reference.
