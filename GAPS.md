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

## Next-session priorities

1. Ship Tn4401 a–h references (download the known variants from NCBI / ISfinder).
2. Model the IS26 single-copy translocatable unit as a low-confidence alongside the pair-merged island.
3. Investigate a lightweight plasmid replicon caller (PlasmidFinder FASTA vs BLAST) as a drop-in replacement for mobsuite.
4. Add a confidence-score field per element based on evidence stack (BLAST + TSD + IR + inference rule).
5. Add ter-site / origin-of-transfer motif library for rolling-circle elements.
