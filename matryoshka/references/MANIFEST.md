# Reference sequence manifest

Reference FASTAs bundled with Matryoshka. All are derived from public databases
(NCBI GenBank, PlasmidFinder) and redistributed under the source databases'
respective terms — NCBI records are in the public domain; PlasmidFinder DB is
licensed under Apache-2.0 (Camacho / Carattoli et al.).

To re-fetch or update:

```bash
NCBI_EMAIL=you@example.org pixi run python scripts/fetch_mge_references.py --force
```

| File | Description | Source |
|---|---|---|
| `tn3_res_site.fasta` | Tn3-family *res* site consensus (~120bp, V00613 origin) | Local motif (hardcoded) |
| `tn4401.fasta` | Tn4401 a + b variants, clean element regions extracted from GenBank features | EU176011, GU386376, CP069050 |
| `tn1546.fasta` | Tn1546 (vanA vancomycin-resistance cluster) | M97297 |
| `tn21.fasta` | Tn21 (class 1 integron + mer operon) | AF071413 |
| `tn1331.fasta` | Tn1331 (aac(6')-Ib + blaOXA-9 + aadA1) | AF479774 |
| `tn5393.fasta` | Tn5393 (strAB) | AF262622 |
| `tn7.fasta` | Tn7 reference | AP002527 |
| `tn552.fasta` | Tn552 (Staphylococcal beta-lactamase) | X52734 |
| `acinetobacter_islands.fasta` | Tn6022 / Tn6022Δ / Tn6172 / Tn6019 backbones | CP012952, JN247441, KU744946, FJ172370 |
| `mcr1_exemplars.fasta` | mcr-1 plasmid with single upstream ISApl1 | CP016184 |
| `integron_archetypes.fasta` | In2 class 1 integron archetype | U67194 |
| `gi_sul2.fasta` | GIsul2 (sul2 + tet(31)) genomic island | KX709966 |
| `isecp1.fasta` | ISEcp1 reference | FJ621588 |
| `plasmidfinder_enterobacteriales.fasta` | PlasmidFinder Enterobacteriales replicon DB (159 records) | [PlasmidFinder DB](https://bitbucket.org/genomicepidemiology/plasmidfinder_db/) (Apache-2.0) |
| `isaba125.fasta` | ISAba125 (IS30 family) reference — ISEScan fails to detect this element | AY751533 |
| `is26.fasta` | IS26 (IS6 family) reference — ISEScan fails at sequence boundaries | X00011 (via KC964607) |
| `tn3_family_extras.fasta` | Tn3-family extras — Tn3, Tn2, Tn1696, Tn1721, Tn6452, Tn1412 | V00613, AY123253, U12338, X61367, KY807920, L36547 |
| `rolling_circle_ter_sites.fasta` | terIS91, terISCR1, terISCR2, oriIS91 consensus motifs (~30bp) | Local motif (experimental) |

All accessions are cited from Partridge et al. 2018 (Clin Microbiol Rev,
PMC6148190) — the paper whose framework Matryoshka implements.

## Versioning

Current reference set: 2026-04 snapshot. Downloads pinned to these accession
versions. To bump, re-run `fetch_mge_references.py --force` and update this
file with the new date.

## Licensing of redistributed data

- **NCBI GenBank records** (`EU176011`, `GU386376`, etc.) are not subject to
  copyright within the US and are freely redistributable.
- **PlasmidFinder DB** is under Apache-2.0 — see their repository for the
  full licence. Attribution: Carattoli et al. *In Silico Detection and Typing
  of Plasmids using PlasmidFinder and Plasmid Multilocus Sequence Typing*
  (Antimicrob Agents Chemother 2014).
- **Local motifs** (Tn3 res, rolling-circle ter sites) are derived from the
  consensus annotations in the references cited above and are released
  under the same MIT licence as the rest of Matryoshka.
