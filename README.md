# Matryoshka

Detect and annotate the nested hierarchy of mobile genetic elements (MGEs) in complete bacterial genome sequences.

Takes a resolved assembly (complete genome or hybrid assembly) and identifies:
- Plasmid replicons (Inc type, MOB type)
- Transposons (Tn3/Tn21 family, Tn7-like, composites)
- Insertion sequences (IS6/IS26, IS1380/ISEcp1, IS30, IS91/ISCR)
- Integrons (class 1, 2, 3) and gene cassettes
- AMR gene cargo

Crucially, it builds a **containment hierarchy** showing what is nested inside what, and confirms element boundaries using **target site duplications (TSDs)** and **inverted repeats (IRs)**.

## Scope

- Input: complete/resolved FASTA (single chromosome or plasmid)
- Target organisms: *Escherichia coli*, *Klebsiella pneumoniae*
- Output: GFF3, GenBank, JSON hierarchy, SVG diagram

## Reference plasmids for validation

- pKpQIL (*K. pneumoniae*) — blaKPC on Tn4401
- pOXA-48a (*K. pneumoniae*) — blaOXA-48 on Tn1999.2
- pEK499 (*E. coli*) — blaCTX-M-15, IS26-heavy

## Development status

Prototype in progress.
