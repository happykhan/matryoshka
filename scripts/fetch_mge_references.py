"""
fetch_mge_references.py — Download MGE reference sequences from NCBI.

Fetches the paper-cited accessions (Partridge et al. 2018, PMC6148190)
and writes them into matryoshka/references/ with headers containing the
metadata that reference_scan.py uses to map hits to MGEFeatures.

Requires: network + biopython. Set NCBI_EMAIL for Entrez.
"""

from __future__ import annotations

import os
import sys
import time
from pathlib import Path

from Bio import Entrez

REFERENCES_DIR = Path(__file__).resolve().parents[1] / "matryoshka" / "references"


# (output_filename, records) where records = (accession, (start,end)|None, metadata)
REFERENCES: list[tuple[str, list[tuple[str, tuple[int, int] | None, dict[str, str]]]]] = [
    (
        "tn4401.fasta",
        [
            ("EU176011", None,
             {"name": "Tn4401b", "element_type": "transposon",
              "family": "Tn4401", "variant": "b"}),
        ],
    ),
    (
        "tn1546.fasta",
        [
            ("M97297", None,
             {"name": "Tn1546", "element_type": "transposon",
              "family": "Tn1546", "cargo": "vanA"}),
        ],
    ),
    (
        "tn21.fasta",
        [
            ("AF071413", None,
             {"name": "Tn21", "element_type": "transposon",
              "family": "Tn21", "cargo": "class1_integron+mer"}),
        ],
    ),
    (
        "tn1331.fasta",
        [
            ("AF479774", None,
             {"name": "Tn1331", "element_type": "transposon",
              "family": "Tn1331"}),
        ],
    ),
    (
        "tn5393.fasta",
        [
            ("AF262622", None,
             {"name": "Tn5393", "element_type": "transposon",
              "family": "Tn5393", "cargo": "strAB"}),
        ],
    ),
    (
        "acinetobacter_islands.fasta",
        [
            ("CP012952", None,
             {"name": "Tn6022_context", "element_type": "transposon",
              "family": "Tn6022"}),
            ("JN247441", None,
             {"name": "Tn6022_delta", "element_type": "transposon",
              "family": "Tn6022"}),
            ("KU744946", None,
             {"name": "Tn6172_context", "element_type": "transposon",
              "family": "Tn6172"}),
            ("FJ172370", None,
             {"name": "Tn6019", "element_type": "transposon",
              "family": "Tn6019"}),
        ],
    ),
    (
        "mcr1_exemplars.fasta",
        [
            ("CP016184", None,
             {"name": "mcr1_plasmid_single_ISApl1",
              "element_type": "plasmid",
              "note": "mcr-1_with_upstream_ISApl1_only"}),
        ],
    ),
    (
        "integron_archetypes.fasta",
        [
            ("U67194", None,
             {"name": "In2", "element_type": "integron",
              "family": "class1_integron", "subtype": "In2"}),
        ],
    ),
    (
        "gi_sul2.fasta",
        [
            ("KX709966", None,
             {"name": "GIsul2", "element_type": "genomic_island",
              "family": "GIsul2", "cargo": "sul2+tet(31)"}),
        ],
    ),
    (
        "isecp1.fasta",
        [
            ("FJ621588", None,
             {"name": "ISEcp1", "element_type": "IS",
              "family": "IS1380"}),
        ],
    ),
    (
        "tn7.fasta",
        [
            ("AP002527", None,
             {"name": "Tn7", "element_type": "transposon",
              "family": "Tn7"}),
        ],
    ),
    (
        "tn552.fasta",
        [
            ("X52734", None,
             {"name": "Tn552", "element_type": "transposon",
              "family": "Tn552"}),
        ],
    ),
]


# Hardcoded local reference — Tn3-family res site consensus (~120bp).
# Source: consensus derived from Tn3/Tn21/Tn4401 res alignments.
TN3_RES_SITE_FASTA = """\
>Tn3_res_site element_type=res_site family=Tn3 name=Tn3_res_site note=consensus
TTATAAATTATCTCTGGCGGTGTTGACATAAATACCACTGGCGGTGATACTGAGCACATC
AGCAGGACGCACTGACCACCATGAAGGTGACGCTCTTAAAAATTAAGCCCTGAAGAAGGG
"""


def write_res_site_reference() -> None:
    REFERENCES_DIR.mkdir(parents=True, exist_ok=True)
    path = REFERENCES_DIR / "tn3_res_site.fasta"
    path.write_text(TN3_RES_SITE_FASTA)
    print(f"  wrote {path} (local motif)")


def format_header(acc: str, meta: dict[str, str]) -> str:
    name = meta.get("name", acc)
    kv = " ".join(f"{k}={v}" for k, v in meta.items() if k != "name")
    return f"{name} {acc} {kv}"


def fetch_records(email: str, force: bool = False) -> None:
    REFERENCES_DIR.mkdir(parents=True, exist_ok=True)
    Entrez.email = email

    write_res_site_reference()

    for filename, records in REFERENCES:
        out_path = REFERENCES_DIR / filename
        if out_path.exists() and not force:
            print(f"  skip {filename} (exists — use --force to refresh)")
            continue
        blocks: list[str] = []
        for acc, region, meta in records:
            try:
                print(f"  fetching {acc} …", end=" ", flush=True)
                kwargs = {"db": "nuccore", "id": acc,
                          "rettype": "fasta", "retmode": "text"}
                if region:
                    kwargs["seq_start"] = region[0]
                    kwargs["seq_stop"] = region[1]
                handle = Entrez.efetch(**kwargs)
                raw = handle.read().strip()
                handle.close()
                if not raw.startswith(">"):
                    print("no FASTA returned")
                    continue
                # Rewrite header with our metadata
                header_line, *seq_lines = raw.splitlines()
                seq_lines = [l for l in seq_lines if not l.startswith(">")]
                new_header = ">" + format_header(acc, meta)
                blocks.append(new_header + "\n" + "\n".join(seq_lines))
                print(f"{sum(len(l) for l in seq_lines):,} bp")
                time.sleep(0.34)  # NCBI rate limit
            except Exception as e:
                print(f"FAILED: {e}")
        if blocks:
            out_path.write_text("\n".join(blocks) + "\n")
            print(f"  wrote {out_path}")


def main() -> int:
    email = os.environ.get("NCBI_EMAIL", "nabil@happykhan.com")
    Entrez.email = email
    force = "--force" in sys.argv
    fetch_records(email, force=force)
    print("Done.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
