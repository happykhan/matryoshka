"""
fetch_mge_references.py — Download MGE reference sequences from NCBI.

Fetches the paper-cited accessions (Partridge et al. 2018, PMC6148190)
and writes them into matryoshka/references/ with headers containing the
metadata that reference_scan.py uses to map hits to MGEFeatures.

Some references (Tn4401 variants) are sub-regions within larger plasmids;
for those, we fetch the GenBank record, locate the mobile_element feature
matching a target annotation, and write only the element sequence.

Requires: network + biopython. Set NCBI_EMAIL for Entrez.
"""

from __future__ import annotations

import os
import re
import sys
import time
from io import StringIO
from pathlib import Path

from Bio import Entrez, SeqIO

REFERENCES_DIR = Path(__file__).resolve().parents[1] / "matryoshka" / "references"


# (output_filename, records) where records = (accession, (start,end)|None, metadata)
REFERENCES: list[tuple[str, list[tuple[str, tuple[int, int] | None, dict[str, str]]]]] = [
    # tn4401.fasta is handled separately by fetch_tn4401_variants() — it
    # extracts just the element region from GenBank features rather than
    # shipping whole plasmids.
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
    # Additional Tn3-family references called out in Partridge 2018 Fig 3 legend.
    (
        "tn3_family_extras.fasta",
        [
            # Tn3 archetype — Tn3 itself
            ("V00613", None,
             {"name": "Tn3", "element_type": "transposon",
              "family": "Tn3"}),
            # Tn2 — blaTEM progenitor (hybrid of Tn1/Tn2/Tn3 family)
            ("AY123253", None,
             {"name": "Tn2", "element_type": "transposon",
              "family": "Tn3", "cargo": "blaTEM"}),
            # Tn1696 — Tn21 subfamily
            ("U12338", None,
             {"name": "Tn1696", "element_type": "transposon",
              "family": "Tn3", "cargo": "mer+integron"}),
            # Tn1721 — Tn3 family tet(A) composite
            ("X61367", None,
             {"name": "Tn1721", "element_type": "transposon",
              "family": "Tn3", "cargo": "tetA"}),
            # Tn6452 — mcr-5 on Tn21 subfamily
            ("KY807920", None,
             {"name": "Tn6452", "element_type": "transposon",
              "family": "Tn3", "cargo": "mcr-5"}),
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
            # Tn6021 / Tn6164 — AbaR variants on CP012005 (large A. baumannii chromosome)
            ("CP012005", None,
             {"name": "Tn6021_Tn6164_context",
              "element_type": "transposon",
              "family": "Tn6021",
              "note": "A_baumannii_chromosome"}),
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


# ---------------------------------------------------------------------------
# Tn4401 variant sub-region extractor
# ---------------------------------------------------------------------------

TN4401_VARIANT_SOURCES: list[tuple[str, str, dict[str, str]]] = [
    # (accession, variant-name-to-find-in-feature, metadata)
    # EU176011 contains Tn4401a at 4054..13959 (GenBank annotation).
    ("EU176011",  "Tn4401a",
     {"variant": "a", "deletion_bp": "100",
      "note": "Naas 2008 original, 100bp Pout deletion"}),
    ("GU386376",  "Tn4401a",
     {"variant": "a", "deletion_bp": "100",
      "note": "K. pneumoniae pRYCKPC3.1"}),
    ("KY930325",  "Tn4401d",
     {"variant": "d", "deletion_bp": "68",
      "note": "Kitchel 2010 variant"}),
    ("FJ223607",  "Tn4401e",
     {"variant": "e", "deletion_bp": "255",
      "note": "Leavitt 2010 variant"}),
    ("KC788405",  "Tn4401f",
     {"variant": "f", "deletion_bp": "188",
      "note": "Rodriguez-Bano variant"}),
    # Full-length Tn4401b from a plasmid where it's explicitly labeled
    ("CP069050",  "Tn4401b",
     {"variant": "b", "deletion_bp": "0",
      "note": "Pseudomonas aeruginosa p33Kpn22-KPC, full length"}),
]


def _find_tn4401_region(record, variant_label: str) -> tuple[int, int] | None:
    """Locate a Tn4401 feature in a GenBank record matching variant_label.

    Falls back to any 'Tn4401' mobile_element if the variant-specific
    label isn't present (common when paper names appear only in the
    publication, not the GenBank annotation).
    """
    pattern = re.compile(variant_label, re.IGNORECASE)
    generic = re.compile(r"Tn4401", re.IGNORECASE)

    def _plausible(feat) -> bool:
        length = int(feat.location.end) - int(feat.location.start)
        return 9_000 <= length <= 12_000

    # First pass: variant-specific match
    for feat in record.features:
        if feat.type != "mobile_element":
            continue
        for qkey in ("mobile_element_type", "note", "label", "product"):
            for val in feat.qualifiers.get(qkey, []):
                if pattern.search(val):
                    return int(feat.location.start) + 1, int(feat.location.end)

    # Second pass: any Tn4401 mobile_element feature of plausible length
    for feat in record.features:
        if feat.type != "mobile_element":
            continue
        for qkey in ("mobile_element_type", "note", "label", "product"):
            for val in feat.qualifiers.get(qkey, []):
                if generic.search(val) and _plausible(feat):
                    return int(feat.location.start) + 1, int(feat.location.end)

    # Third pass: any feature of any type mentioning Tn4401 with plausible length
    for feat in record.features:
        if not _plausible(feat):
            continue
        for val in (v for vs in feat.qualifiers.values() for v in vs):
            if generic.search(val):
                return int(feat.location.start) + 1, int(feat.location.end)

    return None


def fetch_tn4401_variants(force: bool = False) -> None:
    """Extract clean Tn4401 regions from GenBank records."""
    out_path = REFERENCES_DIR / "tn4401.fasta"
    if out_path.exists() and not force:
        print(f"  skip tn4401.fasta (exists — use --force to refresh)")
        return

    blocks: list[str] = []
    for acc, variant_label, meta in TN4401_VARIANT_SOURCES:
        try:
            print(f"  fetching {acc} GenBank for {variant_label} …", end=" ", flush=True)
            h = Entrez.efetch(db="nuccore", id=acc, rettype="gb", retmode="text")
            raw = h.read()
            h.close()
            rec = next(SeqIO.parse(StringIO(raw), "genbank"))
            region = _find_tn4401_region(rec, variant_label)
            if not region:
                print(f"no {variant_label} feature — skipping")
                continue
            start, end = region
            subseq = str(rec.seq[start - 1: end])
            header_meta = dict(meta)
            header_meta["source_accession"] = f"{acc}:{start}-{end}"
            header_meta["element_type"] = "transposon"
            header_meta["family"] = "Tn4401"
            header_meta["name"] = variant_label
            # Include the accession in the seq id so makeblastdb -parse_seqids
            # sees unique identifiers even when two records are the same variant.
            unique_id = f"{variant_label}_{acc}"
            kv = " ".join(f"{k}={v}" for k, v in header_meta.items())
            hdr = f">{unique_id} {acc}:{start}-{end} {kv}"
            # Wrap at 60
            wrapped = "\n".join(subseq[i:i+60] for i in range(0, len(subseq), 60))
            blocks.append(f"{hdr}\n{wrapped}")
            print(f"{len(subseq):,} bp ({start}-{end})")
            time.sleep(0.34)
        except Exception as e:
            print(f"FAILED: {e}")

    if blocks:
        out_path.write_text("\n".join(blocks) + "\n")
        print(f"  wrote {out_path}")


def format_header(acc: str, meta: dict[str, str]) -> str:
    name = meta.get("name", acc)
    kv = " ".join(f"{k}={v}" for k, v in meta.items() if k != "name")
    return f"{name} {acc} {kv}"


def fetch_records(email: str, force: bool = False) -> None:
    REFERENCES_DIR.mkdir(parents=True, exist_ok=True)
    Entrez.email = email

    write_res_site_reference()
    fetch_tn4401_variants(force=force)

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
