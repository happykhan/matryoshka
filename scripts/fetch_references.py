"""
fetch_references.py — Download reference plasmids from NCBI for validation.

Plasmids:
  NC_014016  pKpQIL    K. pneumoniae, blaKPC-2 on Tn4401 (Tn3, 5bp TSD)
  JN626286   pOXA-48a  K. pneumoniae, blaOXA-48 on Tn1999.2
  EU935739   pEK499    E. coli, blaCTX-M-15 mobilised by ISEcp1, IS26-heavy
"""

from __future__ import annotations

import time
from pathlib import Path

from Bio import Entrez, SeqIO

Entrez.email = "nabil@happykhan.com"

ACCESSIONS = {
    "NC_014016": "pKpQIL",
    "JN626286": "pOXA-48a",
    "EU935739": "pEK499",
}

OUTPUT_DIR = Path(__file__).parent.parent / "data" / "reference_plasmids"


def fetch_plasmid(accession: str, name: str, outdir: Path) -> Path:
    outfile = outdir / f"{name}.fasta"
    if outfile.exists():
        # Quick size check — reject anything under 10 kb (probably a gene hit)
        size = outfile.stat().st_size
        if size > 10_000:
            print(f"  {name} ({accession}) already cached ({size:,} bytes) — skipping")
            return outfile
        else:
            print(f"  {name} cached file too small ({size} bytes) — re-fetching")

    print(f"  Fetching {name} ({accession}) from NCBI...")
    handle = Entrez.efetch(
        db="nucleotide",
        id=accession,
        rettype="fasta",
        retmode="text",
    )
    record = SeqIO.read(handle, "fasta")
    handle.close()

    # Rename the record for clarity
    record.id = accession
    record.description = f"{name} | {record.description}"

    SeqIO.write(record, outfile, "fasta")
    print(f"  Saved {len(record.seq):,} bp → {outfile}")
    return outfile


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    print(f"Fetching reference plasmids to {OUTPUT_DIR}")

    for accession, name in ACCESSIONS.items():
        fetch_plasmid(accession, name, OUTPUT_DIR)
        time.sleep(0.5)  # NCBI rate limit

    print("Done.")


if __name__ == "__main__":
    main()
