#!/usr/bin/env python3
"""Fetch the small public plasmid acceptance corpus from NCBI."""

from __future__ import annotations

import argparse
import hashlib
import urllib.parse
import urllib.request
from pathlib import Path

from Bio import SeqIO

ACCESSIONS = {
    "pKpQIL": "NC_014016",
    "pOXA-48a": "JN626286",
    "pEK499": "EU935739",
}


def fetch(accession: str) -> bytes:
    query = urllib.parse.urlencode({
        "db": "nuccore",
        "id": accession,
        "rettype": "fasta",
        "retmode": "text",
    })
    request = urllib.request.Request(
        f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?{query}",
        headers={"User-Agent": "matryoshka-reference-fetch/0.1"},
    )
    with urllib.request.urlopen(request, timeout=60) as response:
        return response.read()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--out",
        type=Path,
        default=Path(__file__).parents[1] / "tests" / "test-data" / "acceptance",
    )
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    checksums = []
    metadata = ["file\taccession\trecord_id\tlength"]
    for name, accession in ACCESSIONS.items():
        path = args.out / f"{name}.fasta"
        # EFetch may add one or more blank lines; pin one canonical final newline.
        data = fetch(accession).rstrip() + b"\n"
        path.write_bytes(data)
        record = SeqIO.read(path, "fasta")
        if accession not in record.id:
            raise ValueError(f"NCBI returned {record.id!r} for {accession}")
        checksums.append(f"{hashlib.sha256(data).hexdigest()}  {path.name}")
        metadata.append(f"{path.name}\t{accession}\t{record.id}\t{len(record.seq)}")
    (args.out / "SHA256SUMS").write_text("\n".join(checksums) + "\n", encoding="utf-8")
    (args.out / "METADATA.tsv").write_text("\n".join(metadata) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
