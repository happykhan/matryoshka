#!/usr/bin/env python3
"""Build the exact-IS reference set used by MARA structural rules."""

from __future__ import annotations

import argparse
import urllib.parse
import urllib.request
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def fetch_ncbi_interval(accession: str, start: int, end: int) -> str:
    query = urllib.parse.urlencode({
        "db": "nuccore",
        "id": accession,
        "rettype": "fasta",
        "retmode": "text",
        "seq_start": start,
        "seq_stop": end,
    })
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?{query}"
    with urllib.request.urlopen(url, timeout=30) as response:  # noqa: S310
        lines = response.read().decode().splitlines()
    return "".join(line.strip() for line in lines if not line.startswith(">"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path)
    args = parser.parse_args()

    root = Path(__file__).parents[1]
    supplied = root / "tests" / "test-data" / "partridge-examples" / "other_Tn.fasta"
    supplied_records = {record.id: record for record in SeqIO.parse(supplied, "fasta")}

    records = [
        SeqRecord(
            supplied_records["ISKpn6_ISfinder"].seq,
            id="ISKpn6_EU176011_12118_13657",
            description=(
                "element_type=IS family=IS1182 name=ISKpn6 "
                "source_accession=EU176011 coordinates=12118..13657 strand=reverse "
                "provenance=Sally_Partridge_worked_example"
            ),
        ),
        SeqRecord(
            supplied_records["ISKpn7_GU595196"].seq,
            id="ISKpn7_GU595196",
            description=(
                "element_type=IS family=IS21 name=ISKpn7 "
                "source_accession=GU595196 provenance=Sally_Partridge"
            ),
        ),
        SeqRecord(
            Seq(fetch_ncbi_interval("JN626286", 2704, 4031)),
            id="IS1999_JN626286_2704_4031",
            description=(
                "element_type=IS family=IS4 name=IS1999 "
                "source_accession=JN626286 coordinates=2704..4031 tsd_length=5"
            ),
        ),
    ]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(records, args.output, "fasta")


if __name__ == "__main__":
    main()
