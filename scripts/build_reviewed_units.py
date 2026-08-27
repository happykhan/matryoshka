#!/usr/bin/env python3
"""Build the expert-reviewed curated whole-unit reference collection."""

from __future__ import annotations

import argparse
from pathlib import Path

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

UNIT_METADATA = {
    "Tn1721_AB366441": ("Tn1721", "Tn1721", "AB366441", 5),
    "Tn1722_AB366441": ("Tn1722", "Tn1722", "AB366441", 5),
    "Tn4401_GU595196": ("Tn4401", "Tn4401", "GU595196", 5),
    "Tn5393_AF262622": ("Tn5393", "Tn5393", "AF262622", 5),
    "Tn5403_X75779": ("Tn5403", "Tn5403", "X75779", 5),
}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("source", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()

    supplied = {record.id: record for record in SeqIO.parse(args.source, "fasta")}
    records: list[SeqRecord] = []
    for record_id, (family, name, accession, tsd_length) in UNIT_METADATA.items():
        source = supplied[record_id]
        fid = f"LOCUS-REF-{record_id}"
        records.append(SeqRecord(
            source.seq,
            id=record_id,
            description=(
                f"element_type=transposon family={family} name={name} "
                f"fid={fid} source_accession={accession} tsd_length={tsd_length} "
                "provenance=expert_reviewed feature_db_version=2026-08-26 "
                "requires_complete=true"
            ),
        ))

    args.output.parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(records, args.output, "fasta")


if __name__ == "__main__":
    main()
