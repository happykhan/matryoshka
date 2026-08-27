#!/usr/bin/env python3
"""Build the exact ISEcp1 reference from non-Tn1/Tn2/Tn3 expert-reviewed sequence."""

from __future__ import annotations

import argparse
from pathlib import Path

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("source", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()

    supplied = next(
        (
            record
            for record in SeqIO.parse(args.source, "fasta")
            if record.id == "ISEcp1_ISfinder"
        ),
        None,
    )
    if supplied is None:
        raise SystemExit("ISEcp1_ISfinder record not found")

    record = SeqRecord(
        supplied.seq,
        id="ISEcp1_FJ621588",
        description=(
            "element_type=IS family=IS1380 name=ISEcp1 "
            "source_accession=FJ621588 reference_position=1..1656 tsd_length=5 "
            "provenance=expert_reviewed_worked_example"
        ),
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write([record], args.output, "fasta")


if __name__ == "__main__":
    main()
