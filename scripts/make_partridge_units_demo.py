#!/usr/bin/env python3
"""Create the reproducible Sally Partridge unit-transposon demo FASTA."""

from __future__ import annotations

import argparse
from pathlib import Path

from Bio import SeqIO

UNIT_IDS = {
    "Tn1721_AB366441",
    "Tn1722_AB366441",
    "Tn4401_GU595196",
    "Tn5393_AF262622",
    "Tn5403_X75779",
}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    root = Path(__file__).parents[1]
    records = list(SeqIO.parse(root / "matryoshka" / "references" / "tn21.fasta", "fasta"))
    records.extend(
        record
        for record in SeqIO.parse(
            root / "tests" / "test-data" / "partridge-examples" / "other_Tn.fasta",
            "fasta",
        )
        if record.id in UNIT_IDS
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(records, args.output, "fasta")


if __name__ == "__main__":
    main()
