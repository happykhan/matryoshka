#!/usr/bin/env python3
"""Build the exact In2 IRi-to-IRt region from the bundled Tn21 reference."""

from __future__ import annotations

import argparse
from pathlib import Path

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

IN2_START = 4_040
IN2_END = 15_039


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("tn21_reference", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()

    tn21 = next(SeqIO.parse(args.tn21_reference, "fasta"))
    sequence = tn21.seq[IN2_START - 1:IN2_END]
    record = SeqRecord(
        sequence,
        id="In2_AF071413_4040_15039",
        description=(
            "element_type=integron family=class1_integron name=In2 "
            "fid=LOCUS-REF-In2 source_accession=AF071413 "
            "coordinates=4040..15039 ir_length=25 "
            "provenance=expert_reviewed feature_db_version=2026-08-26"
        ),
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write([record], args.output, "fasta")


if __name__ == "__main__":
    main()
