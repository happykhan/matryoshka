#!/usr/bin/env python3
"""Build curated ISEcp1 transposition-unit references from Sally's examples."""

from __future__ import annotations

import argparse
from pathlib import Path

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def _cargo_name(record_id: str) -> str:
    fields = record_id.split("_")
    if len(fields) < 3 or fields[0] != "TPU":
        raise ValueError(f"unexpected TPU identifier: {record_id}")
    return fields[1]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("source", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()

    records: list[SeqRecord] = []
    for supplied in SeqIO.parse(args.source, "fasta"):
        if not supplied.id.startswith("TPU_"):
            continue
        cargo = _cargo_name(supplied.id)
        accession = supplied.id.rsplit("_", 1)[-1]
        records.append(SeqRecord(
            supplied.seq,
            id=supplied.id,
            description=(
                "element_type=transposition_unit family=ISEcp1_TPU "
                f"name=ISEcp1-bla{cargo}_TPU cargo=bla{cargo} "
                f"source_accession={accession} tsd_length=5 "
                "orientation=ISEcp1_IRL_to_IRR_to_IRalt "
                "provenance=Sally_Partridge requires_complete=true"
            ),
        ))

    if not records:
        raise SystemExit("no TPU records found")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(records, args.output, "fasta")


if __name__ == "__main__":
    main()
