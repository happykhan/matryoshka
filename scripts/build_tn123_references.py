#!/usr/bin/env python3
"""Rebuild the reviewed Tn1/Tn2/Tn3 reference FASTA from public accessions."""

from __future__ import annotations

import argparse
import io
from dataclasses import dataclass
from pathlib import Path

import requests
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
CANONICAL_IDS = ("Tn1_NC_008357", "Tn2_AY123253", "Tn3_HM749966")


@dataclass(frozen=True)
class VariantSource:
    reference_id: str
    accession: str
    start: int
    end: int
    reverse_complement: bool
    expected_length: int


VARIANTS = (
    VariantSource("Tn2c_HM749967", "HM749967.1", 1287, 6236, True, 4950),
    VariantSource("Tn2_1_CP028717", "CP028717.1", 125891, 134869, True, 8979),
    VariantSource("Tn1Mer_GQ160960", "GQ160960.2", 28, 10518, True, 10491),
    VariantSource("Tn3_V00613", "V00613.1", 1, 4957, True, 4957),
)


def _fetch(accession: str) -> SeqRecord:
    response = requests.get(
        EFETCH,
        params={
            "db": "nuccore",
            "id": accession,
            "rettype": "gbwithparts",
            "retmode": "text",
        },
        timeout=60,
    )
    response.raise_for_status()
    return SeqIO.read(io.StringIO(response.text), "genbank")


def _extract(source: VariantSource) -> SeqRecord:
    record = _fetch(source.accession)
    sequence: Seq = record.seq[source.start - 1:source.end]
    if source.reverse_complement:
        sequence = sequence.reverse_complement()
    if len(sequence) != source.expected_length:
        raise ValueError(
            f"{source.reference_id}: expected {source.expected_length} bp, got {len(sequence)}"
        )
    return SeqRecord(
        sequence,
        id=source.reference_id,
        description=f"source_accession={source.accession} extracted={source.start}..{source.end}",
    )


def build(source_fasta: Path, output_fasta: Path) -> None:
    existing = {
        record.id: record
        for record in SeqIO.parse(source_fasta, "fasta")
        if record.id in CANONICAL_IDS
    }
    missing = set(CANONICAL_IDS) - set(existing)
    if missing:
        raise ValueError(f"canonical records missing from {source_fasta}: {sorted(missing)}")
    records = [existing[record_id] for record_id in CANONICAL_IDS]
    records.extend(_extract(source) for source in VARIANTS)
    output_fasta.parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(records, output_fasta, "fasta")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    default = Path(__file__).parents[1] / "matryoshka" / "references" / "tn1_tn2_tn3.fasta"
    parser.add_argument("--source", type=Path, default=default)
    parser.add_argument("--out", type=Path, default=default)
    args = parser.parse_args()
    build(args.source, args.out)


if __name__ == "__main__":
    main()
