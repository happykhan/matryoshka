#!/usr/bin/env python3
"""Build MARA class-1-integron component references from Sally's corpus."""

from __future__ import annotations

import argparse
from pathlib import Path

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

RECORDS = {
    "segments": {
        "3-CS_2239bp_U12338": (
            "integron_segment", "class1_integron", "3'-CS", "MARA-SP-3CS", "U12338",
        ),
        "5-CS_GGG_AF071413": (
            "integron_segment", "class1_integron", "5'-CS", "MARA-SP-5CS-GGG", "AF071413",
        ),
        "5-CS_std": (
            "integron_segment", "class1_integron", "5'-CS", "MARA-SP-5CS-standard", "MARA",
        ),
    },
    "ends": {
        "IRi_AF071413": ("IR", "Tn402", "IRi", "MARA-SP-IRi", "AF071413"),
        "IRt_AF071413": ("IR", "Tn402", "IRt", "MARA-SP-IRt", "AF071413"),
    },
    "tni": {
        "Tn402tni_U67194": (
            "transposon_component", "Tn402", "tni", "MARA-SP-Tn402-tni", "U67194",
        ),
    },
}


def _record(source: SeqRecord, metadata: tuple[str, str, str, str, str]) -> SeqRecord:
    element_type, family, name, fid, accession = metadata
    return SeqRecord(
        source.seq,
        id=source.id,
        description=(
            f"element_type={element_type} family={family} name={name} fid={fid} "
            f"source_accession={accession} provenance=Sally_Partridge "
            "feature_db_version=2026-08-26"
        ),
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("source", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    supplied = {record.id: record for record in SeqIO.parse(args.source, "fasta")}
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for group, mapping in RECORDS.items():
        records = [_record(supplied[record_id], metadata) for record_id, metadata in mapping.items()]
        SeqIO.write(records, args.output_dir / f"mara_integron_{group}.fasta", "fasta")


if __name__ == "__main__":
    main()
