#!/usr/bin/env python3
"""Build the checksummed machine manifest for bundled reference FASTAs."""

from __future__ import annotations

import argparse
import hashlib
import re
from pathlib import Path

import yaml
from Bio import SeqIO

from matryoshka.reference_scan import (
    REFERENCE_PARAMS,
    VALIDATED_REFERENCE_FILES,
)

DATABASE_VERSION = "2026.08.26-alpha.1"
LOCAL_MOTIFS = {"rolling_circle_ter_sites.fasta", "tn3_res_site.fasta"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    digest.update(path.read_bytes())
    return digest.hexdigest()


def _metadata(path: Path) -> tuple[list[str], list[str]]:
    accessions: set[str] = set()
    record_ids: list[str] = []
    for record in SeqIO.parse(path, "fasta"):
        record_ids.append(record.id)
        accessions.update(re.findall(r"[A-Z]{1,4}_?\d{5,}(?:\.\d+)?", record.description))
        for token in record.description.split():
            if token.startswith("source_accession="):
                accessions.add(token.split("=", 1)[1])
    return sorted(accessions), record_ids


def build(references: Path) -> dict:
    entries = []
    for path in sorted(references.glob("*.fasta")):
        accessions, record_ids = _metadata(path)
        if path.name == "plasmidfinder_enterobacteriales.fasta":
            licence = "Apache-2.0"
            provenance = "PlasmidFinder database"
            redistribution_status = "redistributable"
        elif path.name in LOCAL_MOTIFS:
            licence = "project-generated-derived-consensus"
            provenance = "locally derived motif; see MANIFEST.md"
            redistribution_status = "redistributable_with_project"
        else:
            licence = "source-database-terms"
            provenance = "public accession sequence with curated annotation"
            redistribution_status = "public_accession_traceable"
        entries.append({
            "file": path.name,
            "sha256": sha256(path),
            "bytes": path.stat().st_size,
            "records": len(record_ids),
            "record_ids": record_ids,
            "source_accessions": accessions,
            "profile": "validated" if path.name in VALIDATED_REFERENCE_FILES else "all",
            "status": "validated" if path.name in VALIDATED_REFERENCE_FILES else "experimental_or_legacy",
            "licence": licence,
            "redistribution_status": redistribution_status,
            "provenance": provenance,
            "scan_parameters": REFERENCE_PARAMS.get(
                path.name,
                {"min_identity": 85.0, "min_length": 500},
            ),
        })
    return {
        "schema_version": "1.0",
        "database_version": DATABASE_VERSION,
        "reference_count": len(entries),
        "references": entries,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--references",
        type=Path,
        default=Path(__file__).parents[1] / "matryoshka" / "references",
    )
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    out = args.out or args.references / "manifest.yaml"
    out.write_text(
        yaml.safe_dump(build(args.references), sort_keys=False, allow_unicode=True),
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
