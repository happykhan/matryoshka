#!/usr/bin/env python3
"""Write the reproducible Tn1/Tn2/Tn3 definition-validation ledger."""

from __future__ import annotations

import argparse
import csv
import json
import random
import tempfile
from pathlib import Path
from typing import Any

from Bio import SeqIO

from matryoshka.detect import MGEFeature
from matryoshka.reference_scan import REFERENCES_DIR, scan, scan_tn123_components
from matryoshka.tn123 import assemble_tn123_components

TN123_REFERENCE = REFERENCES_DIR / "tn1_tn2_tn3.fasta"


def _row(
    scenario: str,
    record: str,
    accession: str,
    expected: str,
    feature: MGEFeature | None,
    result: str,
) -> dict[str, Any]:
    attributes = feature.attributes if feature is not None else {}
    return {
        "scenario": scenario,
        "record": record,
        "accession": accession,
        "expected": expected,
        "observed_call": feature.name if feature is not None else "no Tn1/2/3 call",
        "classification_basis": attributes.get("classification_basis", ""),
        "rule_based_type_call": attributes.get("rule_based_type_call", ""),
        "defined_type": attributes.get("defined_type", ""),
        "defined_subtype": attributes.get("defined_subtype", ""),
        "closest_definition": attributes.get("closest_definition", ""),
        "start": feature.start if feature is not None else "",
        "end": feature.end if feature is not None else "",
        "identity_percent": attributes.get("blast_identity", ""),
        "reference_coverage_percent": attributes.get("blast_coverage", ""),
        "mismatch_bases": attributes.get("reference_mismatch_bases", ""),
        "inserted_bases": attributes.get("reference_inserted_bases", ""),
        "deleted_bases": attributes.get("reference_deleted_bases", ""),
        "variant_status": attributes.get("variant_status", ""),
        "result": result,
    }


def _mutate(sequence: str, count: int = 25) -> str:
    bases = list(sequence)
    positions = random.Random(71).sample(range(50, len(bases) - 50), count)
    replacement = {"A": "C", "C": "G", "G": "T", "T": "A"}
    for position in positions:
        bases[position] = replacement[bases[position]]
    return "".join(bases)


def _annotate(query: Path) -> list[MGEFeature]:
    """Run component-rule classification plus secondary reference comparison."""
    features = scan(query, TN123_REFERENCE, 95.0, 100)
    features.extend(scan_tn123_components(query))
    features.extend(assemble_tn123_components(features))
    return [feature for feature in features if feature.element_type == "transposon"]


def build_rows(repo: Path) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    exact_expected = {
        "Tn1_NC_008357": ("NC_008357", "Tn1"),
        "Tn2_AY123253": ("AY123253", "Tn2"),
        "Tn3_HM749966": ("HM749966", "Tn3"),
        "Tn2c_HM749967": ("HM749967", "Tn2c"),
        "Tn2_1_CP028717": ("CP028717", "Tn2.1"),
        "Tn1Mer_GQ160960": ("GQ160960", "Tn1Mer"),
        "Tn3_V00613": ("V00613", "Tn3 (legacy 9 bp duplication subtype)"),
    }
    exact_hits = {
        str(hit.attributes["seqid"]): hit
        for hit in _annotate(TN123_REFERENCE)
    }
    for record, (accession, expected) in exact_expected.items():
        feature = exact_hits.get(record)
        passed = feature is not None and feature.name == expected.split(" (")[0]
        rows.append(_row(
            "reviewed real accession",
            record,
            accession,
            expected,
            feature,
            "PASS" if passed else "FAIL",
        ))

    pek499 = repo / "tests" / "test-data" / "acceptance" / "pEK499.fasta"
    for feature in _annotate(pek499):
        expected = {
            (38_747, 40_671): "generic Tn1/2/3 fragment; closest Tn2",
            (60_316, 62_561): "generic Tn1/2/3 fragment; closest Tn2",
        }.get((feature.start, feature.end), "unexpected call")
        passed = expected != "unexpected call" and feature.name == "Tn1/2/3 fragment"
        rows.append(_row(
            "reviewed natural partial fragments",
            "EU935739.1",
            "EU935739",
            expected,
            feature,
            "PASS" if passed else "FAIL",
        ))

    related = REFERENCES_DIR / "tn3_family_extras.fasta"
    related_hits = _annotate(related)
    for record in ("Tn1696", "Tn1721", "Tn5403"):
        feature = next(
            (hit for hit in related_hits if hit.attributes.get("seqid") == record),
            None,
        )
        passed = feature is None or feature.name == "Tn1/2/3 fragment"
        rows.append(_row(
            "related but different Tn3-family element",
            record,
            "",
            "no named Tn1/Tn2/Tn3 call; supported fragment permitted",
            feature,
            "PASS" if passed else "FAIL",
        ))

    sequences = {
        record.id: str(record.seq)
        for record in SeqIO.parse(TN123_REFERENCE, "fasta")
    }
    with tempfile.TemporaryDirectory() as temporary:
        query = Path(temporary) / "edge-cases.fasta"
        tn1 = sequences["Tn1_NC_008357"]
        query.write_text(
            ">novel_complete_variant\n"
            + _mutate(tn1)
            + "\n>left_truncated_fragment\n"
            + tn1[700:]
            + "\n",
            encoding="utf-8",
        )
        edge_hits = {
            str(hit.attributes["seqid"]): hit
            for hit in _annotate(query)
        }
    rows.append(_row(
        "synthetic minor variation",
        "novel_complete_variant",
        "",
        "Tn1-like, not exact Tn1",
        edge_hits.get("novel_complete_variant"),
        "PASS" if edge_hits.get("novel_complete_variant") and edge_hits[
            "novel_complete_variant"
        ].name == "Tn1-like" else "FAIL",
    ))
    rows.append(_row(
        "synthetic partial sequence",
        "left_truncated_fragment",
        "",
        "generic Tn1/2/3 fragment; closest Tn1",
        edge_hits.get("left_truncated_fragment"),
        "PASS" if edge_hits.get("left_truncated_fragment") and edge_hits[
            "left_truncated_fragment"
        ].name == "Tn1/2/3 fragment" else "FAIL",
    ))
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path(__file__).parents[1] / "docs" / "validation",
    )
    args = parser.parse_args()
    repo = Path(__file__).parents[1]
    rows = build_rows(repo)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0])
    with (args.out_dir / "tn123-real-accession-results.tsv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)
    (args.out_dir / "tn123-real-accession-results.json").write_text(
        json.dumps({
            "definition_version": "2026-08-27",
            "status": "PASS" if all(row["result"] == "PASS" for row in rows) else "FAIL",
            "checks": len(rows),
            "rows": rows,
        }, indent=2) + "\n",
        encoding="utf-8",
    )
    print(f"{len(rows)} checks; {sum(row['result'] == 'PASS' for row in rows)} passed")


if __name__ == "__main__":
    main()
