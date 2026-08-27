from pathlib import Path

import pytest
from Bio import SeqIO
from click.testing import CliRunner

from matryoshka.__main__ import cli
from matryoshka.reference_scan import blast_available, scan_reviewed_unit_components
from matryoshka.unit_components import (
    assemble_reviewed_unit_components,
    unit_component_references,
)
from matryoshka.unit_definitions import load_unit_definitions

REFERENCES = Path(__file__).parents[1] / "matryoshka" / "references"


def test_reviewed_unit_definitions_are_complete_and_extractable() -> None:
    document = load_unit_definitions()
    assert set(document["units"]) == {"Tn21", "Tn1721", "Tn1722"}
    assert len(unit_component_references()) == 26
    assert document["units"]["Tn1721"]["required_counts"]["terminal_IR"] == 3
    assert "tet_tnp_junction" in document["units"]["Tn1721"]["required_order"]


def test_reviewed_unit_definitions_can_be_exported_for_review() -> None:
    result = CliRunner().invoke(
        cli,
        ["definitions", "--set", "reviewed-units", "--format", "markdown"],
    )
    assert result.exit_code == 0
    assert "## Tn21" in result.output
    assert "## Tn1721" in result.output
    assert "Required order" in result.output


@pytest.mark.skipif(not blast_available(), reason="BLAST+ is required")
@pytest.mark.parametrize(
    ("filename", "record_id", "expected"),
    [
        ("tn21.fasta", "Tn21", "Tn21-like"),
        ("reviewed_units.fasta", "Tn1721_AB366441", "Tn1721-like"),
        ("reviewed_units.fasta", "Tn1722_AB366441", "Tn1722-like"),
    ],
)
def test_component_only_scan_assembles_reviewed_units(
    filename: str,
    record_id: str,
    expected: str,
) -> None:
    features = scan_reviewed_unit_components(REFERENCES / filename)
    parents = assemble_reviewed_unit_components(features)
    calls = [
        parent.name
        for parent in parents
        if parent.attributes.get("seqid") == record_id
    ]
    assert calls == [expected]


@pytest.mark.skipif(not blast_available(), reason="BLAST+ is required")
def test_related_controls_are_not_promoted_to_reviewed_units() -> None:
    features = scan_reviewed_unit_components(REFERENCES / "reviewed_units.fasta")
    parents = assemble_reviewed_unit_components(features)
    called_contigs = {str(parent.attributes.get("seqid")) for parent in parents}
    assert "Tn4401_GU595196" not in called_contigs
    assert "Tn5393_AF262622" not in called_contigs
    assert "Tn5403_X75779" not in called_contigs


@pytest.mark.skipif(not blast_available(), reason="BLAST+ is required")
def test_minor_sequence_change_and_intergenic_insertion_remain_family_like(
    tmp_path: Path,
) -> None:
    record = next(
        record
        for record in SeqIO.parse(REFERENCES / "reviewed_units.fasta", "fasta")
        if record.id == "Tn1722_AB366441"
    )
    sequence = list(str(record.seq))
    for position in range(250, 5400, 197):
        sequence[position] = "A" if sequence[position] != "A" else "C"
    varied = "".join(sequence[:1850]) + "ACGT" * 30 + "".join(sequence[1850:])
    fasta = tmp_path / "varied.fasta"
    fasta.write_text(f">varied_Tn1722\n{varied}\n")

    parents = assemble_reviewed_unit_components(
        scan_reviewed_unit_components(fasta)
    )
    assert [parent.name for parent in parents] == ["Tn1722-like"]
    assert parents[0].end - parents[0].start + 1 == len(varied)


@pytest.mark.skipif(not blast_available(), reason="BLAST+ is required")
def test_run_bundle_contains_gff3_and_genbank(tmp_path: Path) -> None:
    record = next(
        record
        for record in SeqIO.parse(REFERENCES / "reviewed_units.fasta", "fasta")
        if record.id == "Tn1722_AB366441"
    )
    fasta = tmp_path / "tn1722.fasta"
    SeqIO.write(record, fasta, "fasta")
    out = tmp_path / "results"
    result = CliRunner().invoke(
        cli,
        [
            "run", str(fasta), "--profile", "component-rules",
            "--detectors", "none", "--out", str(out),
        ],
    )
    assert result.exit_code == 0, result.output
    assert (out / "annotation.gff3").read_text().startswith("##gff-version 3")
    records = list(SeqIO.parse(out / "annotation.gbk", "genbank"))
    assert len(records) == 1
    assert any(feature.type == "mobile_element" for feature in records[0].features)
