"""Contract tests for reviewable expert element definitions."""

from __future__ import annotations

import json
from pathlib import Path

from Bio import SeqIO
from Bio.Align import PairwiseAligner
from click.testing import CliRunner

from matryoshka.__main__ import cli
from matryoshka.element_definitions import (
    definitions_as_markdown,
    load_tn123_definitions,
    tn123_reference_metadata,
)


def test_tn123_definitions_cover_types_subtypes_and_review_rules():
    document = load_tn123_definitions()
    assert document["schema_version"] == "1.0"
    assert set(document["types"]) == {"Tn1", "Tn2", "Tn3"}
    assert {
        "Tn1_NC_008357",
        "Tn2_AY123253",
        "Tn3_HM749966",
        "Tn2c_HM749967",
        "Tn2_1_CP028717",
        "Tn1Mer_GQ160960",
        "Tn3_V00613",
    } == set(document["definitions"])
    exact = document["classification"]["exact_definition_match"]
    assert exact["identity_percent"] == 100.0
    assert exact["require_zero_mismatches"] is True
    assert exact["require_zero_inserted_bases"] is True
    assert exact["require_zero_deleted_bases"] is True
    assert document["grammar"]["forward_order"] == [
        "terminal_IR", "blaTEM", "tnpR", "res", "tnpA", "terminal_IR",
    ]
    assignment = document["classification"]["type_assignment"]
    assert assignment["primary_discriminator_roles"] == ["tnpR", "tnpA"]
    assert assignment["supporting_discriminator_roles"] == []
    assert assignment["unweighted_context_roles"] == ["res"]
    assert assignment["non_discriminator_roles"] == ["blaTEM"]
    assert assignment["discriminator_role_weights"] == {
        "tnpR": 5.0,
        "tnpA": 5.0,
    }


def test_reference_metadata_is_generated_from_definitions():
    metadata = tn123_reference_metadata()
    assert metadata["Tn2c_HM749967"]["name"] == "Tn2c"
    assert metadata["Tn2c_HM749967"]["type"] == "Tn2"
    assert metadata["Tn3_V00613"]["subtype"] == "V00613_legacy_9bp_duplication"


def test_declared_component_divergence_matches_canonical_sequences() -> None:
    document = load_tn123_definitions()
    references = Path(__file__).parents[1] / "matryoshka" / "references" / "tn1_tn2_tn3.fasta"
    records = {record.id: str(record.seq) for record in SeqIO.parse(references, "fasta")}
    sequences: dict[str, dict[str, str]] = {}
    for type_name, definition in document["types"].items():
        sequence = records[definition["canonical_reference"]]
        sequences[type_name] = {}
        for component in definition["components"]:
            role = str(component["role"])
            if role == "terminal_IR":
                continue
            start = int(component["start"])
            end = int(component["end"])
            if start < 0:
                start = len(sequence) + start + 1
            if end < 0:
                end = len(sequence) + end + 1
            sequences[type_name][role] = sequence[start - 1:end]

    aligner = PairwiseAligner(
        mode="global",
        match_score=2,
        mismatch_score=-1,
        open_gap_score=-5,
        extend_gap_score=-1,
    )
    evidence = document["classification"]["type_assignment"][
        "discriminator_evidence"
    ]["pairwise_identity_percent"]
    pairs = [("Tn1", "Tn2"), ("Tn1", "Tn3"), ("Tn2", "Tn3")]
    for role, expected_by_pair in evidence.items():
        for left, right in pairs:
            counts = aligner.align(sequences[left][role], sequences[right][role])[0].counts()
            denominator = counts.identities + counts.mismatches + counts.gaps
            observed = round(100 * counts.identities / denominator, 3)
            assert observed == expected_by_pair[f"{left}_vs_{right}"]


def test_markdown_export_contains_worked_expert_logic():
    report = definitions_as_markdown()
    assert "Required component order" in report
    assert "terminal_IR -> blaTEM -> tnpR -> res -> tnpA -> terminal_IR" in report
    assert "Tn2c" in report
    assert "9 bp duplication" in report
    assert "Related but different elements" in report
    assert "blaTEM is required cargo" in report
    assert "it contributes no weight" in report


def test_definitions_cli_exports_reviewable_yaml_json_and_markdown(tmp_path: Path):
    runner = CliRunner()
    markdown = runner.invoke(cli, ["definitions", "--format", "markdown"])
    assert markdown.exit_code == 0
    assert markdown.output.startswith("# Tn1/Tn2/Tn3 expert definitions")

    target = tmp_path / "definitions.json"
    result = runner.invoke(
        cli,
        ["definitions", "--format", "json", "--out", str(target)],
    )
    assert result.exit_code == 0
    assert json.loads(target.read_text())["definitions"]["Tn1_NC_008357"][
        "display_name"
    ] == "Tn1"

    schema = Path(__file__).parents[1] / "docs" / "schema" / "tn123-definitions-v1.schema.json"
    assert json.loads(schema.read_text())["properties"]["schema_version"]["const"] == "1.0"
