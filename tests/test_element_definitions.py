"""Contract tests for reviewable expert element definitions."""

from __future__ import annotations

import json
from pathlib import Path

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


def test_reference_metadata_is_generated_from_definitions():
    metadata = tn123_reference_metadata()
    assert metadata["Tn2c_HM749967"]["name"] == "Tn2c"
    assert metadata["Tn2c_HM749967"]["type"] == "Tn2"
    assert metadata["Tn3_V00613"]["subtype"] == "V00613_legacy_9bp_duplication"


def test_markdown_export_contains_worked_expert_logic():
    report = definitions_as_markdown()
    assert "Required component order" in report
    assert "terminal_IR -> blaTEM -> tnpR -> res -> tnpA -> terminal_IR" in report
    assert "Tn2c" in report
    assert "9 bp duplication" in report
    assert "Related but different elements" in report


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
