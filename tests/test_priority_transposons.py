from click.testing import CliRunner

from matryoshka.__main__ import cli
from matryoshka.priority_transposons import (
    EXPECTED_TARGETS,
    load_priority_transposons,
)


def test_priority_roadmap_contains_exactly_the_17_reviewed_targets() -> None:
    document = load_priority_transposons()
    assert {target["name"] for target in document["targets"]} == EXPECTED_TARGETS
    assert len(document["targets"]) == 17


def test_priority_roadmap_keeps_reference_only_coverage_explicit() -> None:
    document = load_priority_transposons()
    by_name = {target["name"]: target for target in document["targets"]}
    assert by_name["Tn1"]["status"] == "validated_component_rule"
    assert by_name["Tn7"]["status"] == "reference_supported"
    assert by_name["Tn6029"]["status"] == "definition_required"
    assert by_name["Tn2670"]["status"] == "definition_required"


def test_roadmap_is_exportable_for_review() -> None:
    result = CliRunner().invoke(cli, ["roadmap", "--format", "markdown"])
    assert result.exit_code == 0, result.output
    assert "Priority transposon coverage roadmap" in result.output
    assert "| Tn2670 |" in result.output
