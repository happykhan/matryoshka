"""Review-report contract generated from executable expert definitions."""

from __future__ import annotations

from pathlib import Path

from click.testing import CliRunner

from matryoshka.__main__ import cli
from matryoshka.expert_report import expert_rules_as_html, expert_rules_as_markdown


def test_plain_language_report_covers_every_executable_rule_family() -> None:
    report = expert_rules_as_markdown()
    assert "Detect parts independently" in report
    assert "whole-element match is checked afterwards" in report
    assert "reference-supported names must not be mistaken" in report
    for family in ["Tn1", "Tn2", "Tn3", "Tn21", "Tn1721", "Tn1722"]:
        assert family in report
    assert "terminal_IR → blaTEM → tnpR → res → tnpA → terminal_IR" in report
    assert "There is no cross-gene weighting and no averaged type score" in report
    assert "Tn1/Tn3-like tnpR" in report
    assert "Tn1/Tn2/Tn3-group mosaic" in report
    assert "Primarily separates Tn2" in report
    assert "Primarily separates Tn3" in report
    assert "Questions for expert review" in report
    assert "matryoshka/tn123_definitions.yaml" in report


def test_html_report_is_standalone_printable_and_contains_tables() -> None:
    report = expert_rules_as_html()
    assert report.startswith("<!doctype html>")
    assert "<title>Expert rule review</title>" in report
    assert "@media print" in report
    assert "<table>" in report
    assert "Tn1721-like" in report


def test_expert_rules_cli_writes_both_review_formats(tmp_path: Path) -> None:
    runner = CliRunner()
    markdown = tmp_path / "rules.md"
    html = tmp_path / "rules.html"
    for fmt, target in [("markdown", markdown), ("html", html)]:
        result = runner.invoke(
            cli,
            ["expert-rules", "--format", fmt, "--out", str(target)],
        )
        assert result.exit_code == 0, result.output
        assert target.is_file()
    assert markdown.read_text().startswith("# Expert rule review")
    assert html.read_text().startswith("<!doctype html>")
