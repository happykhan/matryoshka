"""End-to-end smoke tests for the public CLI workflows."""

import json
from pathlib import Path

import pytest
from click.testing import CliRunner

from matryoshka.__main__ import cli

DATA = Path(__file__).parent / "test-data"
FASTA = DATA / "acceptance" / "pOXA-48a.fasta"
ISESCAN = DATA / "detector-output" / "pOXA-48a.isescan.tsv"
AMRFINDER = DATA / "detector-output" / "pOXA-48a.amrfinder.tsv"
TN123 = DATA / "reviewed-examples" / "Tn1-Tn2-Tn3.fasta"
ARBITRARY_TN123 = (
    Path(__file__).parents[1]
    / "demo-output"
    / "arbitrary-tn123"
    / "arbitrary-demo.fasta"
)

SKIP_BLAST = pytest.mark.skipif(
    __import__("shutil").which("blastn") is None,
    reason="blastn not on PATH",
)


class TestAnnotateCLI:
    def _run(self, fmt: str, tmp_path: Path) -> str:
        out = tmp_path / f"out.{fmt}"
        runner = CliRunner()
        result = runner.invoke(
            cli,
            [
                "annotate",
                str(FASTA),
                "--isescan", str(ISESCAN),
                "--amrfinder", str(AMRFINDER),
                "--format", fmt,
                "-o", str(out),
            ],
        )
        assert result.exit_code == 0, result.output
        assert out.exists()
        return out.read_text()

    def test_json(self, tmp_path):
        data = json.loads(self._run("json", tmp_path))
        assert isinstance(data, list)
        # Tn1999 transposon should be a root
        families = [f["family"] for f in data]
        assert "Tn1999" in families

    def test_gff3(self, tmp_path):
        text = self._run("gff3", tmp_path)
        assert text.startswith("##gff-version 3")
        assert "Tn1999" in text
        assert "blaOXA-48" in text
        assert "tsd=TGCTG" in text

    def test_genbank(self, tmp_path):
        text = self._run("genbank", tmp_path)
        assert "LOCUS" in text
        assert "mobile_element" in text
        assert "Tn1999" in text

    def test_wolvercote(self, tmp_path):
        text = self._run("wolvercote", tmp_path)
        assert "Tn1999" in text
        assert "blaOXA-48" in text

    def test_linear_svg(self, tmp_path):
        text = self._run("linear", tmp_path)
        assert text.startswith("<svg")
        assert "Tn1999" in text

    def test_missing_tool_output_errors(self, tmp_path):
        out = tmp_path / "out.json"
        runner = CliRunner()
        result = runner.invoke(
            cli,
            ["annotate", str(FASTA), "--no-reference-scan", "-o", str(out)],
        )
        assert result.exit_code != 0
        assert "detection tool output or enable the reference scan" in result.output.lower()

    def test_version(self):
        runner = CliRunner()
        result = runner.invoke(cli, ["--version"])
        assert result.exit_code == 0
        assert "matryoshka" in result.output


@SKIP_BLAST
def test_run_writes_versioned_result_directory(tmp_path):
    out = tmp_path / "result"
    result = CliRunner().invoke(
        cli,
        [
            "run", str(TN123), "--threads", "2", "--detectors", "none",
            "--out", str(out),
        ],
    )
    assert result.exit_code == 0, result.output
    document = json.loads((out / "annotation.json").read_text())
    summary = json.loads((out / "run.json").read_text())
    assert document["schema_version"] == "1.0"
    assert document["reference_database"]["profile"] == "validated"
    assert summary["locus_views"] == 3
    assert [locus["call"] for locus in summary["locus_outputs"]] == [
        "Tn1", "Tn2", "Tn3",
    ]
    for locus in summary["locus_outputs"]:
        assert (out / locus["locus_map"]).is_file()
        assert (out / locus["locus_table"]).is_file()
        assert (out / locus["hierarchy"]).is_file()
    assert summary["proof_status"] == "PASS"
    assert {item["status"] for item in summary["detectors"]} == {"disabled"}
    assert len(list((out / "locus-map").glob("*.svg"))) == 3
    assert len(list((out / "locus-table").glob("*.svg"))) == 3
    assert len(list((out / "hierarchy").glob("*.svg"))) == 3
    assert (out / "annotation.cell").is_file()
    assert (out / "proof" / "proof.json").is_file()
    assert (out / "proof" / "components.tsv").is_file()
    assert (out / "proof" / "matches.tsv").is_file()
    assert (out / "proof" / "report.html").is_file()
    assert (out / "expert-rules.md").is_file()
    assert (out / "expert-rules.html").is_file()
    assert summary["outputs"]["expert_rules_html"] == "expert-rules.html"
    assert "Matryoshka analysis" in result.output
    assert "Analysis complete" in result.output


@SKIP_BLAST
def test_run_quiet_suppresses_rich_progress(tmp_path: Path) -> None:
    out = tmp_path / "quiet-result"
    result = CliRunner().invoke(
        cli,
        [
            "run", str(TN123), "--threads", "2", "--detectors", "none",
            "--quiet", "--out", str(out),
        ],
    )
    assert result.exit_code == 0, result.output
    assert result.output == ""
    assert (out / "run.json").is_file()


@SKIP_BLAST
def test_arbitrary_sequence_proof_connects_detection_matching_and_outputs(tmp_path):
    out = tmp_path / "arbitrary-proof"
    result = CliRunner().invoke(
        cli,
        [
            "run", str(ARBITRARY_TN123), "--threads", "2",
            "--detectors", "none", "--out", str(out),
        ],
    )
    assert result.exit_code == 0, result.output

    proof = json.loads((out / "proof" / "proof.json").read_text())
    assert proof["summary"] == {
        "status": "PASS",
        "tn123_loci": 3,
        "passed": 3,
        "partial": 0,
        "failed": 0,
    }
    loci = proof["records"][0]["loci"]
    assert [locus["call"] for locus in loci] == ["Tn1-like", "Tn2-like", "Tn3"]
    assert [locus["known_element_match"]["best_match"] for locus in loci] == [
        "Tn1", "Tn2", "Tn3",
    ]
    assert [locus["strand"] for locus in loci] == ["+", "+", "-"]
    for locus in loci:
        assert locus["verdict"] == "PASS"
        assert all(locus["checks"].values())
        assert locus["assembly"]["detected_component_count"] == 6
        assert locus["assembly"]["reference_projected_component_count"] == 0
        assert {component["role"] for component in locus["components"]} == {
            "terminal_IR", "blaTEM", "tnpR", "res", "tnpA",
        }
        for relative_path in locus["outputs"].values():
            assert (out / relative_path).exists()

    tn2 = loci[1]
    tnpa = next(
        component for component in tn2["components"] if component["role"] == "tnpA"
    )
    assert tnpa["status"] == "complete_with_insertion"
    assert 790 <= tnpa["inserted_bases"] <= 810

    cell_format = (out / "annotation.cell").read_text()
    assert all(call in cell_format for call in ["Tn1-like", "Tn2-like", "Tn3"])
    hierarchy = (out / "hierarchy" / "arbitrary_demo_contig.svg").read_text()
    assert all(component in hierarchy for component in ["blaTEM", "tnpR", "res", "tnpA"])
    report = (out / "proof" / "report.html").read_text()
    assert "This report is generated from the annotation evidence" in report
    assert report.count('class="verdict pass"') == 3
    embedded = report.split('id="matryoshka-proof">', 1)[1].split("</script>", 1)[0]
    assert json.loads(embedded)["summary"]["status"] == "PASS"


def test_preflight_reports_core_and_optional_detectors() -> None:
    result = CliRunner().invoke(cli, ["preflight", "--format", "json"])
    assert result.exit_code == 0, result.output
    document = json.loads(result.output)
    assert document["core"]["ready"] is True
    assert [item["name"] for item in document["detectors"]] == [
        "amrfinder", "isescan", "integron",
    ]
