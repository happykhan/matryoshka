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
TN123 = DATA / "partridge-examples" / "Tn1-Tn2-Tn3.fasta"

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
        ["run", str(TN123), "--threads", "2", "--out", str(out)],
    )
    assert result.exit_code == 0, result.output
    document = json.loads((out / "annotation.json").read_text())
    summary = json.loads((out / "run.json").read_text())
    assert document["schema_version"] == "1.0"
    assert document["reference_database"]["profile"] == "validated"
    assert summary["mara_loci"] == 3
    assert len(list((out / "mara").glob("*.svg"))) == 3
    assert len(list((out / "mara-table").glob("*.svg"))) == 3
