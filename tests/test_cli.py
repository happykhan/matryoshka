"""End-to-end smoke tests for the `matryoshka annotate` CLI."""

import json
from pathlib import Path

import pytest
from click.testing import CliRunner

from matryoshka.__main__ import cli

DATA = Path(__file__).parent.parent / "data"
FASTA = DATA / "reference_plasmids/pOXA-48a.fasta"
ISESCAN = DATA / "isescan_test/reference_plasmids/pOXA-48a.fasta.tsv"
AMRFINDER = DATA / "amrfinder_test/pOXA-48a.tsv"

SKIP = pytest.mark.skipif(
    not (FASTA.exists() and ISESCAN.exists() and AMRFINDER.exists()),
    reason="pOXA-48a test data not present",
)


@SKIP
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
        result = runner.invoke(cli, ["annotate", str(FASTA), "-o", str(out)])
        assert result.exit_code != 0
        assert "at least one detection tool" in result.output.lower()

    def test_version(self):
        runner = CliRunner()
        result = runner.invoke(cli, ["--version"])
        assert result.exit_code == 0
        assert "matryoshka" in result.output
