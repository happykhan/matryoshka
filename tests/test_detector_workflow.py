from pathlib import Path

import pytest

from matryoshka import detector_workflow
from matryoshka.detect import DetectorRuntime, DetectorVersion


def _runners(tmp_path: Path, failing: str | None = None):
    def make(name: str):
        def run(_fasta, outdir, threads):
            assert threads == 3
            if name == failing:
                raise OSError("synthetic failure")
            output = Path(outdir) / f"{name}.tsv"
            output.parent.mkdir(parents=True, exist_ok=True)
            output.write_text("header\n")
            return output

        return run

    return {name: make(name) for name in detector_workflow.DETECTOR_ORDER}


def test_available_runs_managed_detectors_and_records_outputs(monkeypatch, tmp_path):
    monkeypatch.setattr(
        detector_workflow,
        "detector_runtime",
        lambda name: DetectorRuntime(
            available=True,
            source="pixi",
            executable=name,
            platform="linux-64",
            manifest="/project/pixi.toml",
            reason="test",
        ),
    )
    paths, records = detector_workflow.execute_detectors(
        tmp_path / "input.fasta",
        tmp_path / "detectors",
        mode="available",
        threads=3,
        supplied={},
        runners=_runners(tmp_path),
        version_provider=lambda _name: DetectorVersion("test", "db-test"),
    )
    assert all(paths.values())
    assert [record.status for record in records] == ["completed"] * 3
    assert records[0].output == "detectors/amrfinder/amrfinder.tsv"
    assert records[0].tool_version == "test"
    assert records[0].database_version == "db-test"


def test_available_retains_failure_but_all_makes_it_fatal(monkeypatch, tmp_path):
    monkeypatch.setattr(
        detector_workflow,
        "detector_runtime",
        lambda name: DetectorRuntime(
            available=True,
            source="path",
            executable=name,
            platform="linux-64",
            reason="test",
        ),
    )
    paths, records = detector_workflow.execute_detectors(
        tmp_path / "input.fasta",
        tmp_path / "detectors",
        mode="available",
        threads=3,
        supplied={},
        runners=_runners(tmp_path, failing="isescan"),
        version_provider=lambda _name: DetectorVersion(),
    )
    assert paths["isescan"] is None
    assert next(record for record in records if record.name == "isescan").status == "failed"

    with pytest.raises(detector_workflow.DetectorExecutionError, match="ISEScan failed"):
        detector_workflow.execute_detectors(
            tmp_path / "input.fasta",
            tmp_path / "required",
            mode="all",
            threads=3,
            supplied={},
            runners=_runners(tmp_path, failing="isescan"),
            version_provider=lambda _name: DetectorVersion(),
        )


def test_supplied_output_is_used_even_when_execution_is_disabled(tmp_path):
    supplied = tmp_path / "amr.tsv"
    supplied.write_text("header\n")
    paths, records = detector_workflow.execute_detectors(
        tmp_path / "input.fasta",
        tmp_path / "detectors",
        mode="none",
        threads=3,
        supplied={"amrfinder": str(supplied)},
        runners=_runners(tmp_path),
        version_provider=lambda _name: DetectorVersion(),
    )
    assert paths["amrfinder"] == str(supplied)
    assert records[0].status == "provided"
