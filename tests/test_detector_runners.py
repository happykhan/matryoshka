from pathlib import Path
from subprocess import CompletedProcess

import pytest

from matryoshka import detect


def test_detector_command_prefers_an_installed_executable(monkeypatch):
    monkeypatch.setattr(detect.shutil, "which", lambda name: f"/tools/{name}")
    assert detect._detector_command("amrfinder") == ["/tools/amrfinder"]


def test_detector_command_has_actionable_installed_package_error(monkeypatch):
    monkeypatch.setattr(detect.shutil, "which", lambda _name: None)
    monkeypatch.setattr(detect.platform, "system", lambda: "Darwin")
    with pytest.raises(FileNotFoundError, match="precomputed output"):
        detect._detector_command("integron")


def test_detector_available_does_not_require_a_sibling_checkout(monkeypatch, tmp_path):
    monkeypatch.setattr(detect.shutil, "which", lambda _name: None)
    monkeypatch.setattr(detect.platform, "system", lambda: "Darwin")
    monkeypatch.setattr(detect, "Path", lambda *_args: Path(tmp_path))
    assert not detect.detector_available("isescan")


def test_amrfinder_uses_the_bundled_osx64_environment_on_apple_silicon(
    monkeypatch, tmp_path,
):
    manifest = tmp_path / "pixi.toml"
    manifest.write_text("[workspace]\nname='test'\n")
    monkeypatch.setattr(
        detect.shutil,
        "which",
        lambda name: "/tools/pixi" if name == "pixi" else None,
    )
    monkeypatch.setattr(detect.platform, "system", lambda: "Darwin")
    monkeypatch.setattr(detect.platform, "machine", lambda: "arm64")
    monkeypatch.setattr(detect, "_project_manifest", lambda: manifest)

    command = detect._detector_command("amrfinder")
    assert command == [
        "/tools/pixi", "run", "--manifest-path", str(manifest),
        "--environment", "amrfinder", "--platform", "osx-64", "amrfinder",
    ]
    assert not detect.detector_available("isescan")


def test_amrfinder_database_is_initialised_only_when_missing(monkeypatch):
    calls = []
    responses = iter([
        CompletedProcess(["amrfinder", "-V"], 1, "", "No valid AMRFinder database"),
        CompletedProcess(["amrfinder", "-u"], 0, "", ""),
        CompletedProcess(["amrfinder", "-V"], 0, "Database version: test", ""),
    ])

    def fake_run(command, **kwargs):
        calls.append((command, kwargs))
        return next(responses)

    monkeypatch.setattr(detect.subprocess, "run", fake_run)
    detect._ensure_amrfinder_database(["amrfinder"])
    assert [command[-1] for command, _ in calls] == ["-V", "-u", "-V"]


def test_amrfinder_database_is_not_updated_when_valid(monkeypatch):
    calls = []

    def fake_run(command, **kwargs):
        calls.append(command)
        return CompletedProcess(command, 0, "Database version: test", "")

    monkeypatch.setattr(detect.subprocess, "run", fake_run)
    detect._ensure_amrfinder_database(["amrfinder"])
    assert calls == [["amrfinder", "-V"]]


def test_amrfinder_versions_are_parsed_for_run_provenance(monkeypatch):
    monkeypatch.setattr(detect, "_detector_command", lambda _name: ["amrfinder"])
    monkeypatch.setattr(
        detect.subprocess,
        "run",
        lambda *args, **kwargs: CompletedProcess(
            args[0],
            0,
            "Software version: 4.2.7\nDatabase version: 2026-08-07.1\n",
            "",
        ),
    )
    versions = detect.detector_version("amrfinder")
    assert versions.tool_version == "4.2.7"
    assert versions.database_version == "2026-08-07.1"
