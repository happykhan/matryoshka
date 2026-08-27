from pathlib import Path

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
