"""Read and validate the bundled reference-database manifest."""

from __future__ import annotations

import hashlib
from importlib.resources import files
from pathlib import Path
from typing import Any

import yaml


def load_reference_manifest() -> dict[str, Any]:
    resource = files("matryoshka").joinpath("references", "manifest.yaml")
    manifest = yaml.safe_load(resource.read_text(encoding="utf-8"))
    if not isinstance(manifest, dict) or "database_version" not in manifest:
        raise ValueError("invalid bundled reference manifest")
    return manifest


def sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


REFERENCE_DATABASE_VERSION = str(load_reference_manifest()["database_version"])
