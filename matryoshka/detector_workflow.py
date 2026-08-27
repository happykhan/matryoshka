"""Orchestrate optional third-party detectors for the public workflow."""

from __future__ import annotations

import subprocess
from collections.abc import Callable
from dataclasses import asdict, dataclass
from pathlib import Path

from .detect import (
    DetectorRuntime,
    DetectorVersion,
    detector_runtime,
    detector_version,
    run_amrfinder,
    run_integron_finder,
    run_isescan,
)

DETECTOR_ORDER = ("amrfinder", "isescan", "integron")
DETECTOR_LABELS = {
    "amrfinder": "AMRFinderPlus",
    "isescan": "ISEScan",
    "integron": "IntegronFinder",
}


@dataclass
class DetectorExecution:
    name: str
    label: str
    status: str
    source: str
    runtime: str
    output: str | None = None
    message: str | None = None
    tool_version: str | None = None
    database_version: str | None = None

    def document(self) -> dict[str, str | None]:
        return asdict(self)


class DetectorExecutionError(RuntimeError):
    """Raised when a required third-party detector cannot complete."""


Runner = Callable[[str | Path, str | Path, int], Path]
VersionProvider = Callable[[str], DetectorVersion]


def _runtime_label(info: DetectorRuntime) -> str:
    if info.source == "path":
        return info.executable
    if info.source == "pixi":
        return f"pixi:{info.platform}/{info.executable}"
    return info.reason


def _provenance_path(output: Path, result_root: Path) -> str:
    resolved = output.resolve()
    try:
        return str(resolved.relative_to(result_root.resolve()))
    except ValueError:
        return str(resolved)


def execute_detectors(
    fasta: str | Path,
    detector_dir: str | Path,
    *,
    mode: str,
    threads: int,
    supplied: dict[str, str | None],
    announce: Callable[[str], None] | None = None,
    runners: dict[str, Runner] | None = None,
    version_provider: VersionProvider = detector_version,
) -> tuple[dict[str, str | None], list[DetectorExecution]]:
    """Resolve supplied outputs and optionally run every supported detector.

    ``available`` is deliberately best-effort: a detector failure is retained in
    provenance and the component/reference workflow continues. ``all`` makes all
    three tools mandatory and fails the run on unavailable software or execution
    errors.
    """
    if mode not in {"none", "available", "all"}:
        raise ValueError(f"unknown detector mode: {mode}")
    output_paths = {name: supplied.get(name) for name in DETECTOR_ORDER}
    records: list[DetectorExecution] = []
    runner_map: dict[str, Runner] = runners or {
        "amrfinder": run_amrfinder,
        "isescan": run_isescan,
        "integron": run_integron_finder,
    }
    destination = Path(detector_dir)

    for name in DETECTOR_ORDER:
        label = DETECTOR_LABELS[name]
        provided = output_paths[name]
        if provided:
            records.append(DetectorExecution(
                name=name,
                label=label,
                status="provided",
                source="precomputed",
                runtime="not executed",
                output=str(Path(provided).resolve()),
            ))
            continue
        if mode == "none":
            records.append(DetectorExecution(
                name=name,
                label=label,
                status="disabled",
                source="none",
                runtime="not requested",
            ))
            continue

        info = detector_runtime(name)
        runtime = _runtime_label(info)
        if not info.available:
            record = DetectorExecution(
                name=name,
                label=label,
                status="unavailable",
                source=info.source,
                runtime=runtime,
                message=info.reason,
            )
            records.append(record)
            if mode == "all":
                raise DetectorExecutionError(
                    f"{label} is required by --detectors all but is unavailable: "
                    f"{record.message}"
                )
            if announce:
                announce(f"Detector unavailable; skipping {label}: {record.message}")
            continue

        if announce:
            announce(f"Running {label} ({runtime})…")
        try:
            output = runner_map[name](fasta, destination / name, threads)
        except (OSError, subprocess.SubprocessError) as error:
            record = DetectorExecution(
                name=name,
                label=label,
                status="failed",
                source=info.source,
                runtime=runtime,
                message=str(error),
            )
            records.append(record)
            if mode == "all":
                raise DetectorExecutionError(
                    f"{label} failed and is required by --detectors all: {error}"
                ) from error
            if announce:
                announce(f"Detector failed; continuing without {label}: {error}")
            continue

        resolved = str(Path(output).resolve())
        provenance_output = _provenance_path(Path(output), destination.parent)
        versions = version_provider(name)
        output_paths[name] = resolved
        records.append(DetectorExecution(
            name=name,
            label=label,
            status="completed",
            source=info.source,
            runtime=runtime,
            output=provenance_output,
            tool_version=versions.tool_version,
            database_version=versions.database_version,
        ))

    return output_paths, records


def preflight_runtimes() -> list[tuple[str, str, DetectorRuntime]]:
    """Return typed, non-mutating detector/platform readiness records."""
    return [
        (name, DETECTOR_LABELS[name], detector_runtime(name))
        for name in DETECTOR_ORDER
    ]
