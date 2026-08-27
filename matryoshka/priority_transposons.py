"""Load and export the priority transposon coverage roadmap."""

from __future__ import annotations

import json
from importlib.resources import files
from typing import Any

import yaml

EXPECTED_TARGETS = frozenset({
    "Tn1", "Tn2", "Tn3", "Tn7", "Tn10", "Tn21", "Tn402", "Tn1331",
    "Tn1696", "Tn1721", "Tn1722", "Tn1999", "Tn2670", "Tn4401",
    "Tn5393", "Tn5403", "Tn6029",
})


def load_priority_transposons() -> dict[str, Any]:
    resource = files("matryoshka").joinpath("priority_transposons.yaml")
    document = yaml.safe_load(resource.read_text(encoding="utf-8"))
    if not isinstance(document, dict) or document.get("schema_version") != "1.0":
        raise ValueError("unsupported priority-transposon roadmap")
    targets = document.get("targets")
    if not isinstance(targets, list):
        raise ValueError("priority-transposon targets must be a list")
    names = [target.get("name") for target in targets if isinstance(target, dict)]
    if len(names) != len(targets) or set(names) != EXPECTED_TARGETS:
        raise ValueError("priority-transposon roadmap must contain the 17 reviewed targets")
    if len(names) != len(set(names)):
        raise ValueError("priority-transposon target names must be unique")
    allowed_status = set(document.get("status_levels", {}))
    for target in targets:
        required = {"name", "architecture", "status", "required_components", "remaining_work"}
        missing = required - set(target)
        if missing:
            raise ValueError(f"{target.get('name', 'target')} is missing: {', '.join(sorted(missing))}")
        if target["status"] not in allowed_status:
            raise ValueError(f"{target['name']} has an unknown status")
    return document


def priority_transposons_as_yaml() -> str:
    resource = files("matryoshka").joinpath("priority_transposons.yaml")
    text = resource.read_text(encoding="utf-8")
    return text if text.endswith("\n") else text + "\n"


def priority_transposons_as_json() -> str:
    return json.dumps(load_priority_transposons(), indent=2) + "\n"


def priority_transposons_as_markdown() -> str:
    document = load_priority_transposons()
    lines = [
        "# Priority transposon coverage roadmap",
        "",
        str(document["scope"]),
        "",
        "| Target | Architecture | Current status | Remaining work |",
        "| --- | --- | --- | --- |",
    ]
    for target in document["targets"]:
        lines.append(
            f"| {target['name']} | {target['architecture']} | {target['status']} | "
            f"{target['remaining_work']} |"
        )
    lines.extend(["", "## Acceptance contract", ""])
    lines.extend(f"- {case}" for case in document["acceptance_contract"]["required_cases"])
    lines.append("")
    return "\n".join(lines)
