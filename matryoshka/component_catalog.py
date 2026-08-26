"""Load and export the source-backed MARA component catalogue."""

from __future__ import annotations

import csv
import io
import json
from importlib.resources import files
from typing import Any

import yaml

REQUIRED_COMPONENT_FIELDS = frozenset({
    "level",
    "definition",
    "detection",
    "orientation",
    "containment",
    "render",
    "current_status",
})


def load_component_catalog() -> dict[str, Any]:
    """Return the bundled catalogue after validating its public contract."""
    resource = files("matryoshka").joinpath("mara_component_catalog.yaml")
    catalog = yaml.safe_load(resource.read_text(encoding="utf-8"))
    if not isinstance(catalog, dict):
        raise ValueError("MARA component catalogue must be a mapping")
    for top_level in ("schema_version", "component_types", "assembly_grammars"):
        if top_level not in catalog:
            raise ValueError(f"MARA component catalogue is missing {top_level}")
    for name, component in catalog["component_types"].items():
        missing = REQUIRED_COMPONENT_FIELDS - set(component)
        if missing:
            fields = ", ".join(sorted(missing))
            raise ValueError(f"component {name} is missing: {fields}")
        if "symbol" not in component["render"]:
            raise ValueError(f"component {name} has no render symbol")
    return catalog


def catalog_as_json() -> str:
    """Return the complete catalogue as stable, pretty JSON."""
    return json.dumps(load_component_catalog(), indent=2, sort_keys=False) + "\n"


def catalog_as_tsv() -> str:
    """Return the atomic component inventory as a spreadsheet-friendly table."""
    out = io.StringIO()
    fieldnames = (
        "component",
        "level",
        "definition",
        "detection",
        "orientation",
        "containment",
        "symbol",
        "status",
    )
    writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    for name, component in load_component_catalog()["component_types"].items():
        writer.writerow({
            "component": name,
            "level": component["level"],
            "definition": component["definition"],
            "detection": component["detection"],
            "orientation": component["orientation"],
            "containment": component["containment"],
            "symbol": component["render"]["symbol"],
            "status": component["current_status"],
        })
    return out.getvalue()
