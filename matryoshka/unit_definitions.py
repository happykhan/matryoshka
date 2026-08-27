"""Load component rules and feature maps for reviewed unit transposons."""

from __future__ import annotations

import json
from importlib.resources import files
from typing import Any

import yaml

_REQUIRED_TOP_LEVEL = frozenset({
    "schema_version",
    "definition_version",
    "component_detection",
    "classification",
    "units",
    "feature_maps",
})
_REQUIRED_UNIT_FIELDS = frozenset({
    "display_name",
    "canonical_reference",
    "terminal_ir_length",
    "tsd_length",
    "candidate_span_bp",
    "anchor_role",
    "required_order",
    "required_counts",
    "expert_rule",
    "components",
})
_REQUIRED_COMPONENT_FIELDS = frozenset({
    "role",
    "name",
    "element_type",
    "start",
    "end",
    "strand",
    "category",
})


def load_unit_definitions() -> dict[str, Any]:
    """Return the reviewed unit rules after validating their public contract."""
    resource = files("matryoshka").joinpath("unit_definitions.yaml")
    document = yaml.safe_load(resource.read_text(encoding="utf-8"))
    if not isinstance(document, dict):
        raise ValueError("unit definitions must be a mapping")
    missing = _REQUIRED_TOP_LEVEL - set(document)
    if missing:
        raise ValueError(f"unit definitions are missing: {', '.join(sorted(missing))}")
    if document["schema_version"] != "1.0":
        raise ValueError("unsupported unit-definition schema")

    units = document["units"]
    if not isinstance(units, dict) or not {"Tn21", "Tn1721", "Tn1722"} <= set(units):
        raise ValueError("Tn21/Tn1721/Tn1722 definitions are incomplete")
    for family, definition in units.items():
        if not isinstance(definition, dict):
            raise ValueError(f"unit {family} must be a mapping")
        fields = _REQUIRED_UNIT_FIELDS - set(definition)
        if fields:
            raise ValueError(f"unit {family} is missing: {', '.join(sorted(fields))}")
        reference = definition["canonical_reference"]
        if not isinstance(reference, dict) or not {"file", "record", "accession"} <= set(reference):
            raise ValueError(f"unit {family} has an invalid canonical reference")
        components = definition["components"]
        if not isinstance(components, list) or not components:
            raise ValueError(f"unit {family} has no component definitions")
        roles: list[str] = []
        for component in components:
            if not isinstance(component, dict):
                raise ValueError(f"unit {family} has an invalid component")
            component_fields = _REQUIRED_COMPONENT_FIELDS - set(component)
            if component_fields:
                raise ValueError(
                    f"unit {family} component is missing: "
                    f"{', '.join(sorted(component_fields))}"
                )
            if int(component["start"]) > int(component["end"]):
                raise ValueError(f"unit {family} component has reversed coordinates")
            roles.append(str(component["role"]))
        required_counts = {
            str(role): int(count)
            for role, count in definition["required_counts"].items()
        }
        for role, count in required_counts.items():
            if roles.count(role) < count:
                raise ValueError(f"unit {family} lacks {count} definitions for role {role}")
        if len(definition["required_order"]) != sum(required_counts.values()):
            raise ValueError(f"unit {family} order and required counts disagree")

    feature_maps = document["feature_maps"]
    if not isinstance(feature_maps, dict):
        raise ValueError("unit feature maps must be a mapping")
    for map_name, parts in feature_maps.items():
        if not isinstance(parts, list):
            raise ValueError(f"feature map {map_name} must be a list")
        for part in parts:
            if not isinstance(part, dict):
                raise ValueError(f"feature map {map_name} has an invalid part")
            required = {"name", "element_type", "family", "start", "end", "fid"}
            missing_map_fields = required - set(part)
            if missing_map_fields:
                raise ValueError(
                    f"feature map {map_name} part is missing: "
                    f"{', '.join(sorted(missing_map_fields))}"
                )
    return document


def unit_rules() -> dict[str, Any]:
    """Return the component-detection and classification rules."""
    document = load_unit_definitions()
    return {
        "definition_version": document["definition_version"],
        "component_detection": document["component_detection"],
        "classification": document["classification"],
        "units": document["units"],
    }


def unit_definitions_as_yaml() -> str:
    resource = files("matryoshka").joinpath("unit_definitions.yaml")
    text = resource.read_text(encoding="utf-8")
    return text if text.endswith("\n") else text + "\n"


def unit_definitions_as_json() -> str:
    return json.dumps(load_unit_definitions(), indent=2, sort_keys=False) + "\n"


def unit_definitions_as_markdown() -> str:
    document = load_unit_definitions()
    sections = [
        f"# {document['name']}",
        "",
        f"Definition version: `{document['definition_version']}`",
        "",
        str(document["purpose"]),
    ]
    for family, definition in document["units"].items():
        reference = definition["canonical_reference"]
        sections.extend([
            "",
            f"## {family}",
            "",
            str(definition["expert_rule"]),
            "",
            f"Reviewed reference: `{reference['record']}` from accession "
            f"`{reference['accession']}`.",
            "",
            f"Required order: `{' -> '.join(definition['required_order'])}`",
            "",
            "| Role | Name | Type | Coordinates | Strand | Required |",
            "| --- | --- | --- | --- | :---: | :---: |",
        ])
        required = set(definition["required_counts"])
        for component in definition["components"]:
            sections.append(
                f"| {component['role']} | {component['name']} | "
                f"{component['element_type']} | {component['start']}..{component['end']} | "
                f"{component['strand']} | {'yes' if component['role'] in required else 'no'} |"
            )
    sections.append("")
    return "\n".join(sections)
