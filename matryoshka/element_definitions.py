"""Load, validate and export reviewable expert element definitions."""

from __future__ import annotations

import json
from importlib.resources import files
from typing import Any

import yaml

_REQUIRED_TOP_LEVEL = frozenset({
    "schema_version",
    "definition_version",
    "family",
    "grammar",
    "classification",
    "component_detection",
    "types",
    "definitions",
    "related_element_policy",
})

_REQUIRED_DEFINITION_FIELDS = frozenset({
    "display_name",
    "type",
    "subtype",
    "source_accession",
    "reference_kind",
    "component_reference",
    "review_status",
    "expert_rule",
})

_REQUIRED_COMPONENT_FIELDS = frozenset({
    "role",
    "name",
    "element_type",
    "start",
    "end",
    "strand",
})


def _validate_type_definition(
    type_name: str,
    type_definition: object,
    definitions: dict[str, Any],
) -> None:
    if not isinstance(type_definition, dict):
        raise ValueError(f"type {type_name} must be a mapping")
    components = type_definition.get("components")
    if not isinstance(components, list) or not components:
        raise ValueError(f"type {type_name} has no component definitions")
    for component in components:
        if not isinstance(component, dict):
            raise ValueError(f"type {type_name} has an invalid component")
        missing = _REQUIRED_COMPONENT_FIELDS - set(component)
        if missing:
            raise ValueError(
                f"type {type_name} component is missing: {', '.join(sorted(missing))}"
            )
    canonical = type_definition.get("canonical_reference")
    if canonical not in definitions:
        raise ValueError(f"type {type_name} canonical reference is undefined: {canonical}")


def _validate_named_definition(
    reference_id: str,
    definition: object,
    types: dict[str, Any],
) -> None:
    if not isinstance(definition, dict):
        raise ValueError(f"definition {reference_id} must be a mapping")
    missing = _REQUIRED_DEFINITION_FIELDS - set(definition)
    if missing:
        raise ValueError(f"definition {reference_id} is missing: {', '.join(sorted(missing))}")
    if definition["type"] not in types:
        raise ValueError(f"definition {reference_id} uses unknown type {definition['type']}")
    required_feature_fields = {"element_type", "family", "name", "start", "end", "strand"}
    for feature in definition.get("additional_features", []):
        if not isinstance(feature, dict):
            raise ValueError(f"definition {reference_id} has an invalid additional feature")
        feature_missing = required_feature_fields - set(feature)
        if feature_missing:
            raise ValueError(
                f"definition {reference_id} additional feature is missing: "
                f"{', '.join(sorted(feature_missing))}"
            )


def load_tn123_definitions() -> dict[str, Any]:
    """Return the bundled Tn1/Tn2/Tn3 definitions after contract validation."""
    resource = files("matryoshka").joinpath("tn123_definitions.yaml")
    document = yaml.safe_load(resource.read_text(encoding="utf-8"))
    if not isinstance(document, dict):
        raise ValueError("Tn1/Tn2/Tn3 definitions must be a mapping")
    missing = _REQUIRED_TOP_LEVEL - set(document)
    if missing:
        raise ValueError(f"Tn1/Tn2/Tn3 definitions are missing: {', '.join(sorted(missing))}")
    if document["schema_version"] != "1.0":
        raise ValueError("unsupported Tn1/Tn2/Tn3 definition schema")

    types = document["types"]
    definitions = document["definitions"]
    if not isinstance(types, dict) or not {"Tn1", "Tn2", "Tn3"} <= set(types):
        raise ValueError("Tn1/Tn2/Tn3 type definitions are incomplete")
    if not isinstance(definitions, dict):
        raise ValueError("Tn1/Tn2/Tn3 named definitions must be a mapping")

    for type_name, type_definition in types.items():
        _validate_type_definition(type_name, type_definition, definitions)

    for reference_id, definition in definitions.items():
        _validate_named_definition(reference_id, definition, types)
    return document


def tn123_definition(reference_id: str) -> dict[str, Any]:
    """Return one named whole-locus definition."""
    definitions = load_tn123_definitions()["definitions"]
    try:
        return definitions[reference_id]
    except KeyError as exc:
        raise KeyError(f"unknown Tn1/Tn2/Tn3 definition: {reference_id}") from exc


def tn123_type(type_name: str) -> dict[str, Any]:
    """Return one expert type definition."""
    types = load_tn123_definitions()["types"]
    try:
        return types[type_name]
    except KeyError as exc:
        raise KeyError(f"unknown Tn1/Tn2/Tn3 type: {type_name}") from exc


def tn123_reference_metadata() -> dict[str, dict[str, str]]:
    """Return FASTA-reference metadata derived from the expert definitions."""
    document = load_tn123_definitions()
    metadata: dict[str, dict[str, str]] = {}
    for reference_id, definition in document["definitions"].items():
        type_name = str(definition["type"])
        metadata[reference_id] = {
            "element_type": "transposon",
            "family": type_name,
            "name": str(definition["display_name"]),
            "type": type_name,
            "subtype": str(definition["subtype"]),
            "source_accession": str(definition["source_accession"]),
            "reference_kind": str(definition["reference_kind"]),
            "definition_id": reference_id,
            "definition_version": str(document["definition_version"]),
            "classification_group": "tn123",
            "tn123_canonical": "true",
        }
    return metadata


def tn123_components(reference_id: str) -> list[dict[str, Any]]:
    """Return expanded components for a reference, including subtype overrides."""
    definition = tn123_definition(reference_id)
    components = [dict(item) for item in tn123_type(str(definition["type"]))["components"]]
    bla_allele = definition.get("bla_allele")
    if bla_allele:
        for component in components:
            if component["role"] == "blaTEM":
                component["name"] = str(bla_allele)
    return components


def tn123_rules() -> dict[str, Any]:
    """Return the grammar, classification and detection rules used by the engine."""
    document = load_tn123_definitions()
    return {
        "family": document["family"],
        "grammar": document["grammar"],
        "classification": document["classification"],
        "component_detection": document["component_detection"],
        "related_element_policy": document["related_element_policy"],
    }


def definitions_as_json() -> str:
    """Return the exact expert-definition document as stable JSON."""
    return json.dumps(load_tn123_definitions(), indent=2, sort_keys=False) + "\n"


def definitions_as_yaml() -> str:
    """Return the exact human-maintained YAML used by the detector."""
    resource = files("matryoshka").joinpath("tn123_definitions.yaml")
    text = resource.read_text(encoding="utf-8")
    return text if text.endswith("\n") else text + "\n"


def _lines(values: list[str]) -> str:
    return "\n".join(f"- {value}" for value in values)


def definitions_as_markdown() -> str:
    """Render a review document from the same rules used by the detector."""
    document = load_tn123_definitions()
    family = document["family"]
    grammar = document["grammar"]
    classification = document["classification"]
    sections = [
        f"# {document['name']}",
        "",
        f"Definition version: `{document['definition_version']}`",
        "",
        str(document["purpose"]),
        "",
        "## Family rule",
        "",
        str(family["expert_description"]),
        "",
        "Required component order:",
        "",
        f"`{' -> '.join(grammar['forward_order'])}`",
        "",
        str(grammar["expert_rule"]),
        "",
        "## Classification rules",
        "",
        f"Type assignment: {classification['type_assignment']['expert_rule']}",
        "",
        f"Secondary reference comparison: {classification['reference_comparison']['expert_rule']}",
        "",
        f"Exact definition match: {classification['exact_definition_match']['expert_rule']}",
        "",
        f"Close variant: {classification['close_variant']['expert_rule']}",
        "",
        f"Fragment: {classification['fragment']['expert_rule']}",
    ]
    for type_name, type_definition in document["types"].items():
        sections.extend([
            "",
            f"## Type: {type_name}",
            "",
            str(type_definition["expert_description"]),
            "",
            "| Order | Role | Name | Coordinates | Strand |",
            "| ---: | --- | --- | --- | :---: |",
        ])
        for index, component in enumerate(type_definition["components"], start=1):
            sections.append(
                f"| {index} | {component['role']} | {component['name']} | "
                f"{component['start']}..{component['end']} | {component['strand']} |"
            )
    sections.extend(["", "## Named definitions and subtypes", ""])
    for reference_id, definition in document["definitions"].items():
        sections.extend([
            f"### {definition['display_name']} (`{reference_id}`)",
            "",
            f"Type: `{definition['type']}`. Subtype: `{definition['subtype']}`. "
            f"Source accession: `{definition['source_accession']}`.",
            "",
            str(definition["expert_rule"]),
            "",
            f"Review status: {definition['review_status']}.",
        ])
        differences = definition.get("known_differences_from_parent", [])
        if differences:
            sections.extend(["", "Known differences from the parent definition:", ""])
            sections.extend(f"- {difference}" for difference in differences)
        sections.append("")
    sections.extend([
        "## Related but different elements",
        "",
        str(document["related_element_policy"]["expert_rule"]),
        "",
        _lines(document["related_element_policy"]["examples_expected_not_to_be_named_tn123"]),
        "",
    ])
    return "\n".join(sections)
