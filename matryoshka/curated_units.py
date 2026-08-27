"""Project YAML-defined feature maps into reviewed whole-unit calls."""

from __future__ import annotations

from dataclasses import dataclass

from .detect import MGEFeature
from .tn123 import _project_reference_interval, _project_strand, inserted_sequence_features
from .unit_definitions import load_unit_definitions


@dataclass(frozen=True)
class CuratedPart:
    name: str
    element_type: str
    family: str
    start: int
    end: int
    strand: str = "."
    fid: str = ""
    role: str = ""
    note: str = ""
    fragment: bool = False


def _unit_maps() -> dict[str, tuple[CuratedPart, ...]]:
    maps: dict[str, tuple[CuratedPart, ...]] = {}
    for map_name, parts in load_unit_definitions()["feature_maps"].items():
        maps[str(map_name)] = tuple(
            CuratedPart(
                name=str(part["name"]),
                element_type=str(part["element_type"]),
                family=str(part["family"]),
                start=int(part["start"]),
                end=int(part["end"]),
                strand=str(part.get("strand", ".")),
                fid=str(part["fid"]),
                role=str(part.get("role", "")),
                note=str(part.get("note", "")),
                fragment=bool(part.get("fragment", False)),
            )
            for part in parts
        )
    return maps


UNIT_MAPS = _unit_maps()


def _project(parent: MGEFeature, part: CuratedPart) -> MGEFeature | None:
    projected = _project_reference_interval(parent, part.start, part.end)
    if projected is None:
        return None
    start, end = projected
    strand = _project_strand(parent, part.strand)

    attrs: dict[str, object] = {
        "seqid": parent.attributes.get("seqid", ""),
        "source": "curated_reference",
        "evidence_class": "reference_projected",
        "fid": part.fid,
        "source_accession": parent.attributes.get("source_accession", ""),
        "reference_parent": parent.attributes.get("reference_id", parent.name),
    }
    if part.role:
        attrs["component_role"] = part.role
    if part.note:
        attrs["note"] = part.note
    if part.fragment:
        attrs["fragment"] = True
        attrs["type"] = "p"
    elif part.element_type == "IS":
        attrs["type"] = "c"
    return MGEFeature(
        element_type=part.element_type,
        family=part.family,
        name=part.name,
        start=start,
        end=end,
        strand=strand,
        attributes=attrs,
    )


def _detected_roles(parent: MGEFeature, features: list[MGEFeature]) -> set[str]:
    return {
        str(feature.attributes.get("component_role"))
        for feature in features
        if feature.attributes.get("source") == "reviewed_unit_component_scan"
        and feature.attributes.get("seqid") == parent.attributes.get("seqid")
        and parent.start <= feature.start
        and feature.end <= parent.end
        and feature.attributes.get("component_role")
    }


def annotate_curated_units(features: list[MGEFeature]) -> list[MGEFeature]:
    """Project only map parts not already supported by component detection."""
    children: list[MGEFeature] = []
    parents = [feature for feature in features if feature.element_type == "transposon"]
    definition_version = str(load_unit_definitions()["definition_version"])
    for parent in parents:
        reference_id = str(parent.attributes.get("reference_id", ""))
        map_key = reference_id if reference_id in UNIT_MAPS else parent.family
        parts = UNIT_MAPS.get(map_key)
        if not parts:
            continue
        coverage = float(
            parent.attributes.get("blast_subject_coverage")
            or parent.attributes.get("blast_coverage")
            or 0
        )
        if coverage < 95.0:
            continue
        detected_roles = _detected_roles(parent, features)
        projected = [
            part
            for part in parts
            if not part.role or part.role not in detected_roles
        ]
        parent.attributes["curated_internal_features"] = bool(projected)
        parent.attributes["independently_detected_internal_features"] = bool(
            detected_roles
        )
        parent.attributes["feature_db_version"] = definition_version
        for part in projected:
            child = _project(parent, part)
            if child is not None:
                children.append(child)
        children.extend(inserted_sequence_features(parent))
    return children
