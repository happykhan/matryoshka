"""Curated annotations for the canonical Tn1, Tn2 and Tn3 references.

The reference choices and feature layouts follow the material supplied by
Sally Partridge. Coordinates are 1-based and relative to each complete
transposon reference sequence.
"""

from __future__ import annotations

from dataclasses import dataclass

from .detect import MGEFeature
from .element_definitions import (
    tn123_components,
    tn123_definition,
    tn123_reference_metadata,
    tn123_rules,
    tn123_type,
)


@dataclass(frozen=True)
class InternalFeature:
    element_type: str
    name: str
    start: int
    end: int
    strand: str


@dataclass(frozen=True)
class Tn123ComponentReference:
    """Independently searchable component extracted from a canonical locus."""

    reference_id: str
    parent_reference: str
    parent_family: str
    role: str
    element_type: str
    name: str
    start: int
    end: int
    strand: str
    sequence: str


REFERENCE_METADATA: dict[str, dict[str, str]] = tn123_reference_metadata()


def _internal_features(reference_id: str) -> tuple[InternalFeature, ...]:
    return tuple(
        InternalFeature(
            str(component["element_type"]),
            str(component["name"]),
            int(component["start"]),
            int(component["end"]),
            str(component["strand"]),
        )
        for component in tn123_components(reference_id)
    )


def component_references(
    reference_id: str,
    sequence: str,
) -> list[Tn123ComponentReference]:
    """Return canonical components as independent nucleotide references."""
    metadata = REFERENCE_METADATA[reference_id]
    definition = tn123_definition(reference_id)
    if not definition["component_reference"]:
        return []
    family = metadata["family"]
    references: list[Tn123ComponentReference] = []
    definition_components = tn123_components(reference_id)
    for item, component in zip(
        _internal_features(reference_id), definition_components, strict=True
    ):
        start = _absolute_coordinate(item.start, len(sequence))
        end = _absolute_coordinate(item.end, len(sequence))
        name = item.name
        role = str(component["role"])
        safe_name = "".join(character if character.isalnum() else "_" for character in name)
        references.append(Tn123ComponentReference(
            reference_id=f"{reference_id}__{safe_name}__{start}_{end}",
            parent_reference=reference_id,
            parent_family=family,
            role=role,
            element_type=item.element_type,
            name=name,
            start=start,
            end=end,
            strand=item.strand,
            sequence=sequence[start - 1:end],
        ))
    return references


def _absolute_coordinate(value: int, reference_length: int) -> int:
    return value if value > 0 else reference_length + value + 1


def _project(parent: MGEFeature, start: int, end: int) -> tuple[int, int]:
    if parent.strand == "-":
        return parent.end - end + 1, parent.end - start + 1
    return parent.start + start - 1, parent.start + end - 1


def _map_reference_position(parent: MGEFeature, position: int) -> int | None:
    segments = parent.attributes.get("reference_segments")
    if not isinstance(segments, list):
        return None
    for segment in segments:
        if not isinstance(segment, dict):
            continue
        try:
            qstart = int(segment["qstart"])
            qend = int(segment["qend"])
            sstart = int(segment["sstart"])
            send = int(segment["send"])
        except (KeyError, TypeError, ValueError):
            continue
        if not min(sstart, send) <= position <= max(sstart, send):
            continue
        if send == sstart:
            return qstart
        slope = (qend - qstart) / (send - sstart)
        return round(qstart + (position - sstart) * slope)
    return None


def _project_reference_interval(
    parent: MGEFeature,
    start: int,
    end: int,
) -> tuple[int, int] | None:
    mapped_start = _map_reference_position(parent, start)
    mapped_end = _map_reference_position(parent, end)
    if mapped_start is not None and mapped_end is not None:
        interval_start, interval_end = sorted((mapped_start, mapped_end))
        return interval_start, interval_end
    if parent.attributes.get("reference_segments"):
        return None
    return _project(parent, start, end)


def _interruption_size(parent: MGEFeature, start: int, end: int) -> int:
    segments = parent.attributes.get("reference_segments")
    if not isinstance(segments, list) or len(segments) < 2:
        return 0
    parsed: list[tuple[int, int, int, int]] = []
    reference_length = int(parent.attributes.get("reference_length", 0))
    for segment in segments:
        if not isinstance(segment, dict):
            continue
        try:
            qstart, qend = sorted((int(segment["qstart"]), int(segment["qend"])))
            sstart, send = sorted((int(segment["sstart"]), int(segment["send"])))
        except (KeyError, TypeError, ValueError):
            continue
        if parent.strand == "-":
            sstart, send = reference_length - send + 1, reference_length - sstart + 1
        parsed.append((qstart, qend, sstart, send))
    parsed.sort()
    largest = 0
    for left, right in zip(parsed, parsed[1:], strict=False):
        query_gap = max(0, right[0] - left[1] - 1)
        reference_gap = max(0, right[2] - left[3] - 1)
        inserted = max(0, query_gap - reference_gap)
        junction = left[3] if parent.strand != "-" else reference_length - left[3] + 1
        if start <= junction <= end:
            largest = max(largest, inserted)
    return largest


def _project_strand(parent: MGEFeature, strand: str) -> str:
    if strand == "." or parent.strand != "-":
        return strand
    return "+" if strand == "-" else "-"


def curated_internal_features(parent: MGEFeature) -> list[MGEFeature]:
    """Project canonical internal features onto an exact full-length call."""
    reference_id = str(parent.attributes.get("reference_id", ""))
    if reference_id in REFERENCE_METADATA:
        definition = tn123_definition(reference_id)
        layout = _internal_features(reference_id)
    elif parent.family in {"Tn1", "Tn2", "Tn3"}:
        type_definition = tn123_type(parent.family)
        reference_id = str(type_definition["canonical_reference"])
        definition = tn123_definition(reference_id)
        layout = _internal_features(reference_id)
    else:
        return []
    try:
        coverage = float(parent.attributes.get("blast_subject_coverage", 0.0))
        identity = float(parent.attributes.get("blast_identity", 0.0))
        reference_length = int(parent.attributes["reference_length"])
    except (KeyError, TypeError, ValueError):
        return []
    complete = bool(parent.attributes.get("left_end_covered", True)) and bool(
        parent.attributes.get("right_end_covered", True)
    )
    if coverage < 95.0 or identity < 95.0 or not complete:
        return []

    seqid = parent.attributes.get("seqid", "")
    source_accession = parent.attributes.get("source_accession", "")
    children: list[MGEFeature] = []
    for item in layout:
        rel_start = _absolute_coordinate(item.start, reference_length)
        rel_end = _absolute_coordinate(item.end, reference_length)
        projected = _project_reference_interval(parent, rel_start, rel_end)
        if projected is None:
            continue
        start, end = projected
        name = item.name
        interruption = _interruption_size(parent, rel_start, rel_end)
        attributes: dict[str, object] = {
            "seqid": seqid,
            "source": "curated_reference",
            "source_accession": source_accession,
            "parent_transposon": parent.name,
        }
        if interruption:
            attributes["interrupted"] = True
            attributes["inserted_bases"] = interruption
        children.append(MGEFeature(
            element_type=item.element_type,
            family=parent.family,
            name=name,
            start=start,
            end=end,
            strand=_project_strand(parent, item.strand),
            attributes=attributes,
        ))
    for item in definition.get("additional_features", []):
        rel_start = _absolute_coordinate(int(item["start"]), reference_length)
        rel_end = _absolute_coordinate(int(item["end"]), reference_length)
        projected = _project_reference_interval(parent, rel_start, rel_end)
        if projected is None:
            continue
        start, end = projected
        attributes: dict[str, object] = {
            "seqid": seqid,
            "source": "expert_definition",
            "source_accession": source_accession,
            "parent_transposon": parent.name,
            "evidence_class": item.get("evidence", "definition_projection"),
            "definition_id": reference_id,
        }
        if item.get("status"):
            attributes["structural_status"] = item["status"]
            if item["status"] == "partial":
                attributes["fragment"] = True
        children.append(MGEFeature(
            element_type=str(item["element_type"]),
            family=str(item["family"]),
            name=str(item["name"]),
            start=start,
            end=end,
            strand=_project_strand(parent, str(item["strand"])),
            attributes=attributes,
        ))
    parent.attributes["curated_internal_features"] = True
    rules = tn123_rules()["family"]
    parent.attributes["ir_length"] = int(rules["terminal_ir_length"])
    parent.tsd_length = int(rules["tsd_length"])
    return children


def inserted_sequence_features(parent: MGEFeature) -> list[MGEFeature]:
    """Emit query gaps between collinear reference HSPs as inserted sequence."""
    segments = parent.attributes.get("reference_segments")
    if not isinstance(segments, list) or len(segments) < 2:
        return []
    parsed: list[tuple[int, int, int, int]] = []
    reference_length = int(parent.attributes.get("reference_length", 0))
    for segment in segments:
        if not isinstance(segment, dict):
            continue
        try:
            qstart, qend = sorted((int(segment["qstart"]), int(segment["qend"])))
            sstart, send = sorted((int(segment["sstart"]), int(segment["send"])))
        except (KeyError, TypeError, ValueError):
            continue
        if parent.strand == "-":
            sstart, send = reference_length - send + 1, reference_length - sstart + 1
        parsed.append((qstart, qend, sstart, send))
    parsed.sort()

    out: list[MGEFeature] = []
    for index, (left, right) in enumerate(zip(parsed, parsed[1:], strict=False), start=1):
        query_gap = right[0] - left[1] - 1
        reference_gap = max(0, right[2] - left[3] - 1)
        inserted = query_gap - max(0, reference_gap)
        if inserted <= 20:
            continue
        start = left[1] + 1
        end = right[0] - 1
        out.append(MGEFeature(
            element_type="inserted_sequence",
            family="insertion",
            name=f"inserted sequence {index}",
            start=start,
            end=end,
            strand=".",
            attributes={
                "seqid": parent.attributes.get("seqid", ""),
                "source": "reference_hsp_gap",
                "parent_transposon": parent.name,
                "inserted_bases": inserted,
                "note": "sequence present between collinear Tn reference matches",
            },
        ))
    return out


def _required_component_counts() -> dict[str, int]:
    grammar = tn123_rules()["grammar"]
    return {
        str(component["role"]): int(component["minimum_count"])
        for component in grammar["required_components"]
    }


def _component_role(feature: MGEFeature) -> str:
    role = feature.attributes.get("component_role")
    if isinstance(role, str):
        return role
    if feature.element_type == "IR":
        return "terminal_IR"
    if feature.element_type == "AMR" and feature.name.startswith("blaTEM"):
        return "blaTEM"
    if feature.element_type == "res_site":
        return "res"
    return feature.name if feature.name in {"tnpA", "tnpR"} else ""


def _overlap_fraction(left: MGEFeature, right: MGEFeature) -> float:
    overlap = max(0, min(left.end, right.end) - max(left.start, right.start) + 1)
    shorter = min(left.end - left.start + 1, right.end - right.start + 1)
    return overlap / shorter if shorter else 0.0


def _components_for_parent(
    parent: MGEFeature,
    features: list[MGEFeature],
) -> list[MGEFeature]:
    return [
        feature
        for feature in features
        if feature.attributes.get("source") == "tn123_component_scan"
        and feature.attributes.get("seqid") == parent.attributes.get("seqid")
        and parent.start <= feature.start
        and feature.end <= parent.end
    ]


def _select_component_path(
    parent: MGEFeature,
    components: list[MGEFeature],
) -> tuple[list[MGEFeature], dict[str, bool], bool]:
    by_role: dict[str, list[MGEFeature]] = {}
    for component in components:
        role = _component_role(component)
        if role:
            by_role.setdefault(role, []).append(component)

    selected: list[MGEFeature] = []
    terminal_irs = sorted(by_role.get("terminal_IR", []), key=lambda item: item.start)
    if len(terminal_irs) >= 2:
        selected.extend([terminal_irs[0], terminal_irs[-1]])
    elif terminal_irs:
        selected.extend(terminal_irs)

    for role in ("blaTEM", "tnpR", "res", "tnpA"):
        candidates = by_role.get(role, [])
        if candidates:
            selected.append(max(
                candidates,
                key=lambda item: (
                    float(item.attributes.get("blast_coverage", 0)),
                    float(item.attributes.get("blast_identity", 0)),
                ),
            ))

    selected.sort(key=lambda item: (item.start, item.end))
    roles = [_component_role(feature) for feature in selected]
    grammar = tn123_rules()["grammar"]
    forward = [str(role) for role in grammar["forward_order"]]
    expected = forward if parent.strand != "-" else list(reversed(forward))
    required_counts = _required_component_counts()
    requirements = {
        role: len(by_role.get(role, [])) >= count
        for role, count in required_counts.items()
    }
    order_valid = roles == expected
    return selected, requirements, order_valid


def _record_component_assembly(
    parent: MGEFeature,
    components: list[MGEFeature],
) -> bool:
    selected, requirements, order_valid = _select_component_path(parent, components)
    component_complete = all(requirements.values()) and order_valid
    if selected:
        irs = [feature for feature in selected if _component_role(feature) == "terminal_IR"]
        if len(irs) == 2:
            if parent.strand == "-":
                irs[0].name, irs[1].name = "IRR", "IRL"
            else:
                irs[0].name, irs[1].name = "IRL", "IRR"
    parent.attributes.update({
        "component_assembly_status": (
            "complete" if component_complete else "partial_or_conflicting"
        ),
        "component_order_valid": order_valid,
        "component_requirements": requirements,
        "detected_component_count": len(selected),
        "component_evidence": [
            {
                "role": _component_role(feature),
                "name": feature.name,
                "start": feature.start,
                "end": feature.end,
                "strand": feature.strand,
                "evidence_class": feature.attributes.get("evidence_class"),
                "identity": feature.attributes.get("blast_identity"),
                "coverage": feature.attributes.get("blast_coverage"),
                "status": feature.attributes.get("structural_status"),
                "component_status": feature.attributes.get("component_status"),
                "reference": feature.attributes.get("component_reference"),
                "reference_segments": feature.attributes.get("reference_segments", []),
                "inserted_bases": feature.attributes.get("inserted_bases", 0),
            }
            for feature in selected
        ],
        "naming_evidence": [
            "component_grammar",
            *(
                ["whole_locus_alignment"]
                if parent.attributes.get("source") == "reference_scan"
                else []
            ),
        ],
    })
    if component_complete:
        parent.attributes["evidence_class"] = "sequence_detected_and_assembled"
    return component_complete


def assemble_tn123_components(features: list[MGEFeature]) -> list[MGEFeature]:
    """Validate reference parents and emit reference-independent candidates."""
    parents = [
        feature
        for feature in features
        if feature.element_type == "transposon"
        and (
            feature.family in {"Tn1", "Tn2", "Tn3", "Tn3_family"}
            or feature.attributes.get("tn123_canonical") == "true"
        )
    ]
    components = [
        feature
        for feature in features
        if feature.attributes.get("source") == "tn123_component_scan"
    ]
    for parent in parents:
        complete = _record_component_assembly(
            parent,
            _components_for_parent(parent, components),
        )
        if (
            not complete
            and parent.attributes.get("variant_status") == "exact_reference"
            and not parent.attributes.get("fragment")
        ):
            parent.attributes["reference_assigned_name"] = parent.name
            parent.name = f"{parent.name} reference-match candidate"
            parent.family = "Tn3_family"
            parent.attributes["note"] = (
                "whole-locus reference match lacks a complete independently "
                "detected component grammar"
            )

    emitted: list[MGEFeature] = []
    for tnpa in [feature for feature in components if _component_role(feature) == "tnpA"]:
        if any(_overlap_fraction(parent, tnpa) >= 0.8 for parent in parents):
            continue
        seqid = tnpa.attributes.get("seqid")
        nearby = [
            component
            for component in components
            if component.attributes.get("seqid") == seqid
            and component.start >= tnpa.start - int(
                tn123_rules()["grammar"]["candidate_component_window_bp"]
            )
            and component.end <= tnpa.end + int(
                tn123_rules()["grammar"]["candidate_component_window_bp"]
            )
        ]
        irs = sorted(
            (feature for feature in nearby if _component_role(feature) == "terminal_IR"),
            key=lambda item: item.start,
        )
        if len(irs) < 2:
            continue
        candidate = MGEFeature(
            element_type="transposon",
            family="Tn3_family",
            name="Tn3-family unit",
            start=irs[0].start,
            end=irs[-1].end,
            strand=tnpa.strand,
            tsd_length=5,
            attributes={
                "seqid": seqid,
                "source": "component_assembly",
                "best_match": "unresolved",
                "variant_status": "component_assembled",
            },
        )
        if _record_component_assembly(candidate, nearby):
            emitted.append(candidate)
    return emitted


def annotate_tn123(features: list[MGEFeature]) -> list[MGEFeature]:
    """Project only components that were not independently sequence-detected."""
    out: list[MGEFeature] = []
    parents = [feature for feature in features if feature.element_type == "transposon"]
    for parent in parents:
        projected = curated_internal_features(parent)
        detected = _components_for_parent(parent, features)
        detected_roles = {_component_role(feature) for feature in detected}
        retained = [
            feature
            for feature in projected
            if _component_role(feature) not in detected_roles
        ]
        for feature in retained:
            feature.attributes["evidence_class"] = "reference_projected"
        required_roles = set(_required_component_counts())
        projected_components = [
            feature for feature in retained
            if _component_role(feature) in required_roles
        ]
        projected_context = [
            feature for feature in retained
            if _component_role(feature) not in required_roles
        ]
        parent.attributes["curated_internal_features"] = bool(retained)
        parent.attributes["independently_detected_internal_features"] = bool(detected)
        # Only required grammar parts count as projected components. Reviewed
        # subtype context (for example the ISEcp1-associated insertion in
        # Tn2.1) is deliberately separate: it enriches the drawing but is not
        # used to manufacture the core Tn1/2/3 component proof.
        parent.attributes["reference_projected_component_count"] = len(
            projected_components
        )
        parent.attributes["reference_projected_context_feature_count"] = len(
            projected_context
        )
        out.extend(retained)
        out.extend(inserted_sequence_features(parent))
    return out
