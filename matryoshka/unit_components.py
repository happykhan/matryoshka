"""Component references and expert-rule assembly for reviewed unit transposons."""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from importlib.resources import files
from typing import Any

from Bio import SeqIO

from .detect import MGEFeature
from .unit_definitions import unit_rules


@dataclass(frozen=True)
class UnitComponentReference:
    reference_id: str
    parent_reference: str
    parent_family: str
    role: str
    element_type: str
    name: str
    strand: str
    category: str
    sequence: str


def _reference_records(filename: str) -> dict[str, str]:
    resource = files("matryoshka").joinpath("references", filename)
    with resource.open("r", encoding="utf-8") as stream:
        return {record.id: str(record.seq) for record in SeqIO.parse(stream, "fasta")}


def unit_component_references() -> list[UnitComponentReference]:
    """Extract independently searchable components from reviewed definitions."""
    references: list[UnitComponentReference] = []
    record_cache: dict[str, dict[str, str]] = {}
    for family, definition in unit_rules()["units"].items():
        reference = definition["canonical_reference"]
        filename = str(reference["file"])
        if filename not in record_cache:
            record_cache[filename] = _reference_records(filename)
        record_id = str(reference["record"])
        sequence = record_cache[filename][record_id]
        for index, component in enumerate(definition["components"], start=1):
            start, end = int(component["start"]), int(component["end"])
            references.append(UnitComponentReference(
                reference_id=f"{family}__{component['role']}__{index}__{start}_{end}",
                parent_reference=record_id,
                parent_family=str(family),
                role=str(component["role"]),
                element_type=str(component["element_type"]),
                name=str(component["name"]),
                strand=str(component["strand"]),
                category=str(component["category"]),
                sequence=sequence[start - 1:end],
            ))
    return references


def _profile_score(feature: MGEFeature, family: str) -> float:
    matches = feature.attributes.get("component_profile_matches", [])
    if not isinstance(matches, list):
        return 0.0
    return max(
        (
            float(match.get("profile_score", 0.0))
            for match in matches
            if isinstance(match, dict) and match.get("family") == family
        ),
        default=0.0,
    )


def _orientation(features: list[MGEFeature]) -> str:
    orientations = [
        str(feature.attributes.get("component_alignment_strand"))
        for feature in features
        if feature.attributes.get("component_alignment_strand") in {"+", "-"}
    ]
    return Counter(orientations).most_common(1)[0][0] if orientations else "+"


def _best_role_feature(
    features: list[MGEFeature],
    family: str,
    role: str,
) -> MGEFeature | None:
    candidates = [
        feature
        for feature in features
        if feature.attributes.get("component_role") == role
        and _profile_score(feature, family) > 0
    ]
    return max(
        candidates,
        key=lambda feature: (
            _profile_score(feature, family),
            float(feature.attributes.get("blast_coverage", 0)),
        ),
        default=None,
    )


def _candidate(
    family: str,
    definition: dict[str, Any],
    left: MGEFeature,
    right: MGEFeature,
    features: list[MGEFeature],
) -> tuple[MGEFeature, list[MGEFeature], float] | None:
    contained = [
        feature
        for feature in features
        if left.start <= feature.start and feature.end <= right.end
        and _profile_score(feature, family) > 0
    ]
    required_counts = {
        str(role): int(count)
        for role, count in definition["required_counts"].items()
    }
    terminal_irs = sorted(
        (feature for feature in contained if feature.attributes.get("component_role") == "terminal_IR"),
        key=lambda feature: (feature.start, feature.end),
    )
    if len(terminal_irs) != required_counts.get("terminal_IR", 0):
        return None

    selected = list(terminal_irs)
    for role, count in required_counts.items():
        if role == "terminal_IR":
            continue
        if count != 1:
            return None
        component = _best_role_feature(contained, family, role)
        if component is None:
            return None
        selected.append(component)
    selected.sort(key=lambda feature: (feature.start, feature.end))
    orientation = _orientation(selected)
    observed = [str(feature.attributes.get("component_role")) for feature in selected]
    expected = [str(role) for role in definition["required_order"]]
    if orientation == "-":
        expected.reverse()
    if observed != expected:
        return None

    minimum_score = float(
        unit_rules()["classification"]["minimum_component_profile_score_percent"]
    )
    profile_scores = [_profile_score(feature, family) for feature in selected]
    if min(profile_scores) < minimum_score:
        return None

    label = str(unit_rules()["classification"]["complete_label_template"]).format(
        family=family
    )
    parent = MGEFeature(
        element_type="transposon",
        family=family,
        name=label,
        start=left.start,
        end=right.end,
        strand=orientation,
        tsd_length=int(definition["tsd_length"]),
        attributes={
            "seqid": left.attributes.get("seqid", ""),
            "source": "component_assembly",
            "classification_basis": "expert_component_rules",
            "rule_based_family_call": family,
            "rule_based_type_call": family,
            "component_assembly_status": "complete",
            "component_order_valid": True,
            "detected_component_count": len(selected),
            "reference_projected_component_count": 0,
            "component_roles": observed,
            "component_requirements": {
                role: True for role in required_counts
            },
            "component_profile_score": round(sum(profile_scores) / len(profile_scores), 3),
            "expert_rule": definition["expert_rule"],
            "definition_version": unit_rules()["definition_version"],
            "variant_status": "rule_based_family_candidate",
            "naming_evidence": ["component_grammar", "component_profiles"],
        },
    )
    score = sum(profile_scores) - abs(
        (right.end - left.start + 1) - int(definition["canonical_reference"]["length"])
    ) / 1000
    return parent, selected, score


def _family_candidates(
    family: str,
    definition: dict[str, Any],
    features: list[MGEFeature],
) -> list[tuple[MGEFeature, list[MGEFeature], float]]:
    matching = [feature for feature in features if _profile_score(feature, family) > 0]
    irs = sorted(
        (feature for feature in matching if feature.attributes.get("component_role") == "terminal_IR"),
        key=lambda feature: (feature.start, feature.end),
    )
    anchors = [
        feature
        for feature in matching
        if feature.attributes.get("component_role") == definition["anchor_role"]
    ]
    minimum = int(definition["candidate_span_bp"]["minimum"])
    maximum = int(definition["candidate_span_bp"]["maximum"])
    candidates: list[tuple[MGEFeature, list[MGEFeature], float]] = []
    for left_index, left in enumerate(irs):
        for right in irs[left_index + 1:]:
            span = right.end - left.start + 1
            if not minimum <= span <= maximum:
                continue
            if not any(left.start <= anchor.start and anchor.end <= right.end for anchor in anchors):
                continue
            candidate = _candidate(family, definition, left, right, matching)
            if candidate is not None:
                candidates.append(candidate)
    return candidates


def _same_contig(left: MGEFeature, right: MGEFeature) -> bool:
    return left.attributes.get("seqid") == right.attributes.get("seqid")


def _contained(inner: MGEFeature, outer: MGEFeature) -> bool:
    return _same_contig(inner, outer) and outer.start <= inner.start and inner.end <= outer.end


def assemble_reviewed_unit_components(features: list[MGEFeature]) -> list[MGEFeature]:
    """Assemble Tn21/Tn1721/Tn1722 from independently detected components."""
    component_features = [
        feature
        for feature in features
        if feature.attributes.get("source") == "reviewed_unit_component_scan"
    ]
    by_contig: dict[str, list[MGEFeature]] = {}
    for feature in component_features:
        by_contig.setdefault(str(feature.attributes.get("seqid", "")), []).append(feature)

    units = unit_rules()["units"]
    assembled: list[MGEFeature] = []
    for contig_features in by_contig.values():
        candidates: list[tuple[MGEFeature, list[MGEFeature], float]] = []
        for family, definition in units.items():
            candidates.extend(_family_candidates(str(family), definition, contig_features))
        candidates.sort(
            key=lambda item: (
                len(units[item[0].family]["required_order"]),
                item[2],
            ),
            reverse=True,
        )
        selected: list[MGEFeature] = []
        for parent, _, _ in candidates:
            if parent.family == "Tn1722" and any(
                outer.family == "Tn1721" and _contained(parent, outer)
                for outer in selected
            ):
                continue
            if any(
                parent.family == outer.family and _contained(parent, outer)
                for outer in selected
            ):
                continue
            selected.append(parent)
        assembled.extend(selected)
    return assembled
