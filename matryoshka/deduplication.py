"""Resolve overlapping rule, reference and curated component calls."""

from __future__ import annotations

from dataclasses import dataclass

from .detect import MGEFeature


def _overlap(left: MGEFeature, right: MGEFeature) -> bool:
    if left.end < right.start or right.end < left.start:
        return False
    overlap = min(left.end, right.end) - max(left.start, right.start)
    shorter = min(left.end - left.start, right.end - right.start)
    return shorter > 0 and overlap / shorter > 0.5


def _same_reference_unit(left: MGEFeature, right: MGEFeature) -> bool:
    if left.attributes.get("seqid") != right.attributes.get("seqid"):
        return False
    if not _overlap(left, right):
        return False
    return (
        left.family == right.family
        or left.name == right.name
        or (left.family == "Tn3" and left.name in {right.family, right.name})
        or (right.family == "Tn3" and right.name in {left.family, left.name})
    )


def _reference_rank(feature: MGEFeature) -> tuple[float, float, float, int]:
    provenance = 1.0 if feature.attributes.get("provenance") == "expert_reviewed" else 0.0
    coverage = float(
        feature.attributes.get("blast_subject_coverage")
        or feature.attributes.get("blast_coverage")
        or 0
    )
    identity = float(feature.attributes.get("blast_identity") or 0)
    return provenance, coverage, identity, feature.end - feature.start


def _selected_reference_ids(reference_units: list[MGEFeature]) -> set[int]:
    selected: list[MGEFeature] = []
    for candidate in sorted(reference_units, key=_reference_rank, reverse=True):
        if any(_same_reference_unit(candidate, existing) for existing in selected):
            continue
        selected.append(candidate)
    return {id(feature) for feature in selected}


@dataclass(frozen=True)
class _FeatureIndex:
    selected_reference_ids: set[int]
    blast_by_family: dict[str, list[MGEFeature]]
    exact_reference_is: list[MGEFeature]
    curated_components: list[MGEFeature]


def _index_features(features: list[MGEFeature]) -> _FeatureIndex:
    reference_units = [
        feature
        for feature in features
        if feature.element_type == "transposon"
        and feature.attributes.get("source") == "reference_scan"
    ]
    blast_by_family: dict[str, list[MGEFeature]] = {}
    for feature in reference_units:
        blast_by_family.setdefault(feature.family, []).append(feature)
    return _FeatureIndex(
        selected_reference_ids=_selected_reference_ids(reference_units),
        blast_by_family=blast_by_family,
        exact_reference_is=[
            feature
            for feature in features
            if feature.element_type == "IS"
            and feature.attributes.get("source") == "reference_scan"
        ],
        curated_components=[
            feature
            for feature in features
            if feature.attributes.get("source") == "curated_reference"
        ],
    )


def _duplicates_curated_component(
    feature: MGEFeature,
    curated_components: list[MGEFeature],
) -> bool:
    if feature.attributes.get("source") != "reference_scan":
        return False
    feature_fid = feature.attributes.get("fid")
    return any(
        feature.element_type == curated.element_type
        and (
            feature.name == curated.name
            or (feature_fid and feature_fid == curated.attributes.get("fid"))
        )
        and _overlap(feature, curated)
        for curated in curated_components
    )


def _should_drop(feature: MGEFeature, index: _FeatureIndex) -> bool:
    source = feature.attributes.get("source")
    if (
        feature.element_type == "transposon"
        and source == "reference_scan"
        and id(feature) not in index.selected_reference_ids
    ):
        return True
    if _duplicates_curated_component(feature, index.curated_components):
        return True
    if (
        feature.element_type == "IS"
        and source not in {"reference_scan", "curated_reference"}
        and any(_overlap(feature, exact) for exact in index.exact_reference_is)
    ):
        return True
    return (
        feature.element_type == "transposon"
        and source != "reference_scan"
        and any(
            _overlap(feature, reference)
            for reference in index.blast_by_family.get(feature.family, [])
        )
    )


def suppress_redundant_inference(features: list[MGEFeature]) -> list[MGEFeature]:
    """Prefer reviewed/reference evidence over overlapping duplicate calls."""
    index = _index_features(features)
    return [feature for feature in features if not _should_drop(feature, index)]
