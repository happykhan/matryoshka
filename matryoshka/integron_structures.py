"""Assemble locus class-1-integron and Tn402 structures from exact components."""

from __future__ import annotations

from .detect import MGEFeature


def _seqid(feature: MGEFeature) -> str:
    return str(feature.attributes.get("seqid", ""))


def _same_record(left: MGEFeature, right: MGEFeature) -> bool:
    return _seqid(left) == _seqid(right)


def _distance_to_boundary(feature: MGEFeature, position: int) -> int:
    return min(abs(feature.start - position), abs(feature.end - position))


def infer_integron_structures(features: list[MGEFeature]) -> list[MGEFeature]:
    """Infer conservative parents from 5'-CS/3'-CS, IRi/IRt and tni evidence."""
    five_segments = [
        feature
        for feature in features
        if feature.element_type == "integron_segment" and feature.name == "5'-CS"
    ]
    three_segments = [
        feature
        for feature in features
        if feature.element_type == "integron_segment" and feature.name == "3'-CS"
    ]
    iris = [feature for feature in features if feature.element_type == "IR" and feature.name == "IRi"]
    irts = [feature for feature in features if feature.element_type == "IR" and feature.name == "IRt"]
    tni_regions = [
        feature
        for feature in features
        if feature.element_type == "transposon_component" and feature.family == "Tn402"
    ]
    existing_integrons = [feature for feature in features if feature.element_type == "integron"]

    emitted: list[MGEFeature] = []
    for five in five_segments:
        reverse = five.strand == "-"
        candidates = [
            three
            for three in three_segments
            if _same_record(five, three)
            and (
                three.end < five.start
                if reverse
                else three.start > five.end
            )
            and abs(three.start - five.end) <= 30_000
        ]
        if not candidates:
            continue
        three = (
            max(candidates, key=lambda feature: feature.end)
            if reverse
            else min(candidates, key=lambda feature: feature.start)
        )
        core_start = min(five.start, three.start)
        core_end = max(five.end, three.end)

        iri = min(
            (
                feature
                for feature in iris
                if _same_record(five, feature)
                and _distance_to_boundary(feature, five.start if not reverse else five.end) <= 50
            ),
            key=lambda feature: _distance_to_boundary(feature, five.start),
            default=None,
        )
        tni_candidates = [
            feature
            for feature in tni_regions
            if _same_record(five, feature)
            and (
                feature.end < three.start
                if reverse
                else feature.start > three.end
            )
            and min(abs(feature.start - three.end), abs(feature.end - three.start)) <= 15_000
        ]
        tni = (
            max(tni_candidates, key=lambda feature: feature.end)
            if reverse and tni_candidates
            else min(tni_candidates, key=lambda feature: feature.start, default=None)
        )
        distal_position = (
            tni.start if reverse and tni is not None
            else tni.end if tni is not None
            else three.start if reverse
            else three.end
        )
        irt = min(
            (
                feature
                for feature in irts
                if _same_record(five, feature)
                and _distance_to_boundary(feature, distal_position) <= 100
            ),
            key=lambda feature: _distance_to_boundary(feature, distal_position),
            default=None,
        )

        covered_by_existing = any(
            _same_record(five, parent)
            and parent.start <= core_start
            and core_end <= parent.end
            for parent in existing_integrons
        )
        if not covered_by_existing:
            emitted.append(MGEFeature(
                element_type="integron",
                family="class1_integron",
                name="class 1 integron candidate",
                start=core_start,
                end=core_end,
                strand="-" if reverse else "+",
                attributes={
                    "seqid": _seqid(five),
                    "source": "inference",
                    "structural_status": (
                        "conserved_segments_with_IR_ends"
                        if iri is not None and irt is not None
                        else "conserved_segment_pair"
                    ),
                    "iri_status": "confirmed" if iri is not None else "not_detected",
                    "irt_status": "confirmed" if irt is not None else "not_detected",
                    "note": "assembled from exact locus component references",
                },
            ))

        if (
            iri is not None
            and irt is not None
            and tni is not None
            and not tni.attributes.get("fragment")
        ):
            emitted.append(MGEFeature(
                element_type="transposon",
                family="Tn402",
                name="Tn402",
                start=min(iri.start, irt.start),
                end=max(iri.end, irt.end),
                strand="-" if reverse else "+",
                tsd_length=5,
                attributes={
                    "seqid": _seqid(five),
                    "source": "inference",
                    "structural_status": "complete_component_set",
                    "left_component": iri.name,
                    "right_component": irt.name,
                    "tni_component": tni.name,
                    "note": "IRi/IRt, conserved segments and complete tni region detected",
                },
            ))
    return emitted
