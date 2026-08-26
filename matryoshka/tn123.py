"""Curated annotations for the canonical Tn1, Tn2 and Tn3 references.

The reference choices and feature layouts follow the material supplied by
Sally Partridge. Coordinates are 1-based and relative to each complete
transposon reference sequence.
"""

from __future__ import annotations

from dataclasses import dataclass

from .detect import MGEFeature


@dataclass(frozen=True)
class InternalFeature:
    element_type: str
    name: str
    start: int
    end: int
    strand: str


REFERENCE_METADATA: dict[str, dict[str, str]] = {
    "Tn1_NC_008357": {
        "element_type": "transposon",
        "family": "Tn1",
        "name": "Tn1",
        "source_accession": "NC_008357",
        "tn123_canonical": "true",
    },
    "Tn2_AY123253": {
        "element_type": "transposon",
        "family": "Tn2",
        "name": "Tn2",
        "source_accession": "AY123253",
        "tn123_canonical": "true",
    },
    "Tn3_HM749966": {
        "element_type": "transposon",
        "family": "Tn3",
        "name": "Tn3",
        "source_accession": "HM749966",
        "tn123_canonical": "true",
    },
}


_COMMON_LAYOUT = (
    InternalFeature("IR", "IRL", 1, 38, "."),
    InternalFeature("AMR", "{bla}", 148, 1008, "-"),
    InternalFeature("gene", "tnpR", 1191, 1748, "-"),
    InternalFeature("res_site", "res", 1754, 1867, "."),
    InternalFeature("gene", "tnpA", 1911, -34, "+"),
    InternalFeature("IR", "IRR", -38, -1, "."),
)

_BLA_NAMES = {
    "Tn1": "blaTEM-2",
    "Tn2": "blaTEM-1b",
    "Tn3": "blaTEM-1",
}


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
    if parent.family not in _BLA_NAMES:
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
    for item in _COMMON_LAYOUT:
        rel_start = _absolute_coordinate(item.start, reference_length)
        rel_end = _absolute_coordinate(item.end, reference_length)
        projected = _project_reference_interval(parent, rel_start, rel_end)
        if projected is None:
            continue
        start, end = projected
        name = _BLA_NAMES[parent.family] if item.name == "{bla}" else item.name
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
    parent.attributes["curated_internal_features"] = True
    parent.attributes["ir_length"] = 38
    parent.tsd_length = 5
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


def annotate_tn123(features: list[MGEFeature]) -> list[MGEFeature]:
    """Return curated sub-features for exact Tn1/Tn2/Tn3 calls."""
    out: list[MGEFeature] = []
    for feature in features:
        out.extend(curated_internal_features(feature))
        out.extend(inserted_sequence_features(feature))
    return out
