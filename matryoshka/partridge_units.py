"""Curated internal maps for Sally Partridge's supplied MARA unit references."""

from __future__ import annotations

from dataclasses import dataclass

from .detect import MGEFeature
from .tn123 import _project_reference_interval, _project_strand, inserted_sequence_features


@dataclass(frozen=True)
class CuratedPart:
    name: str
    element_type: str
    family: str
    start: int
    end: int
    strand: str = "."
    fid: str = ""
    note: str = ""
    fragment: bool = False


UNIT_MAPS: dict[str, tuple[CuratedPart, ...]] = {
    "Tn21": (
        CuratedPart("IRR", "IR", "Tn21", 1, 38, fid="MARA-SP-Tn21-IRR"),
        CuratedPart("tnp", "transposon_component", "Tn21", 1, 4039,
                    fid="MARA-SP-Tn21-tnp"),
        CuratedPart("IRi", "IR", "Tn402", 4040, 4064, fid="MARA-SP-IRi"),
        CuratedPart("5'-CS", "integron_segment", "class1_integron", 4040, 5394,
                    fid="MARA-SP-5CS-GGG", note="AF071413 three-G variant"),
        CuratedPart("aadA1 cassette region", "cassette", "class1_integron", 5395, 6250,
                    fid="MARA-SP-In2-aadA1",
                    note="curated cassette region between supplied 5'-CS and 3'-CS"),
        CuratedPart("3'-CS", "integron_segment", "class1_integron", 6251, 8275,
                    fid="MARA-SP-3CS", note="2025 bp present of 2239 bp reference",
                    fragment=True),
        CuratedPart("tni", "transposon_component", "Tn402", 12358, 15039,
                    fid="MARA-SP-Tn402-tni",
                    note="2682 bp terminal fragment of the 4733 bp supplied tni region",
                    fragment=True),
        CuratedPart("IRt", "IR", "Tn402", 15015, 15039, strand="-",
                    fid="MARA-SP-IRt"),
        CuratedPart("mer", "transposon_component", "Tn21", 15040, 19672,
                    fid="MARA-SP-Tn21-mer"),
        CuratedPart("IRL", "IR", "Tn21", 19635, 19672,
                    fid="MARA-SP-Tn21-IRL"),
    ),
    "Tn1721_AB366441": (
        CuratedPart("IRL", "IR", "Tn1721", 1, 38, fid="MARA-SP-Tn1721-IRL"),
        CuratedPart("Tn1722-like backbone", "transposon_component", "Tn1721", 1, 5640,
                    fid="MARA-SP-Tn1722-backbone"),
        CuratedPart("IRR", "IR", "Tn1721", 5603, 5640, strand="-",
                    fid="MARA-SP-Tn1721-IRR-internal",
                    note="internal IRR copy at the partial duplication boundary"),
        CuratedPart("tet(A)-tnp duplicated region", "transposon_component", "Tn1721",
                    5641, 11128, fid="MARA-SP-Tn1721-tet"),
        CuratedPart("IRR", "IR", "Tn1721", 11091, 11128, strand="-",
                    fid="MARA-SP-Tn1721-IRR-terminal"),
    ),
    "Tn1722_AB366441": (
        CuratedPart("IRL", "IR", "Tn1722", 1, 38, fid="MARA-SP-Tn1721-IRL"),
        CuratedPart("IRR", "IR", "Tn1722", 5603, 5640, strand="-",
                    fid="MARA-SP-Tn1721-IRR"),
    ),
    "Tn4401_GU595196": (
        CuratedPart("IRL", "IR", "Tn4401", 1, 40, fid="MARA-SP-Tn4401-IRL"),
        CuratedPart("tnp", "transposon_component", "Tn4401", 1, 4927,
                    fid="MARA-SP-Tn4401-tnp"),
        CuratedPart("ISKpn7", "IS", "IS21", 4928, 6883, "+",
                    fid="MARA-SP-ISKpn7", note="99.387% to supplied GU595196 reference"),
        CuratedPart("blaKPC region", "transposon_component", "Tn4401", 6884, 8062,
                    fid="MARA-SP-Tn4401-KPC-region"),
        CuratedPart("ISKpn6", "IS", "IS1182", 8066, 9605, "-",
                    fid="MARA-SP-ISKpn6"),
        CuratedPart("right end", "transposon_component", "Tn4401", 9608, 9907,
                    fid="MARA-SP-Tn4401-right-end"),
        CuratedPart("IRR", "IR", "Tn4401", 9868, 9907, strand="-",
                    fid="MARA-SP-Tn4401-IRR"),
    ),
    "Tn5393_AF262622": (
        CuratedPart("IRL", "IR", "Tn5393", 1, 38, fid="MARA-SP-Tn5393-IRL"),
        CuratedPart("IRR", "IR", "Tn5393", 5433, 5470, strand="-",
                    fid="MARA-SP-Tn5393-IRR"),
    ),
    "Tn5403_X75779": (
        CuratedPart("IRL", "IR", "Tn5403", 1, 38, fid="MARA-SP-Tn5403-IR"),
        CuratedPart("IRR", "IR", "Tn5403", 3626, 3663, strand="-",
                    fid="MARA-SP-Tn5403-IR"),
    ),
}


def _project(parent: MGEFeature, part: CuratedPart) -> MGEFeature | None:
    projected = _project_reference_interval(parent, part.start, part.end)
    if projected is None:
        return None
    start, end = projected
    strand = _project_strand(parent, part.strand)

    attrs: dict[str, object] = {
        "seqid": parent.attributes.get("seqid", ""),
        "source": "curated_reference",
        "fid": part.fid,
        "source_accession": parent.attributes.get("source_accession", ""),
        "reference_parent": parent.attributes.get("reference_id", parent.name),
    }
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


def annotate_partridge_units(features: list[MGEFeature]) -> list[MGEFeature]:
    """Project curated components into complete Sally-reference unit calls."""
    children: list[MGEFeature] = []
    for parent in features:
        if parent.element_type != "transposon":
            continue
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
        parent.attributes["curated_internal_features"] = True
        parent.attributes["feature_db_version"] = "2026-08-26"
        for part in parts:
            projected = _project(parent, part)
            if projected is not None:
                children.append(projected)
        children.extend(inserted_sequence_features(parent))
    return children
