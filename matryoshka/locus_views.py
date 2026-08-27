"""Extract readable viewports around validated mobile-element loci."""

from __future__ import annotations

from dataclasses import dataclass

from .detect import MGEFeature

TN123_FAMILIES = frozenset({"Tn1", "Tn2", "Tn3", "Tn3_family"})
PRIORITY_LOCUS_FAMILIES = TN123_FAMILIES | frozenset({
    "Tn7", "Tn10", "Tn21", "Tn402", "Tn1331", "Tn1696", "Tn1721",
    "Tn1722", "Tn1999", "Tn2670", "Tn4401", "Tn5393", "Tn5403", "Tn6029",
})
STRUCTURAL_LOCUS_TARGETS = frozenset({
    "transposition_unit",
    "composite_transposon",
    "IS26_translocatable_unit",
    "ISCR_capture_unit",
    "integron",
    "multiresistance_region",
    "accessory_region",
})


@dataclass(frozen=True)
class LocusView:
    target: MGEFeature
    view_start: int
    view_end: int
    roots: list[MGEFeature]

    @property
    def view_length(self) -> int:
        return self.view_end - self.view_start + 1

    @property
    def suffix(self) -> str:
        name = self.target.name.replace("/", "-").replace(" ", "_")
        return f"{name}__{self.target.start}-{self.target.end}"


def _walk(features: list[MGEFeature]) -> list[MGEFeature]:
    out: list[MGEFeature] = []
    for feature in features:
        out.append(feature)
        out.extend(_walk(feature.children))
    return out


def _is_locus_target(feature: MGEFeature) -> bool:
    if feature.element_type in STRUCTURAL_LOCUS_TARGETS:
        return True
    return (
        feature.element_type == "transposon"
        and (
            feature.family in PRIORITY_LOCUS_FAMILIES
            or feature.attributes.get("tn123_canonical") == "true"
        )
    )


def _clone_in_window(
    feature: MGEFeature,
    view_start: int,
    view_end: int,
) -> MGEFeature | None:
    if feature.end < view_start or feature.start > view_end:
        return None
    absolute_start = feature.start
    absolute_end = feature.end
    clipped_start = max(feature.start, view_start)
    clipped_end = min(feature.end, view_end)
    attrs = dict(feature.attributes)
    attrs["absolute_start"] = absolute_start
    attrs["absolute_end"] = absolute_end
    if absolute_start < view_start:
        attrs["view_clipped_left"] = True
    if absolute_end > view_end:
        attrs["view_clipped_right"] = True

    clone = MGEFeature(
        element_type=feature.element_type,
        family=feature.family,
        name=feature.name,
        start=clipped_start - view_start + 1,
        end=clipped_end - view_start + 1,
        strand=feature.strand,
        tsd_length=feature.tsd_length,
        tsd_seq=feature.tsd_seq,
        ir_left=feature.ir_left,
        ir_right=feature.ir_right,
        score=feature.score,
        attributes=attrs,
    )
    clone.children = [
        child_clone
        for child in feature.children
        if (child_clone := _clone_in_window(child, view_start, view_end)) is not None
    ]
    return clone


def extract_locus_views(
    roots: list[MGEFeature],
    seq_len: int,
    flank: int = 5_000,
) -> list[LocusView]:
    """Return one clipped hierarchy per validated locus-scale target locus."""
    candidates = [feature for feature in _walk(roots) if _is_locus_target(feature)]
    # An integron already visible inside a named transposon/TPU should not
    # create a second, near-identical locus. Broad resistance regions remain
    # additional context maps while their independently mobile children retain
    # their own readable locus maps.
    targets = sorted(
        (
            feature
            for feature in candidates
            if not (
                feature.element_type == "integron"
                and any(
                    other is not feature
                    and other.element_type
                    in {
                        "transposon",
                        "transposition_unit",
                        "composite_transposon",
                        "IS26_translocatable_unit",
                        "ISCR_capture_unit",
                    }
                    and other.start <= feature.start
                    and feature.end <= other.end
                    for other in candidates
                )
            )
        ),
        key=lambda feature: (feature.start, feature.end),
    )
    loci: list[LocusView] = []
    for target in targets:
        view_start = max(1, target.start - flank)
        view_end = min(seq_len, target.end + flank)
        clipped = [
            root_clone
            for root in roots
            if (root_clone := _clone_in_window(root, view_start, view_end)) is not None
        ]
        loci.append(LocusView(
            target=target,
            view_start=view_start,
            view_end=view_end,
            roots=clipped,
        ))
    return loci
