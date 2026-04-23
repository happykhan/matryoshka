"""
transposon.py — Infer composite transposons from IS element + cargo patterns.

Tn4401 (ISKpn7 + blaKPC + ISKpn6):
  - ISKpn7 = IS21 family, upstream of blaKPC
  - ISKpn6 = IS1182 family, downstream of blaKPC
  - Total span: ~10 kb; promoter gap (ISKpn7 end → blaKPC start) ≤ 200bp

IS26 composite transposons:
  - Two IS26 (IS6 family) copies flanking one or more cargo genes
  - Max span between IS pair: 20 kb (conservative)

Rules produce MGEFeature(element_type="transposon") objects that can be
fed into build_hierarchy() like any other feature.
"""

from __future__ import annotations

from .detect import MGEFeature

# Maximum allowed gap between IS element boundary and cargo gene start
_TN4401_PROMOTER_MAX = 500    # bp; ISKpn7 end → blaKPC start is ~56bp
_TN4401_DOWNSTREAM_MAX = 1000 # bp; blaKPC end → ISKpn6 start is ~344bp
_IS26_MAX_SPAN = 20_000       # bp; IS26 composite transposons rarely exceed this


def _overlaps_or_adjacent(a: MGEFeature, b: MGEFeature, gap: int) -> bool:
    """True if b starts within `gap` bp of a ending."""
    return b.start <= a.end + gap


def infer_tn4401(features: list[MGEFeature]) -> list[MGEFeature]:
    """
    Detect Tn4401: IS21 upstream of blaKPC within promoter distance,
    followed by IS1182 within downstream distance.

    Returns new MGEFeature(transposon, Tn4401) objects; does not modify input.
    """
    is21 = [f for f in features if f.family == "IS21"]
    is1182 = [f for f in features if f.family == "IS1182"]
    blaKPC = [f for f in features if f.element_type == "AMR" and "KPC" in f.name.upper()]

    transposons: list[MGEFeature] = []
    for kpc in blaKPC:
        # Find IS21 ending just before blaKPC
        upstream = [f for f in is21 if f.end <= kpc.start and kpc.start - f.end <= _TN4401_PROMOTER_MAX]
        if not upstream:
            continue
        # Pick closest
        up = max(upstream, key=lambda f: f.end)

        # Find IS1182 starting just after blaKPC
        downstream = [f for f in is1182 if f.start >= kpc.end and f.start - kpc.end <= _TN4401_DOWNSTREAM_MAX]
        if not downstream:
            continue
        down = min(downstream, key=lambda f: f.start)

        tn = MGEFeature(
            element_type="transposon",
            family="Tn4401",
            name="Tn4401",
            start=up.start,
            end=down.end,
            strand=kpc.strand,
        )
        transposons.append(tn)

    return transposons


def infer_is26_composites(features: list[MGEFeature]) -> list[MGEFeature]:
    """
    Detect IS26 composite transposons: pairs of IS6-family elements flanking
    one or more AMR genes, within _IS26_MAX_SPAN bp of each other.

    Returns new MGEFeature(transposon, IS26_composite) objects.
    """
    is26 = sorted(
        [f for f in features if f.family == "IS6"],
        key=lambda f: f.start,
    )
    amr = [f for f in features if f.element_type == "AMR"]

    composites: list[MGEFeature] = []
    seen: set[tuple[int, int]] = set()

    for i, left in enumerate(is26):
        for right in is26[i + 1:]:
            span = right.end - left.start
            if span > _IS26_MAX_SPAN:
                break

            # At least one AMR gene must sit between the two IS26s
            cargo = [
                a for a in amr
                if a.start >= left.end and a.end <= right.start
            ]
            if not cargo:
                continue

            key = (left.start, right.end)
            if key in seen:
                continue
            seen.add(key)

            gene_names = "+".join(c.name for c in cargo)
            tn = MGEFeature(
                element_type="transposon",
                family="IS26_composite",
                name=f"IS26::{gene_names}",
                start=left.start,
                end=right.end,
                strand=".",
            )
            composites.append(tn)

    return composites


def infer_transposons(features: list[MGEFeature]) -> list[MGEFeature]:
    """Run all transposon inference rules and return the new features."""
    return infer_tn4401(features) + infer_is26_composites(features)
