"""
hierarchy.py — Build containment hierarchy from a flat list of MGEFeatures.

A feature A is *strictly contained* in feature B when:
    B.start <= A.start  and  A.end <= B.end  and  (B.start, B.end) != (A.start, A.end)

Each feature is assigned to its smallest strict container; containers with
no parent become roots. Ties (two features with identical coordinates) are
broken deterministically by element_type + name so the output is stable.

Features with no meaningful coordinates (start == end == 0, e.g. mob_typer
replicons) are kept as roots but do not participate in containment.
"""

from __future__ import annotations

from .detect import MGEFeature


def _length(f: MGEFeature) -> int:
    return f.end - f.start


def _has_coords(f: MGEFeature) -> bool:
    return not (f.start == 0 and f.end == 0)


# Biologically-meaningful parent-child pairs. An IS element cannot
# *contain* an AMR cargo unless a higher-level feature (transposon /
# integron / genomic island) explicitly says so — otherwise an IS that
# happens to overlap a cargo gene spuriously becomes its parent.
_FORBIDDEN_NESTING: frozenset[tuple[str, str]] = frozenset({
    ("IS", "AMR"),
    ("IS", "cassette"),
    ("IS", "attC"),
    ("IS", "integrase"),
    ("IS", "integron"),
    ("IS", "transposon"),
    ("AMR", "AMR"),
    ("AMR", "IS"),
    ("res_site", "AMR"),
    ("res_site", "IS"),
})


def contains(parent: MGEFeature, child: MGEFeature) -> bool:
    """True if parent strictly contains child (proper superset span) AND
    the parent/child element types form a biologically meaningful pair.
    """
    if parent is child:
        return False
    if parent.start > child.start or parent.end < child.end:
        return False
    if (parent.start, parent.end) == (child.start, child.end):
        return False
    return (parent.element_type, child.element_type) not in _FORBIDDEN_NESTING


def build_hierarchy(features: list[MGEFeature]) -> list[MGEFeature]:
    """
    Assign each feature to its smallest strict-container parent.
    Returns the list of root-level features (those with no parent).

    Mutates features' `children` lists. Runs in O(n²) time — adequate for
    realistic feature counts (< 10⁴ per contig).
    """
    # Clear any pre-existing children so the function is idempotent
    for f in features:
        f.children = []

    coord_feats = [f for f in features if _has_coords(f)]
    no_coord = [f for f in features if not _has_coords(f)]

    # Sort ascending by length so smaller features search against larger candidates
    ordered = sorted(coord_feats, key=lambda f: (_length(f), f.start, f.element_type, f.name))

    parent_of: dict[int, MGEFeature] = {}
    for i, child in enumerate(ordered):
        best: MGEFeature | None = None
        best_len: int | None = None
        # Candidates are strictly longer features (later in ordered list)
        for cand in ordered[i + 1:]:
            if not contains(cand, child):
                continue
            cl = _length(cand)
            if best is None or cl < best_len:
                best = cand
                best_len = cl
        if best is not None:
            parent_of[id(child)] = best

    for child, parent in ((c, parent_of[id(c)]) for c in ordered if id(c) in parent_of):
        parent.children.append(child)

    # Keep children sorted by start for stable rendering
    for f in coord_feats:
        if f.children:
            f.children.sort(key=lambda c: (c.start, c.end, c.element_type))

    roots = [f for f in ordered if id(f) not in parent_of]
    roots.sort(key=lambda f: (f.start, f.end))
    return roots + no_coord
