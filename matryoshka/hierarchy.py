"""
hierarchy.py — Build containment hierarchy from a flat list of MGEFeatures.

An element A is contained within element B if:
  B.start <= A.start and A.end <= B.end and A != B

Returns a tree (list of root-level elements with children populated).
"""

from .detect import MGEFeature


def build_hierarchy(features: list[MGEFeature]) -> list[MGEFeature]:
    """
    Assign each feature to its smallest containing parent.
    Returns only root-level features (not contained in anything else).
    """
    # Sort by length descending so larger elements are parents
    sorted_feats = sorted(features, key=lambda f: f.end - f.start, reverse=True)

    assigned = set()

    for i, parent in enumerate(sorted_feats):
        for j, child in enumerate(sorted_feats):
            if i == j or j in assigned:
                continue
            if parent.start <= child.start and child.end <= parent.end:
                # Check no smaller element already contains this child
                already_assigned = any(
                    sorted_feats[k].start <= child.start
                    and child.end <= sorted_feats[k].end
                    and (sorted_feats[k].end - sorted_feats[k].start)
                    < (parent.end - parent.start)
                    for k in range(len(sorted_feats))
                    if k != i and k != j
                )
                if not already_assigned:
                    parent.children.append(child)
                    assigned.add(j)

    roots = [f for i, f in enumerate(sorted_feats) if i not in assigned]
    return roots
