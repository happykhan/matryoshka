"""Unit tests for build_hierarchy — the containment algorithm."""

from matryoshka.detect import MGEFeature
from matryoshka.hierarchy import build_hierarchy, contains


def mk(name, start, end, etype="IS", family="IS6") -> MGEFeature:
    return MGEFeature(
        element_type=etype, family=family, name=name,
        start=start, end=end, strand="+",
    )


class TestContains:
    def test_strict_containment(self):
        a = mk("outer", 100, 500)
        b = mk("inner", 200, 400)
        assert contains(a, b)
        assert not contains(b, a)

    def test_identical_coords_not_contained(self):
        a = mk("a", 100, 500)
        b = mk("b", 100, 500)
        assert not contains(a, b)
        assert not contains(b, a)

    def test_partial_overlap(self):
        a = mk("a", 100, 300)
        b = mk("b", 200, 400)
        assert not contains(a, b)
        assert not contains(b, a)


class TestBuildHierarchy:
    def test_empty(self):
        assert build_hierarchy([]) == []

    def test_single_feature(self):
        f = mk("lone", 100, 200)
        roots = build_hierarchy([f])
        assert roots == [f]
        assert f.children == []

    def test_single_level_nesting(self):
        parent = mk("parent", 100, 1000, etype="transposon")
        child = mk("child", 200, 300)
        roots = build_hierarchy([parent, child])
        assert roots == [parent]
        assert parent.children == [child]

    def test_double_nesting(self):
        outer = mk("outer", 100, 2000, etype="integron")
        mid = mk("mid", 200, 1500, etype="transposon")
        inner = mk("inner", 300, 500)
        roots = build_hierarchy([outer, mid, inner])
        assert roots == [outer]
        # Each feature should attach to its smallest container
        assert outer.children == [mid]
        assert mid.children == [inner]
        assert inner.children == []

    def test_smallest_container_wins(self):
        outer = mk("outer", 100, 2000)
        mid = mk("mid", 150, 1500)
        inner = mk("inner", 200, 300)
        roots = build_hierarchy([outer, mid, inner])
        assert roots == [outer]
        # Inner should be child of mid, not outer
        assert inner in mid.children
        assert inner not in outer.children

    def test_siblings_under_same_parent(self):
        parent = mk("parent", 100, 1000)
        a = mk("a", 150, 250)
        b = mk("b", 500, 600)
        roots = build_hierarchy([parent, a, b])
        assert roots == [parent]
        assert a in parent.children
        assert b in parent.children
        assert len(parent.children) == 2
        # Children sorted by start
        assert parent.children[0].start < parent.children[1].start

    def test_partial_overlap_becomes_roots(self):
        a = mk("a", 100, 500)
        b = mk("b", 400, 800)
        roots = build_hierarchy([a, b])
        assert a in roots and b in roots
        assert len(roots) == 2

    def test_identical_coord_duplicates(self):
        a = mk("a", 100, 500)
        b = mk("b", 100, 500)
        roots = build_hierarchy([a, b])
        # Neither strictly contains the other → both roots
        assert a in roots and b in roots
        assert len(roots) == 2

    def test_zero_coord_feature_kept_as_root(self):
        f = mk("replicon", 0, 0, etype="replicon")
        real = mk("real", 100, 500)
        roots = build_hierarchy([f, real])
        assert f in roots
        assert real in roots

    def test_idempotent(self):
        parent = mk("parent", 100, 1000)
        child = mk("child", 200, 300)
        build_hierarchy([parent, child])
        # Running again shouldn't duplicate children
        build_hierarchy([parent, child])
        assert parent.children == [child]
