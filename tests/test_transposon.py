"""Tests for transposon inference against pKpQIL data."""

from pathlib import Path

from matryoshka.detect import parse_amrfinder, parse_isescan
from matryoshka.transposon import infer_is26_composites, infer_tn4401

DATA = Path(__file__).parent / "test-data" / "detector-output"
ISESCAN_TSV = DATA / "pKpQIL.isescan.tsv"
AMRFINDER_TSV = DATA / "pKpQIL.amrfinder.tsv"


def pKpQIL_features():
    return (
        parse_isescan(ISESCAN_TSV) +
        parse_amrfinder(AMRFINDER_TSV)
    )


class TestInferTn4401:
    def test_one_tn4401_found(self):
        feats = pKpQIL_features()
        transposons = infer_tn4401(feats)
        tn4401 = [t for t in transposons if t.family == "Tn4401"]
        assert len(tn4401) == 1

    def test_tn4401_span(self):
        feats = pKpQIL_features()
        tn = infer_tn4401(feats)[0]
        # ISKpn7 start (13998) to ISKpn6 end (20019)
        assert tn.start == 13998
        assert tn.end == 20019

    def test_tn4401_contains_blakpc(self):
        feats = pKpQIL_features()
        from matryoshka.hierarchy import build_hierarchy
        tn = infer_tn4401(feats)
        roots = build_hierarchy(feats + tn)
        tn_node = next(r for r in roots if r.family == "Tn4401")
        child_names = [c.name for c in tn_node.children]
        assert "blaKPC-3" in child_names

    def test_tn4401_contains_iskpn7_and_iskpn6(self):
        feats = pKpQIL_features()
        from matryoshka.hierarchy import build_hierarchy
        tn = infer_tn4401(feats)
        roots = build_hierarchy(feats + tn)
        tn_node = next(r for r in roots if r.family == "Tn4401")
        families = {c.family for c in tn_node.children}
        assert "IS21" in families    # ISKpn7
        assert "IS1182" in families  # ISKpn6


class TestInferIS26Composite:
    def test_opposite_orientation_pair_is_not_a_pseudocomposite(self):
        feats = pKpQIL_features()
        composites = infer_is26_composites(feats)
        assert composites == []
