"""Tests for transposon inference against pKpQIL data."""

from pathlib import Path
import pytest
from matryoshka.detect import parse_isescan, parse_amrfinder, MGEFeature
from matryoshka.transposon import infer_tn4401, infer_is26_composites, infer_transposons

DATA = Path(__file__).parent.parent / "data"
ISESCAN_TSV = DATA / "isescan_test/reference_plasmids/pKpQIL.fasta.tsv"
AMRFINDER_TSV = DATA / "amrfinder_test/pKpQIL.tsv"


def pKpQIL_features():
    return (
        parse_isescan(ISESCAN_TSV) +
        parse_amrfinder(AMRFINDER_TSV)
    )


@pytest.mark.skipif(
    not (ISESCAN_TSV.exists() and AMRFINDER_TSV.exists()),
    reason="pKpQIL test data not present",
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


@pytest.mark.skipif(
    not (ISESCAN_TSV.exists() and AMRFINDER_TSV.exists()),
    reason="pKpQIL test data not present",
)
class TestInferIS26Composite:
    def test_one_composite_found(self):
        feats = pKpQIL_features()
        composites = infer_is26_composites(feats)
        assert len(composites) == 1

    def test_island_cargo_recorded(self):
        feats = pKpQIL_features()
        comp = infer_is26_composites(feats)[0]
        cargo = comp.attributes.get("cargo", "")
        assert "blaOXA" in cargo
        assert "blaTEM" in cargo

    def test_composite_span(self):
        feats = pKpQIL_features()
        comp = infer_is26_composites(feats)[0]
        assert comp.start == 25540   # first IS26
        assert comp.end == 32163     # second IS26

    def test_composite_contains_amr_genes(self):
        feats = pKpQIL_features()
        from matryoshka.hierarchy import build_hierarchy
        composites = infer_is26_composites(feats)
        roots = build_hierarchy(feats + composites)
        comp_node = next(r for r in roots if r.family == "IS26_island")
        amr_children = [c for c in comp_node.children if c.element_type == "AMR"]
        assert len(amr_children) == 2
