"""Tests for detect.py parsers against real tool outputs."""

from pathlib import Path
import pytest
from matryoshka.detect import parse_isescan, parse_amrfinder, parse_integron_finder

DATA = Path(__file__).parent.parent / "data"
ISESCAN_TSV = DATA / "isescan_test/reference_plasmids/pEK499.fasta.tsv"
AMRFINDER_TSV = DATA / "amrfinder_test/pEK499.tsv"
INTEGRONS_FILE = DATA / "integron_test/Results_Integron_Finder_pEK499/pEK499.integrons"


@pytest.mark.skipif(not ISESCAN_TSV.exists(), reason="ISEScan test data not present")
class TestParseISEScan:
    def test_count(self):
        feats = parse_isescan(ISESCAN_TSV)
        assert len(feats) == 12

    def test_element_type(self):
        feats = parse_isescan(ISESCAN_TSV)
        assert all(f.element_type == "IS" for f in feats)

    def test_is26_family(self):
        feats = parse_isescan(ISESCAN_TSV)
        is6 = [f for f in feats if f.family == "IS6"]
        assert len(is6) == 7

    def test_coordinates_1based(self):
        feats = parse_isescan(ISESCAN_TSV)
        first_is26 = next(f for f in feats if f.name == "IS6_292")
        assert first_is26.start == 36287
        assert first_is26.end == 37106

    def test_tir_parsed(self):
        feats = parse_isescan(ISESCAN_TSV)
        first_is26 = next(f for f in feats if f.name == "IS6_292")
        assert first_is26.ir_left == "GGCACTGTTGCAAA"
        assert first_is26.ir_right == "GGCACTGTTGCAAA"

    def test_partial_no_tir(self):
        feats = parse_isescan(ISESCAN_TSV)
        partials = [f for f in feats if f.name == "new_4"]
        assert all(f.ir_left is None for f in partials)


@pytest.mark.skipif(not AMRFINDER_TSV.exists(), reason="AMRFinder test data not present")
class TestParseAMRFinder:
    def test_count(self):
        feats = parse_amrfinder(AMRFINDER_TSV)
        assert len(feats) == 12

    def test_element_type(self):
        feats = parse_amrfinder(AMRFINDER_TSV)
        assert all(f.element_type == "AMR" for f in feats)

    def test_blatem_present(self):
        feats = parse_amrfinder(AMRFINDER_TSV)
        names = [f.name for f in feats]
        assert "blaTEM" in names

    def test_blactxm15_present(self):
        feats = parse_amrfinder(AMRFINDER_TSV)
        names = [f.name for f in feats]
        assert "blaCTX-M-15" in names

    def test_strand(self):
        feats = parse_amrfinder(AMRFINDER_TSV)
        for f in feats:
            assert f.strand in ("+", "-", ".")


@pytest.mark.skipif(not INTEGRONS_FILE.exists(), reason="IntegronFinder test data not present")
class TestParseIntegronFinder:
    def test_one_integron(self):
        feats = parse_integron_finder(INTEGRONS_FILE)
        assert len(feats) == 1

    def test_integron_span(self):
        feats = parse_integron_finder(INTEGRONS_FILE)
        integron = feats[0]
        assert integron.start == 65567
        assert integron.end == 68683

    def test_children_count(self):
        feats = parse_integron_finder(INTEGRONS_FILE)
        assert len(feats[0].children) == 6

    def test_child_types(self):
        feats = parse_integron_finder(INTEGRONS_FILE)
        types = {c.element_type for c in feats[0].children}
        assert "integrase" in types
        assert "attC" in types
        assert "cassette" in types

    def test_integrase_coordinates(self):
        feats = parse_integron_finder(INTEGRONS_FILE)
        integrase = next(c for c in feats[0].children if c.element_type == "integrase")
        assert integrase.start == 65567
        assert integrase.end == 66580
