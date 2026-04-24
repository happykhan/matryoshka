"""Tests for Tn1999 inference and TSD detection on pOXA-48a."""

from pathlib import Path

import pytest
from Bio import SeqIO

from matryoshka.boundaries import confirm_boundaries
from matryoshka.detect import parse_amrfinder, parse_isescan
from matryoshka.hierarchy import build_hierarchy
from matryoshka.transposon import infer_tn1999

DATA = Path(__file__).parent.parent / "data"
ISESCAN_TSV  = DATA / "isescan_test/reference_plasmids/pOXA-48a.fasta.tsv"
AMRFINDER_TSV = DATA / "amrfinder_test/pOXA-48a.tsv"
FASTA = DATA / "reference_plasmids/pOXA-48a.fasta"

SKIP = pytest.mark.skipif(
    not (ISESCAN_TSV.exists() and AMRFINDER_TSV.exists() and FASTA.exists()),
    reason="pOXA-48a test data not present",
)


@SKIP
class TestInferTn1999:
    def _feats(self):
        return parse_isescan(ISESCAN_TSV) + parse_amrfinder(AMRFINDER_TSV)

    def test_one_tn1999_found(self):
        tn = infer_tn1999(self._feats())
        assert len(tn) == 1
        assert tn[0].family == "Tn1999"

    def test_tn1999_span(self):
        tn = infer_tn1999(self._feats())[0]
        # IS4_157 (left) start=2704, IS4_157 (right) end=7596
        assert tn.start == 2704
        assert tn.end == 7596

    def test_tn1999_contains_blaoxa(self):
        feats = self._feats()
        tn = infer_tn1999(feats)
        roots = build_hierarchy(feats + tn)
        tn_node = next(r for r in roots if r.family == "Tn1999")
        child_names = [c.name for c in tn_node.children]
        assert "blaOXA-48" in child_names

    def test_tn1999_contains_is4_elements(self):
        feats = self._feats()
        tn = infer_tn1999(feats)
        roots = build_hierarchy(feats + tn)
        tn_node = next(r for r in roots if r.family == "Tn1999")
        is4_children = [c for c in tn_node.children if c.family == "IS4"]
        assert len(is4_children) == 2


@SKIP
class TestTn1999TSD:
    def test_5bp_tsd_found(self):
        feats = parse_isescan(ISESCAN_TSV) + parse_amrfinder(AMRFINDER_TSV)
        tn = infer_tn1999(feats)
        seq = str(SeqIO.read(FASTA, "fasta").seq)
        confirm_boundaries(seq, tn)
        assert tn[0].tsd_seq is not None
        assert len(tn[0].tsd_seq) == 5

    def test_tsd_sequence(self):
        feats = parse_isescan(ISESCAN_TSV) + parse_amrfinder(AMRFINDER_TSV)
        tn = infer_tn1999(feats)
        seq = str(SeqIO.read(FASTA, "fasta").seq)
        confirm_boundaries(seq, tn)
        assert tn[0].tsd_seq == "TGCTG"
