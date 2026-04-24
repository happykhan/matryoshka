"""Tests for boundaries.py TSD and IR detection."""

from matryoshka.boundaries import confirm_boundaries, find_ir, find_tsd
from matryoshka.detect import MGEFeature


def make_feature(family, start, end, strand="+") -> MGEFeature:
    return MGEFeature(
        element_type="IS",
        family=family,
        name=family + "_test",
        start=start,
        end=end,
        strand=strand,
    )


class TestFindTSD:
    def test_exact_tsd_found(self):
        # IS6 family expects 8bp TSD
        # Build sequence: [left_flank][TSD][element][TSD][right_flank]
        tsd = "ATCGATCG"
        element = "N" * 820
        seq = "XXXXXXXXXXXX" + tsd + element + tsd + "XXXXXXXXXXXX"
        # element starts at offset 12 + 8 = 20 (1-based: 21)
        start = 12 + 8 + 1   # 1-based start of element
        end = start + len(element) - 1
        f = make_feature("IS6", start, end)
        assert find_tsd(seq, f) == tsd

    def test_no_tsd_returns_none(self):
        seq = "AAAAAAAAAAAAAAAAAAAAA" + "N" * 820 + "TTTTTTTTTTTTTTTTTTTTT"
        f = make_feature("IS6", 22, 841)
        assert find_tsd(seq, f) is None

    def test_rolling_circle_skipped(self):
        # IS91 has no TSD — should return None regardless
        seq = "ATCGATCGATCG" + "N" * 820 + "ATCGATCGATCG"
        f = make_feature("IS91", 13, 832)
        assert find_tsd(seq, f) is None

    def test_unknown_family_skipped(self):
        seq = "ATCGATCG" + "N" * 820 + "ATCGATCG"
        f = make_feature("unknown_family", 9, 828)
        assert find_tsd(seq, f) is None


class TestFindIR:
    def test_perfect_ir_found(self):
        # Build element with perfect inverted repeats at termini
        irl = "GGCACTGTTGCAAA"
        irr = str(__import__("Bio.Seq", fromlist=["Seq"]).Seq(irl).reverse_complement())
        element = irl + "N" * 792 + irr
        padding = "X" * 50
        seq = padding + element + padding
        start = len(padding) + 1  # 1-based
        end = start + len(element) - 1
        f = make_feature("IS6", start, end)
        result = find_ir(seq, f, ir_len=14, mismatch=0)
        assert result is not None
        assert result[0] == irl

    def test_ir_with_mismatches(self):
        irl = "GGCACTGTTGCAAA"
        # 2 mismatches at positions 0 and 1
        irr_mutated = "TTCACTGTTGCAAA"
        irr_rc = str(__import__("Bio.Seq", fromlist=["Seq"]).Seq(irr_mutated).reverse_complement())
        element = irl + "N" * 792 + irr_rc
        padding = "X" * 50
        seq = padding + element + padding
        start = len(padding) + 1
        end = start + len(element) - 1
        f = make_feature("IS6", start, end)
        result = find_ir(seq, f, ir_len=14, mismatch=2)
        assert result is not None

    def test_no_ir_returns_none(self):
        # Random sequence — no inverted repeats
        element = "A" * 820
        seq = "X" * 50 + element + "X" * 50
        f = make_feature("IS6", 51, 870)
        result = find_ir(seq, f)
        assert result is None


class TestConfirmBoundaries:
    def test_annotates_in_place(self):
        tsd = "ATCGATCG"
        element = "N" * 820
        seq = "X" * 20 + tsd + element + tsd + "X" * 20
        start = 29   # 1-based: 20 padding + 8 TSD + 1
        end = start + len(element) - 1
        f = make_feature("IS6", start, end)
        result = confirm_boundaries(seq, [f])
        assert result[0].tsd_seq == tsd

    def test_ir_annotated(self):
        from Bio.Seq import Seq
        irl = "GGCACTGTTGCAAA"
        irr = str(Seq(irl).reverse_complement())
        element = irl + "N" * 792 + irr
        seq = "X" * 50 + element + "X" * 50
        start = 51
        end = start + len(element) - 1
        f = make_feature("IS6", start, end)
        confirm_boundaries(seq, [f])
        # ir_left is the first ir_len(=20) bases from element start; IRL occupies first 14
        assert f.ir_left is not None
        assert f.ir_left[:len(irl)] == irl
