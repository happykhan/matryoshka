"""Tests for boundaries.py TSD and IR detection."""

import random

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

    def test_feature_specific_tsd_length_is_used_for_named_family(self):
        tsd = "AACGT"
        element = "N" * 200
        seq = "X" * 20 + tsd + element + tsd + "X" * 20
        f = make_feature("Tn1", 26, 225)
        f.tsd_length = 5
        assert find_tsd(seq, f) == tsd

    def test_non_adjacent_flank_match_is_not_a_tsd(self):
        element = "N" * 200
        seq = "GAGCC" + "A" * 6 + element + "GAGCC" + "A" * 20
        f = make_feature("Tn1", 12, 211)
        f.tsd_length = 5
        assert find_tsd(seq, f) is None

    def test_random_five_base_false_match_rate_is_below_one_percent(self):
        rng = random.Random(20260827)
        f = make_feature("Tn1", 41, 140)
        f.tsd_length = 5
        positives = 0
        trials = 5_000
        for _ in range(trials):
            seq = "".join(rng.choice("ACGT") for _ in range(180))
            positives += find_tsd(seq, f) is not None
        assert positives / trials < 0.01


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
        assert result[0].attributes["tsd_evidence"] == "sequence_matched"
        assert result[0].attributes["tsd_evidence_strength"] == "strong"
        assert result[0].attributes["tsd_left_start"] == 21
        assert result[0].attributes["tsd_right_start"] == 849

    def test_short_repeat_is_weak_candidate_not_confirmed_tsd(self):
        tsd = "AT"
        element = "N" * 200
        seq = "A" * 20 + tsd + element + tsd + "A" * 20
        f = make_feature("IS30", 23, 222)
        confirm_boundaries(seq, [f])
        assert f.tsd_seq is None
        assert f.attributes["tsd_candidate_seq"] == tsd
        assert f.attributes["tsd_evidence"] == "short_repeat_candidate"
        assert f.attributes["tsd_evidence_strength"] == "weak"

    def test_records_when_expected_tsd_was_searched_but_not_found(self):
        seq = "A" * 30 + "N" * 820 + "T" * 30
        f = make_feature("IS6", 31, 850)
        confirm_boundaries(seq, [f])
        assert f.tsd_seq is None
        assert f.attributes["tsd_evidence"] == "searched_not_found"

    def test_records_when_flanks_are_unavailable(self):
        f = make_feature("IS6", 1, 820)
        confirm_boundaries("N" * 820, [f])
        assert f.attributes["tsd_evidence"] == "untestable_missing_flank"

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
