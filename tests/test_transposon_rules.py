"""Unit tests for the new transposon inference rules.

Covers ONE_ENDED_RULES (ISEcp1 capture), ROLLING_CIRCLE_RULES (ISCR),
SIGNATURE_RULES (Tn1546 vanA signature, Tn1331 multi-AMR), ISApl1-mcr-1
composite, and the res-site annotator.
"""

from matryoshka.detect import MGEFeature
from matryoshka.transposon import (
    FLANKED_RULES,
    ONE_ENDED_RULES,
    ROLLING_CIRCLE_RULES,
    SIGNATURE_RULES,
    TN3_FAMILY_MEMBERS,
    annotate_res_sites,
    infer_flanked,
    infer_one_ended,
    infer_signature,
    infer_transposons,
)


def is_elem(family, start, end, name=None, complete=True):
    return MGEFeature(
        element_type="IS", family=family, name=name or f"{family}_elem",
        start=start, end=end, strand="+",
        attributes={"type": "c" if complete else "p"},
    )


def amr(name, start, end):
    return MGEFeature(
        element_type="AMR", family="", name=name,
        start=start, end=end, strand="+",
    )


# ---------------------------------------------------------------------------
# Tn6330 (ISApl1 composite)
# ---------------------------------------------------------------------------

class TestTnMcr1:
    def _rule(self):
        return next(r for r in FLANKED_RULES if r.family == "Tn6330")

    def test_two_isapl1_flanking_mcr1(self):
        feats = [
            is_elem("IS30", 1000, 2000, "ISApl1_a"),
            amr("mcr-1", 2100, 3600),
            is_elem("IS30", 3700, 4700, "ISApl1_b"),
        ]
        out = infer_flanked(self._rule(), feats)
        assert len(out) == 1
        assert out[0].family == "Tn6330"
        assert out[0].start == 1000 and out[0].end == 4700

    def test_only_one_isapl1_no_composite(self):
        feats = [
            is_elem("IS30", 1000, 2000, "ISApl1_a"),
            amr("mcr-1", 2100, 3600),
        ]
        assert infer_flanked(self._rule(), feats) == []

    def test_partial_isapl1_excluded(self):
        feats = [
            is_elem("IS30", 1000, 2000, "ISApl1_a", complete=False),
            amr("mcr-1", 2100, 3600),
            is_elem("IS30", 3700, 4700, "ISApl1_b"),
        ]
        assert infer_flanked(self._rule(), feats) == []


# ---------------------------------------------------------------------------
# Tn4001 (IS256 + AAC(6')-Ie/APH(2'')-Ia) — Staphylococcal aminoglycoside
# ---------------------------------------------------------------------------

class TestTn4001Composite:
    def _rule(self):
        return next(r for r in FLANKED_RULES if r.family == "Tn4001")

    def test_is256_flanking_aac6(self):
        feats = [
            is_elem("IS256", 5, 1319, "IS256_a"),
            amr("aac(6')-Ie/aph(2'')-Ia", 1712, 3148),
            is_elem("IS256", 3235, 4550, "IS256_b"),
        ]
        out = infer_flanked(self._rule(), feats)
        assert len(out) == 1
        assert out[0].family == "Tn4001"
        assert out[0].start == 5 and out[0].end == 4550

    def test_partial_is256_excluded(self):
        feats = [
            is_elem("IS256", 5, 1319, "IS256_a", complete=False),
            amr("aac(6')-Ie/aph(2'')-Ia", 1712, 3148),
            is_elem("IS256", 3235, 4550, "IS256_b"),
        ]
        assert infer_flanked(self._rule(), feats) == []

    def test_only_one_is256_no_composite(self):
        feats = [
            is_elem("IS256", 5, 1319, "IS256_a"),
            amr("aac(6')-Ie/aph(2'')-Ia", 1712, 3148),
        ]
        assert infer_flanked(self._rule(), feats) == []


# ---------------------------------------------------------------------------
# Tn2006 (ISAba1 + blaOXA-23) — Acinetobacter carbapenem resistance
# ---------------------------------------------------------------------------

class TestTn2006:
    def _rule(self):
        return next(r for r in FLANKED_RULES if r.family == "Tn2006")

    def test_isaba1_flanking_oxa23(self):
        feats = [
            is_elem("IS4", 1000, 2300, "ISAba1_a"),
            amr("blaOXA-23", 2500, 3400),
            is_elem("IS4", 3500, 4800, "ISAba1_b"),
        ]
        out = infer_flanked(self._rule(), feats)
        assert len(out) == 1
        assert out[0].family == "Tn2006"


# ---------------------------------------------------------------------------
# Tn125 (ISAba125 + blaNDM) — Acinetobacter NDM
# ---------------------------------------------------------------------------

class TestTn125:
    def _rule(self):
        return next(r for r in FLANKED_RULES if r.family == "Tn125")

    def test_isaba125_flanking_ndm(self):
        feats = [
            is_elem("IS30", 1000, 2000, "ISAba125_a"),
            amr("blaNDM-1", 2500, 3300),
            is_elem("IS30", 3800, 4800, "ISAba125_b"),
        ]
        out = infer_flanked(self._rule(), feats)
        assert len(out) == 1
        assert out[0].family == "Tn125"


# ---------------------------------------------------------------------------
# ISEcp1 one-ended capture
# ---------------------------------------------------------------------------

class TestISEcp1Capture:
    def _rule(self):
        return next(r for r in ONE_ENDED_RULES if r.family == "ISEcp1_capture")

    def test_isecp1_upstream_of_ctx_m(self):
        feats = [
            is_elem("IS1380", 1000, 2600, "ISEcp1"),
            amr("blaCTX-M-15", 2800, 3700),
        ]
        out = infer_one_ended(self._rule(), feats)
        assert len(out) == 1
        assert out[0].family == "ISEcp1_capture"
        assert out[0].start == 1000 and out[0].end == 3700

    def test_isecp1_downstream_also_counts(self):
        feats = [
            amr("blaCTX-M-15", 1000, 1900),
            is_elem("IS1380", 2100, 3700, "ISEcp1"),
        ]
        out = infer_one_ended(self._rule(), feats)
        assert len(out) == 1
        assert out[0].start == 1000 and out[0].end == 3700

    def test_no_match_on_blatem(self):
        # Bug-regression: "LAT" substring was hitting blaTEM
        feats = [
            is_elem("IS1380", 1000, 2600, "ISEcp1"),
            amr("blaTEM", 2800, 3700),
        ]
        assert infer_one_ended(self._rule(), feats) == []

    def test_gap_too_large(self):
        feats = [
            is_elem("IS1380", 1000, 2600, "ISEcp1"),
            amr("blaCTX-M-15", 10_000, 10_800),
        ]
        assert infer_one_ended(self._rule(), feats) == []


# ---------------------------------------------------------------------------
# ISCR rolling-circle capture
# ---------------------------------------------------------------------------

class TestISCRCapture:
    def _rule(self):
        return next(r for r in ROLLING_CIRCLE_RULES if r.family == "ISCR_capture")

    def test_iscr1_with_sul2(self):
        feats = [
            is_elem("ISCR1", 1000, 2700, "ISCR1"),
            amr("sul2", 2800, 3600),
        ]
        out = infer_one_ended(self._rule(), feats)
        assert len(out) == 1
        assert out[0].family == "ISCR_capture"

    def test_iscr2_prefix_matches(self):
        feats = [
            is_elem("ISCR2", 1000, 2700, "ISCR2"),
            amr("sul2", 2800, 3600),
        ]
        out = infer_one_ended(self._rule(), feats)
        assert len(out) == 1


# ---------------------------------------------------------------------------
# Tn1546 / Tn1331 signature rules
# ---------------------------------------------------------------------------

class TestTn1546:
    def _rule(self):
        return next(r for r in SIGNATURE_RULES if r.family == "Tn1546")

    def test_vana_present_emits_tn1546(self):
        feats = [amr("vanA", 5000, 6100)]
        out = infer_signature(self._rule(), feats)
        assert len(out) == 1
        assert out[0].family == "Tn1546"
        assert out[0].attributes.get("signature_genes") == "vanA"

    def test_no_vana_no_tn1546(self):
        feats = [amr("blaTEM", 5000, 5800)]
        assert infer_signature(self._rule(), feats) == []


class TestTn1331Signature:
    def _rule(self):
        return next(r for r in SIGNATURE_RULES if r.family == "Tn1331")

    def test_all_three_genes_cluster(self):
        feats = [
            amr("aac(6')-Ib", 1000, 1500),
            amr("blaOXA-9", 2000, 2800),
            amr("aadA1", 3000, 3800),
        ]
        out = infer_signature(self._rule(), feats)
        assert len(out) == 1
        assert out[0].start == 1000 and out[0].end == 3800

    def test_missing_one_gene_no_match(self):
        feats = [
            amr("aac(6')-Ib", 1000, 1500),
            amr("blaOXA-9", 2000, 2800),
            # no aadA1
        ]
        assert infer_signature(self._rule(), feats) == []

    def test_genes_too_far_apart(self):
        feats = [
            amr("aac(6')-Ib", 1000, 1500),
            amr("blaOXA-9", 2000, 2800),
            amr("aadA1", 50_000, 50_800),   # way beyond max_span
        ]
        assert infer_signature(self._rule(), feats) == []


# ---------------------------------------------------------------------------
# Res-site annotator
# ---------------------------------------------------------------------------

class TestResSite:
    def test_tn4401_gets_res_site(self):
        tn = MGEFeature(
            element_type="transposon", family="Tn4401", name="Tn4401b",
            start=5000, end=15000, strand="+",
        )
        res = annotate_res_sites([tn])
        assert len(res) == 1
        assert res[0].element_type == "res_site"
        assert res[0].family == "Tn4401"
        assert 5000 < res[0].start < res[0].end < 15000

    def test_is26_island_no_res_site(self):
        tn = MGEFeature(
            element_type="transposon", family="IS26_island", name="IS26_island_1",
            start=5000, end=15000, strand=".",
        )
        assert annotate_res_sites([tn]) == []

    def test_tn1999_no_res_site(self):
        # Tn1999 is a composite, not a Tn3-family unit transposon
        tn = MGEFeature(
            element_type="transposon", family="Tn1999", name="Tn1999",
            start=5000, end=15000, strand="-",
        )
        assert "Tn1999" not in TN3_FAMILY_MEMBERS
        assert annotate_res_sites([tn]) == []

    def test_reverse_strand_res_at_end(self):
        tn = MGEFeature(
            element_type="transposon", family="Tn4401", name="Tn4401",
            start=5000, end=15000, strand="-",
        )
        res = annotate_res_sites([tn])[0]
        # For a - strand element the res site should be near the element end
        assert res.end > 14_000


# ---------------------------------------------------------------------------
# Top-level infer_transposons dedup
# ---------------------------------------------------------------------------

class TestInferTransposonsIntegration:
    def test_dedup_across_rules(self):
        # Same Tn4401 structure — should emit exactly one Tn4401 feature
        feats = [
            is_elem("IS21", 1000, 2500, "ISKpn7"),
            amr("blaKPC", 2600, 3500),
            is_elem("IS1182", 3600, 5000, "ISKpn6"),
        ]
        out = infer_transposons(feats)
        tn4401 = [f for f in out if f.family == "Tn4401"]
        assert len(tn4401) == 1

    def test_multiple_rules_fire(self):
        # Tn4401 + vanA signature + IS26 island all present in one feature set
        feats = [
            is_elem("IS21", 1000, 2500, "ISKpn7"),
            amr("blaKPC", 2600, 3500),
            is_elem("IS1182", 3600, 5000, "ISKpn6"),
            amr("vanA", 10_000, 11_000),
            is_elem("IS6", 20_000, 21_000, "IS26_a"),
            amr("blaTEM", 22_000, 23_000),
            is_elem("IS6", 24_000, 25_000, "IS26_b"),
        ]
        families = {f.family for f in infer_transposons(feats)}
        assert "Tn4401" in families
        assert "Tn1546" in families
        assert "IS26_island" in families
