"""Tests for reference_scan.py — BLAST-based MGE detection."""

from pathlib import Path

import pytest

from matryoshka.reference_scan import (
    REFERENCES_DIR,
    blast_available,
    scan,
    scan_all,
)

DATA = Path(__file__).parent.parent / "data"
PKPQIL = DATA / "reference_plasmids" / "pKpQIL.fasta"
POXA48A = DATA / "reference_plasmids" / "pOXA-48a.fasta"
TN4401_REF = REFERENCES_DIR / "tn4401.fasta"

SKIP_BLAST = pytest.mark.skipif(
    not blast_available(), reason="blastn not on PATH",
)
SKIP_REFS = pytest.mark.skipif(
    not TN4401_REF.exists(),
    reason="reference sequences not downloaded — run scripts/fetch_mge_references.py",
)
SKIP_QUERY = pytest.mark.skipif(
    not (PKPQIL.exists() and POXA48A.exists()),
    reason="query plasmid FASTAs not present (data/ is gitignored)",
)


@SKIP_BLAST
@SKIP_REFS
@SKIP_QUERY
class TestScan:
    def test_tn4401_detected_on_pkpqil(self):
        hits = scan(PKPQIL, TN4401_REF, min_identity=95.0, min_length=4_000)
        assert len(hits) >= 1
        h = hits[0]
        assert h.family == "Tn4401"
        assert h.attributes["source"] == "reference_scan"
        assert float(h.attributes["blast_identity"]) > 98.0

    def test_no_tn4401_on_poxa48a(self):
        hits = scan(POXA48A, TN4401_REF, min_identity=95.0, min_length=4_000)
        # pOXA-48a carries Tn1999, not Tn4401
        assert hits == []


@SKIP_BLAST
@SKIP_REFS
class TestScanAll:
    def test_pek499_finds_isecp1(self):
        pek499 = DATA / "reference_plasmids" / "pEK499.fasta"
        if not pek499.exists():
            pytest.skip("pEK499 fasta not present")
        hits = scan_all(pek499)
        isecp1 = [h for h in hits if h.name == "ISEcp1"]
        # pEK499 is known to carry multiple ISEcp1 copies (ISEScan misses them)
        assert len(isecp1) >= 1

    def test_no_spurious_acinetobacter_on_ecoli(self):
        pek499 = DATA / "reference_plasmids" / "pEK499.fasta"
        if not pek499.exists():
            pytest.skip("pEK499 fasta not present")
        hits = scan_all(pek499)
        acineto = [h for h in hits if h.family.startswith("Tn60")]
        # With tightened thresholds, shared Tn3 backbones should not
        # appear as Acinetobacter island calls in E. coli plasmids
        assert acineto == []
