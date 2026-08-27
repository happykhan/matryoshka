"""Tests that the acceptance proof fails closed when evidence is missing."""

from matryoshka.detect import MGEFeature
from matryoshka.proof import build_tn123_proof


def test_proof_reports_no_loci_instead_of_false_pass():
    proof = build_tn123_proof([("empty", 1000, [])], {})
    assert proof["summary"] == {
        "status": "NO_TN123_LOCI",
        "tn123_loci": 0,
        "passed": 0,
        "partial": 0,
        "failed": 0,
    }


def test_proof_fails_when_required_component_is_missing():
    parent = MGEFeature(
        element_type="transposon",
        family="Tn1",
        name="Tn1 reference-match candidate",
        start=101,
        end=5049,
        strand="+",
        attributes={
            "tn123_canonical": "true",
            "best_match": "Tn1",
            "naming_evidence": ["whole_locus_alignment"],
            "component_assembly_status": "partial_or_conflicting",
            "component_order_valid": False,
            "component_requirements": {
                "terminal_IR": True,
                "blaTEM": True,
                "tnpR": True,
                "res": False,
                "tnpA": True,
            },
            "detected_component_count": 5,
            "reference_projected_component_count": 1,
        },
    )
    proof = build_tn123_proof([("missing_res", 6000, [parent])], {})
    locus = proof["records"][0]["loci"][0]
    assert proof["summary"]["status"] == "FAIL"
    assert locus["verdict"] == "FAIL"
    assert locus["checks"]["all_required_components_present"] is False
    assert locus["checks"]["component_grammar_complete"] is False
    assert locus["checks"]["no_reference_projected_components"] is False
