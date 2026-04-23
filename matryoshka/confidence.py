"""
confidence.py — Assign a confidence score to each annotated MGEFeature.

Score range: 0.0 – 1.0, higher = stronger evidence.

Scoring rubric:

    1.00  BLAST reference hit with identity ≥98% and subject coverage ≥95%
    0.90  BLAST reference hit with identity ≥95% and coverage ≥80%
    0.85  Flanked-cargo rule with BOTH TSD and IR confirmed
    0.80  Flanked-cargo rule with TSD confirmed
    0.75  Flanked-cargo rule with IR only, no TSD
    0.70  Flanked-cargo rule, no boundary evidence
    0.65  One-ended rule (ISEcp1 / ISCR capture) with some boundary evidence
    0.55  One-ended rule, no boundary evidence
    0.50  Signature rule (gene-pattern only, no IS anchor)
    0.40  IS26 translocatable unit (single-IS, speculative)
    0.30  res site (positional, not sequence-confirmed)
    0.60  IS element with IR confirmed (ISEScan baseline)
    0.45  IS element without IR
    0.70  AMR / cassette / attC / integrase (detection-tool call, trusted)
    0.70  integron (IntegronFinder)

The score is stored in ``feature.attributes['confidence']`` as a float
rounded to 2 dp, and a human-readable label in
``feature.attributes['confidence_label']`` for GFF3 / GenBank output.
"""

from __future__ import annotations

from .detect import MGEFeature


def _label(score: float) -> str:
    if score >= 0.85:
        return "high"
    if score >= 0.65:
        return "medium"
    if score >= 0.45:
        return "low"
    return "speculative"


def _blast_score(f: MGEFeature) -> float:
    ident_s = f.attributes.get("blast_identity")
    cov_s = (
        f.attributes.get("blast_subject_coverage")
        or f.attributes.get("blast_coverage")
    )
    try:
        ident = float(ident_s) if ident_s is not None else 0.0
        cov = float(cov_s) if cov_s is not None else 0.0
    except (TypeError, ValueError):
        return 0.8
    if ident >= 98.0 and cov >= 95.0:
        return 1.00
    if ident >= 95.0 and cov >= 80.0:
        return 0.90
    if ident >= 90.0 and cov >= 60.0:
        return 0.80
    return 0.70


def _rule_score(f: MGEFeature) -> float:
    family = f.family
    has_tsd = bool(f.tsd_seq)
    has_ir = bool(f.ir_left or f.ir_right)

    # Single-IS translocatable unit is the least-confident call
    if family == "IS26_TU":
        return 0.40

    # Signature-only rules (vanA Tn1546, Tn1331 three-gene)
    if family in ("Tn1546", "Tn1331"):
        return 0.50

    # One-ended / rolling-circle capture
    if family in ("ISEcp1_capture", "ISCR_capture", "IS91_capture"):
        return 0.65 if (has_tsd or has_ir) else 0.55

    # Genomic islands emitted by compound detection
    if f.element_type == "genomic_island":
        return 0.75

    # Flanked-cargo composites (Tn4401, Tn1999, Tn_mcr1, IS26_island)
    if f.element_type == "transposon":
        if has_tsd and has_ir:
            return 0.85
        if has_tsd:
            return 0.80
        if has_ir:
            return 0.75
        return 0.70

    return 0.65


def _detection_score(f: MGEFeature) -> float:
    if f.element_type == "IS":
        return 0.60 if (f.ir_left or f.ir_right) else 0.45
    if f.element_type == "res_site":
        return 0.30
    if f.element_type in ("AMR", "cassette", "attC", "integrase",
                          "integron", "replicon"):
        return 0.70
    return 0.50


def assign_confidence(features: list[MGEFeature]) -> list[MGEFeature]:
    """Annotate every feature (recursively including children) with a
    confidence score and label. In-place mutation; returns the input
    for chaining.
    """
    seen: set[int] = set()

    def _visit(f: MGEFeature) -> None:
        fid = id(f)
        if fid in seen:
            return
        seen.add(fid)

        if f.attributes.get("source") == "reference_scan":
            score = _blast_score(f)
        elif f.attributes.get("source") in ("rule", "inference"):
            score = _rule_score(f)
        elif f.element_type == "transposon":
            # Transposons without an explicit source came from an
            # inference rule — use the rule scorer
            score = _rule_score(f)
        else:
            score = _detection_score(f)

        # Respect explicit low-confidence flag on the feature
        if f.attributes.get("confidence") == "low":
            score = min(score, 0.40)

        f.attributes["confidence"] = round(score, 2)
        f.attributes["confidence_label"] = _label(score)
        for c in f.children:
            _visit(c)

    for f in features:
        _visit(f)
    return features
