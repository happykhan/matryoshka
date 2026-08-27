"""
boundaries.py — Confirm element boundaries using TSD and IR detection.

For each predicted MGEFeature, look for:
  1. Target site duplications (TSDs) flanking the element
  2. Inverted repeats (IRs) at element termini

Coordinates are 1-based inclusive (GFF3 convention). Converted to 0-based
internally before slicing the sequence string.

TSD detection scans a small offset window (±OFFSET_WINDOW bp) on both
sides of the predicted boundary to tolerate off-by-a-few errors from
upstream callers like ISEScan. The *best* match (exact, longest) is
returned — ties are broken by the match closest to the predicted end.
"""

from __future__ import annotations

from dataclasses import dataclass

from Bio.Seq import Seq

from .detect import TSD_LENGTHS, MGEFeature

# Default maximum drift for an upstream caller with independent terminal-repeat
# evidence.  Compound elements can supply a small, explicit family/caller
# window when the terminal IS calls are known to exclude a few outer bases.
OFFSET_WINDOW = 3

@dataclass(frozen=True)
class TSDMatch:
    """One exact repeat pair immediately outside candidate element ends."""

    sequence: str
    left_start: int
    left_end: int
    right_start: int
    right_end: int
    start_offset: int
    end_offset: int


def find_tsd_match(
    seq: str,
    feature: MGEFeature,
    offset_window: int = 0,
) -> TSDMatch | None:
    """Return the nearest exact, boundary-adjacent TSD candidate.

    Candidate repeats must sit flush outside a candidate pair of element
    boundaries. ``offset_window`` is only appropriate when independent end
    evidence (for example terminal IRs) justifies correcting an upstream
    caller's coordinates.
    """
    expected_len = feature.tsd_length or TSD_LENGTHS.get(feature.family)
    if expected_len is None or expected_len < 2:
        return None

    start0 = feature.start - 1
    end0 = feature.end
    candidates: list[TSDMatch] = []
    offsets = range(-offset_window, offset_window + 1)
    for start_offset in offsets:
        candidate_start0 = start0 + start_offset
        left_start0 = candidate_start0 - expected_len
        if left_start0 < 0 or candidate_start0 > len(seq):
            continue
        left = seq[left_start0:candidate_start0]
        if len(left) != expected_len or "N" in left.upper():
            continue
        for end_offset in offsets:
            right_start0 = end0 + end_offset
            if right_start0 < 0 or right_start0 + expected_len > len(seq):
                continue
            right = seq[right_start0:right_start0 + expected_len]
            if "N" in right.upper() or left.upper() != right.upper():
                continue
            candidates.append(TSDMatch(
                sequence=left.upper(),
                left_start=left_start0 + 1,
                left_end=candidate_start0,
                right_start=right_start0 + 1,
                right_end=right_start0 + expected_len,
                start_offset=start_offset,
                end_offset=end_offset,
            ))

    if not candidates:
        return None
    return min(
        candidates,
        key=lambda match: (
            abs(match.start_offset) + abs(match.end_offset),
            abs(match.start_offset - match.end_offset),
            match.start_offset,
            match.end_offset,
        ),
    )


def find_tsd(
    seq: str,
    feature: MGEFeature,
    flank: int = 0,
    right_offset: int = 0,
) -> str | None:
    """
    Look for a target site duplication flanking the element.

    Compatibility wrapper returning only the repeat sequence. ``flank`` is
    retained for callers of the former API but no longer broadens the search:
    non-adjacent repeats are not valid target-site duplications. Use
    ``right_offset`` to request a symmetric boundary window when independent
    end evidence supports coordinate refinement.
    """
    del flank
    match = find_tsd_match(seq, feature, offset_window=right_offset)
    return match.sequence if match else None


def find_ir(
    seq: str,
    feature: MGEFeature,
    ir_len: int = 20,
    mismatch: int = 2,
) -> tuple[str, str] | None:
    """
    Look for inverted repeats at element termini.
    Returns (IRL, IRR) sequences if found, None otherwise.
    """
    start0 = feature.start - 1
    end0 = feature.end

    if end0 - start0 < 2 * ir_len:
        # Element too short for meaningful IR at both ends
        return None

    irl = seq[start0: start0 + ir_len]
    irr_region = seq[end0 - ir_len: end0]
    if len(irl) < ir_len or len(irr_region) < ir_len:
        return None
    irr_rc = str(Seq(irr_region).reverse_complement())

    mismatches = sum(a != b for a, b in zip(irl, irr_rc, strict=True))
    if mismatches <= mismatch:
        return irl, irr_region
    return None


def _apply_ml_refinement(seq: str, features: list[MGEFeature]) -> None:
    try:
        from .boundary_refine import refine_boundaries

        for feature in features:
            if (
                feature.start > 0
                and feature.end > 0
                and feature.attributes.get("source") == "reference_scan"
                and feature.end - feature.start > 100
            ):
                feature.start, feature.end = refine_boundaries(
                    seq, feature.start, feature.end, search_window=80, step=3
                )
    except Exception:
        # The experimental model is optional and disabled by default. A missing
        # or incompatible model must not disable deterministic IR/TSD checks.
        return


def _annotate_ir(seq: str, feature: MGEFeature) -> None:
    ir_length = int(feature.attributes.get("ir_length", 20) or 20)
    ir_mismatch = max(2, round(ir_length * 0.05))
    ir = find_ir(seq, feature, ir_len=ir_length, mismatch=ir_mismatch)
    if ir:
        feature.ir_left, feature.ir_right = ir
        feature.attributes["ir_evidence"] = "sequence_matched"


def _offset_window(feature: MGEFeature) -> int:
    if feature.attributes.get("source") == "reference_scan":
        return 0
    if not (feature.ir_left and feature.ir_right):
        return 0
    return int(
        feature.attributes.get("boundary_offset_window", OFFSET_WINDOW)
        or OFFSET_WINDOW
    )


def _record_tsd_match(
    feature: MGEFeature,
    match: TSDMatch,
    expected_tsd: int,
    offset_window: int,
) -> None:
    predicted_start, predicted_end = feature.start, feature.end
    feature.attributes.update({
        "tsd_candidate_seq": match.sequence,
        "tsd_left_start": match.left_start,
        "tsd_left_end": match.left_end,
        "tsd_right_start": match.right_start,
        "tsd_right_end": match.right_end,
        "tsd_start_offset": match.start_offset,
        "tsd_end_offset": match.end_offset,
    })
    if expected_tsd < 4:
        feature.attributes["tsd_evidence"] = "short_repeat_candidate"
        feature.attributes["tsd_evidence_strength"] = "weak"
        return

    feature.tsd_seq = match.sequence
    feature.attributes["tsd_evidence"] = "sequence_matched"
    feature.attributes["tsd_evidence_strength"] = (
        "strong" if expected_tsd >= 8 else "supporting"
    )
    if offset_window and (match.start_offset or match.end_offset):
        feature.attributes.update({
            "predicted_start": predicted_start,
            "predicted_end": predicted_end,
            "boundary_refinement": "sequence_matched_TSD",
        })
        feature.start = predicted_start + match.start_offset
        feature.end = predicted_end + match.end_offset


def _record_missing_tsd(seq: str, feature: MGEFeature, expected_tsd: int) -> None:
    has_flanks = (
        feature.start - 1 >= expected_tsd
        and feature.end + expected_tsd <= len(seq)
    )
    feature.attributes["tsd_evidence"] = (
        "searched_not_found" if has_flanks else "untestable_missing_flank"
    )


def confirm_boundaries(
    seq: str,
    features: list[MGEFeature],
    ml_refine: bool = False,
) -> list[MGEFeature]:
    """Annotate TSD and IR evidence on each feature in-place.

    ML refinement is disabled by default because its evaluation worsened
    boundary accuracy. Deterministic sequence evidence remains authoritative.
    """
    if ml_refine:
        _apply_ml_refinement(seq, features)

    for feature in features:
        if feature.start <= 0 or feature.end <= 0:
            continue
        expected_tsd = feature.tsd_length or TSD_LENGTHS.get(feature.family)
        if expected_tsd:
            feature.tsd_length = expected_tsd
        _annotate_ir(seq, feature)
        offset_window = _offset_window(feature)
        match = find_tsd_match(seq, feature, offset_window=offset_window)
        if match and expected_tsd is not None:
            _record_tsd_match(feature, match, expected_tsd, offset_window)
        elif expected_tsd:
            _record_missing_tsd(seq, feature, expected_tsd)
    return features
