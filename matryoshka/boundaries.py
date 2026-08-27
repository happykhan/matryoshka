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


def confirm_boundaries(
    seq: str,
    features: list[MGEFeature],
    ml_refine: bool = False,
) -> list[MGEFeature]:
    """Annotate TSD and IR evidence on each feature in-place.

    ml_refine is disabled by default: evaluation on TnCentral showed the
    model degraded boundary accuracy (MAE 24.6→39.4bp, 7/9 features worse).
    Enable only if retraining on a larger, more diverse dataset.
    """
    if ml_refine:
        try:
            from .boundary_refine import refine_boundaries
            for f in features:
                if (
                    f.start > 0 and f.end > 0
                    and f.attributes.get("source") == "reference_scan"
                    and f.end - f.start > 100
                ):
                    new_start, new_end = refine_boundaries(
                        seq, f.start, f.end, search_window=80, step=3
                    )
                    if new_start != f.start or new_end != f.end:
                        f.start = new_start
                        f.end = new_end
        except Exception:
            pass  # model absent or error — continue without refinement

    for f in features:
        if f.start <= 0 or f.end <= 0:
            continue
        expected_tsd = f.tsd_length or TSD_LENGTHS.get(f.family)
        if expected_tsd:
            f.tsd_length = expected_tsd
        ir_length = int(f.attributes.get("ir_length", 20) or 20)
        ir_mismatch = max(2, round(ir_length * 0.05))
        ir = find_ir(seq, f, ir_len=ir_length, mismatch=ir_mismatch)
        if ir:
            f.ir_left, f.ir_right = ir
            f.attributes["ir_evidence"] = "sequence_matched"

        # Reference alignments already provide precise ends. Coordinate drift
        # is considered only when terminal-repeat evidence independently
        # supports an upstream caller's boundaries.
        has_independent_ends = bool(f.ir_left and f.ir_right)
        offset_window = 0
        if (
            f.attributes.get("source") != "reference_scan"
            and has_independent_ends
        ):
            offset_window = int(
                f.attributes.get("boundary_offset_window", OFFSET_WINDOW)
                or OFFSET_WINDOW
            )
        match = find_tsd_match(seq, f, offset_window=offset_window)
        if match and expected_tsd is not None:
            predicted_start = f.start
            predicted_end = f.end
            f.attributes.update({
                "tsd_candidate_seq": match.sequence,
                "tsd_left_start": match.left_start,
                "tsd_left_end": match.left_end,
                "tsd_right_start": match.right_start,
                "tsd_right_end": match.right_end,
                "tsd_start_offset": match.start_offset,
                "tsd_end_offset": match.end_offset,
            })
            if expected_tsd >= 4:
                f.tsd_seq = match.sequence
                f.attributes["tsd_evidence"] = "sequence_matched"
                f.attributes["tsd_evidence_strength"] = (
                    "strong" if expected_tsd >= 8 else "supporting"
                )
                # Only evidence-gated searches may alter coordinates.  The
                # exact repeat remains immediately outside the corrected
                # element, so downstream hierarchy and drawing use the
                # resolved biological boundary rather than the raw caller end.
                if offset_window and (
                    match.start_offset != 0 or match.end_offset != 0
                ):
                    f.attributes.update({
                        "predicted_start": predicted_start,
                        "predicted_end": predicted_end,
                        "boundary_refinement": "sequence_matched_TSD",
                    })
                    f.start = predicted_start + match.start_offset
                    f.end = predicted_end + match.end_offset
            else:
                # Equality of a 2-3 bp word is too common to confirm a
                # transposition boundary without additional event evidence.
                f.attributes["tsd_evidence"] = "short_repeat_candidate"
                f.attributes["tsd_evidence_strength"] = "weak"
        elif expected_tsd:
            start0 = f.start - 1
            has_flanks = (
                start0 >= expected_tsd
                and f.end + expected_tsd <= len(seq)
            )
            f.attributes["tsd_evidence"] = (
                "searched_not_found" if has_flanks else "untestable_missing_flank"
            )
    return features
