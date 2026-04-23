"""
boundaries.py — Confirm element boundaries using TSD and IR detection.

For each predicted MGEFeature, look for:
  1. Target site duplications (TSDs) flanking the element
  2. Inverted repeats (IRs) at element termini

TSD detection: family-specific length, exact match required at both flanks.
IR detection: reverse complement match at IRL vs IRR termini.
"""

from Bio.Seq import Seq
from .detect import MGEFeature, TSD_LENGTHS


def find_tsd(seq: str, feature: MGEFeature, flank: int = 20) -> str | None:
    """
    Look for a target site duplication flanking the element.
    Returns the TSD sequence if found, None otherwise.
    """
    expected_len = TSD_LENGTHS.get(feature.family)
    if expected_len is None or expected_len < 3:
        return None

    left_flank = seq[max(0, feature.start - flank): feature.start]
    right_flank = seq[feature.end: feature.end + flank]

    # Check all positions in left flank for a match in right flank
    for i in range(len(left_flank) - expected_len + 1):
        candidate = left_flank[i: i + expected_len]
        if right_flank.startswith(candidate):
            return candidate

    return None


def find_ir(seq: str, feature: MGEFeature, ir_len: int = 20, mismatch: int = 2) -> tuple[str, str] | None:
    """
    Look for inverted repeats at element termini.
    Returns (IRL, IRR) sequences if found, None otherwise.
    """
    irl = seq[feature.start: feature.start + ir_len]
    irr_region = seq[feature.end - ir_len: feature.end]
    irr_rc = str(Seq(irr_region).reverse_complement())

    mismatches = sum(a != b for a, b in zip(irl, irr_rc))
    if mismatches <= mismatch:
        return irl, irr_region

    return None


def confirm_boundaries(seq: str, features: list[MGEFeature]) -> list[MGEFeature]:
    """Annotate TSD and IR evidence on each feature."""
    for f in features:
        tsd = find_tsd(seq, f)
        if tsd:
            f.tsd_seq = tsd
        ir = find_ir(seq, f)
        if ir:
            f.ir_left, f.ir_right = ir
    return features
