"""
boundaries.py — Confirm element boundaries using TSD and IR detection.

For each predicted MGEFeature, look for:
  1. Target site duplications (TSDs) flanking the element
  2. Inverted repeats (IRs) at element termini

Coordinates are 1-based inclusive (GFF3 convention). Converted to 0-based
internally before slicing the sequence string.

TSD detection: family-specific length, exact match in right flank.
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

    # Convert 1-based inclusive to 0-based
    start0 = feature.start - 1
    end0 = feature.end          # 1-based inclusive end == 0-based exclusive end

    left_flank = seq[max(0, start0 - flank): start0]
    right_flank = seq[end0: end0 + flank]

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
    start0 = feature.start - 1
    end0 = feature.end

    irl = seq[start0: start0 + ir_len]
    irr_region = seq[end0 - ir_len: end0]
    irr_rc = str(Seq(irr_region).reverse_complement())

    mismatches = sum(a != b for a, b in zip(irl, irr_rc))
    if mismatches <= mismatch:
        return irl, irr_region

    return None


def confirm_boundaries(seq: str, features: list[MGEFeature]) -> list[MGEFeature]:
    """Annotate TSD and IR evidence on each feature in-place."""
    for f in features:
        tsd = find_tsd(seq, f)
        if tsd:
            f.tsd_seq = tsd
        ir = find_ir(seq, f)
        if ir:
            f.ir_left, f.ir_right = ir
    return features
