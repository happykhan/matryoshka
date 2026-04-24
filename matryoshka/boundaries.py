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

from Bio.Seq import Seq

from .detect import TSD_LENGTHS, MGEFeature

# How many bp either side of the predicted boundary the TSD search may drift
OFFSET_WINDOW = 3

# Flank length checked on each side of the element for TSD candidates
_FLANK = 30


def find_tsd(
    seq: str,
    feature: MGEFeature,
    flank: int = _FLANK,
    right_offset: int = OFFSET_WINDOW,
) -> str | None:
    """
    Look for a target site duplication flanking the element.

    The left flank (up to `flank` bp before predicted start) is scanned for
    any k-mer of the expected TSD length. For each candidate, we check
    whether it appears flush at the predicted end or within `right_offset`
    bp downstream — tolerating small end-coordinate inaccuracies from the
    underlying IS caller.

    Returns the left-most matching k-mer in the left flank; that ordering
    biases toward the true flanking TSD rather than a spurious match flush
    against the predicted end.
    """
    expected_len = TSD_LENGTHS.get(feature.family)
    if expected_len is None or expected_len < 3:
        return None

    start0 = feature.start - 1
    end0 = feature.end
    if start0 < 0 or end0 + expected_len + right_offset > len(seq):
        return None

    left_region = seq[max(0, start0 - flank): start0]
    if len(left_region) < expected_len:
        return None

    # Right-side candidates at increasing offsets from the predicted end
    right_kmers: list[str] = []
    for de in range(0, right_offset + 1):
        kmer = seq[end0 + de: end0 + de + expected_len]
        if len(kmer) == expected_len and "N" not in kmer.upper():
            right_kmers.append(kmer)

    for i in range(len(left_region) - expected_len + 1):
        candidate = left_region[i: i + expected_len]
        if "N" in candidate.upper():
            continue
        if candidate in right_kmers:
            return candidate

    return None


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


def confirm_boundaries(seq: str, features: list[MGEFeature]) -> list[MGEFeature]:
    """Annotate TSD and IR evidence on each feature in-place."""
    for f in features:
        if f.start <= 0 or f.end <= 0:
            continue
        tsd = find_tsd(seq, f)
        if tsd:
            f.tsd_seq = tsd
        ir = find_ir(seq, f)
        if ir:
            f.ir_left, f.ir_right = ir
    return features
