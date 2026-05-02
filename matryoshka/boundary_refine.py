"""
boundary_refine.py — ML-based transposon boundary refinement.

BLAST-detected transposons often have boundaries off by 20-100bp due to
alignment truncation at sequence divergence. This module refines those
boundaries by scanning a window around the approximate position and
scoring candidate boundary positions using a logistic regression trained
on TnCentral ground-truth boundaries.

Features at each candidate position p (for a start boundary):
  - 2-mer composition of 20bp window centred at p (16 features)
  - TSD score: best exact k-mer match (k=3..9) flanking p (1 feature)
  - TSD Hamming: lowest Hamming distance for any k-mer pair (1 feature)
  - IR mismatch: mismatches in 20bp IR at p vs reverse complement at end (1 feature)
  - GC delta: |GC left 10bp - GC right 10bp| across p (1 feature)
  - Homopolymer: longest run in 10bp window at p (1 feature)
Total: 21 features per position.

Usage:
  # Train:
  python3 scripts/train_boundary_refiner.py

  # Inference (called from boundaries.py or __main__.py):
  from matryoshka.boundary_refine import refine_boundary
  new_start = refine_boundary(seq, approx_start, search_window=100,
                              expected_tsd_len=8, ir_seq=feature.ir_left)
"""

from __future__ import annotations

import pickle
from pathlib import Path

import numpy as np

_MODEL_PATH = Path(__file__).parent / "models" / "boundary_model.pkl"
_model = None  # lazy-loaded


def _load_model():
    global _model
    if _model is None and _MODEL_PATH.exists():
        with open(_MODEL_PATH, "rb") as f:
            _model = pickle.load(f)
    return _model


# ---------------------------------------------------------------------------
# Feature extraction
# ---------------------------------------------------------------------------

_DIMERS = [a + b for a in "ACGT" for b in "ACGT"]
_DIMER_IDX = {d: i for i, d in enumerate(_DIMERS)}


def _dimer_freq(seq: str) -> np.ndarray:
    """16-element 2-mer frequency vector over seq."""
    vec = np.zeros(16, dtype=np.float32)
    seq = seq.upper()
    for i in range(len(seq) - 1):
        d = seq[i:i + 2]
        if d in _DIMER_IDX:
            vec[_DIMER_IDX[d]] += 1
    total = vec.sum()
    return vec / total if total > 0 else vec


def _tsd_score_paired(seq: str, start: int, end: int, k_range=(3, 9)) -> tuple[float, float]:
    """TSD score comparing seq[start-k:start] vs seq[end:end+k] (true TSD location)."""
    best_exact = 0.0
    best_hamming = float("inf")
    for k in range(k_range[0], k_range[1] + 1):
        if start < k or end + k > len(seq):
            continue
        l_kmer = seq[start - k: start].upper()
        r_kmer = seq[end: end + k].upper()
        if "N" in l_kmer or "N" in r_kmer:
            continue
        if l_kmer == r_kmer:
            best_exact = max(best_exact, k / k_range[1])
        ham = sum(a != b for a, b in zip(l_kmer, r_kmer))
        best_hamming = min(best_hamming, ham / k)
    return best_exact, best_hamming if best_hamming < float("inf") else 1.0


def _tsd_score(left_flank: str, right_flank: str, k_range=(3, 9)) -> tuple[float, float]:
    """Best exact TSD match comparing left_flank[-k:] vs right_flank[:k]."""
    left = left_flank.upper()
    right = right_flank.upper()
    best_exact = 0.0
    best_hamming = float("inf")
    for k in range(k_range[0], k_range[1] + 1):
        if len(left) < k or len(right) < k:
            continue
        l_kmer = left[-k:]
        r_kmer = right[:k]
        if "N" in l_kmer or "N" in r_kmer:
            continue
        if l_kmer == r_kmer:
            best_exact = max(best_exact, k / k_range[1])
        ham = sum(a != b for a, b in zip(l_kmer, r_kmer))
        best_hamming = min(best_hamming, ham / k)
    return best_exact, best_hamming if best_hamming < float("inf") else 1.0


def _ir_mismatch(seq: str, pos: int, ir_len: int = 20) -> float:
    """Fraction of mismatches between left IR at pos and reverse complement
    of right IR at the far end (uses last ir_len bp of seq as IRR proxy)."""
    from Bio.Seq import Seq
    if pos < 0 or pos + ir_len > len(seq) or ir_len * 2 > len(seq):
        return 1.0
    irl = seq[pos: pos + ir_len].upper()
    irr = seq[-ir_len:].upper()
    irr_rc = str(Seq(irr).reverse_complement())
    mismatches = sum(a != b for a, b in zip(irl, irr_rc, strict=True))
    return mismatches / ir_len


def _gc_delta(seq: str, pos: int, half: int = 10) -> float:
    """Absolute GC content difference across the boundary at pos."""
    left = seq[max(0, pos - half): pos].upper()
    right = seq[pos: pos + half].upper()
    if not left or not right:
        return 0.0
    gc_l = (left.count("G") + left.count("C")) / len(left)
    gc_r = (right.count("G") + right.count("C")) / len(right)
    return abs(gc_l - gc_r)


def _homopolymer(seq: str, pos: int, half: int = 10) -> float:
    """Longest homopolymer run in window as fraction of window length."""
    window = seq[max(0, pos - half): pos + half].upper()
    if not window:
        return 0.0
    max_run = best = 1
    for i in range(1, len(window)):
        if window[i] == window[i - 1]:
            best += 1
            max_run = max(max_run, best)
        else:
            best = 1
    return max_run / len(window)


def extract_features(seq: str, pos: int, ir_len: int = 20) -> np.ndarray:
    """21-element feature vector for a candidate boundary at pos."""
    half = 10
    window = seq[max(0, pos - half): pos + half]
    dimer = _dimer_freq(window)                                         # 16

    left_flank = seq[max(0, pos - 9): pos]
    right_flank = seq[pos: pos + 9]
    tsd_exact, tsd_ham = _tsd_score(left_flank, right_flank)          # 2

    ir_mm = _ir_mismatch(seq, pos, ir_len)                             # 1
    gc_d = _gc_delta(seq, pos)                                         # 1
    hp = _homopolymer(seq, pos)                                        # 1

    return np.concatenate([dimer, [tsd_exact, tsd_ham, ir_mm, gc_d, hp]])


# ---------------------------------------------------------------------------
# Inference
# ---------------------------------------------------------------------------

def _extract_paired_features(seq: str, start: int, end: int, true_len: int) -> np.ndarray:
    """38-element feature vector for a (start, end) boundary pair."""
    tsd_exact, tsd_ham = _tsd_score_paired(seq, start, end)
    ir_mm = _ir_mismatch(seq, start, 20)
    gc_s = _gc_delta(seq, start)
    gc_e = _gc_delta(seq, end)
    win_s = seq[max(0, start - 10): start + 10]
    win_e = seq[max(0, end - 10): end + 10]
    dim_s = _dimer_freq(win_s)
    dim_e = _dimer_freq(win_e)
    len_ratio = (end - start) / true_len if true_len > 0 else 1.0
    return np.concatenate([[tsd_exact, tsd_ham, ir_mm, gc_s, gc_e], dim_s, dim_e, [len_ratio]])


def refine_boundaries(
    seq: str,
    approx_start: int,
    approx_end: int,
    search_window: int = 100,
    step: int = 3,
    min_confidence: float = 0.55,
) -> tuple[int, int]:
    """Return refined (start, end) boundaries near (approx_start, approx_end).

    Scans a grid of (start, end) candidate pairs within search_window of the
    approximate positions and returns the pair with highest model confidence.
    Falls back to (approx_start, approx_end) if the model is unavailable or
    not confident enough.

    step controls scan resolution (default 3bp steps to keep runtime manageable).
    """
    model = _load_model()
    if model is None:
        return approx_start, approx_end

    true_len = approx_end - approx_start
    lo_s = max(20, approx_start - search_window)
    hi_s = min(approx_end - 50, approx_start + search_window)
    lo_e = max(approx_start + 50, approx_end - search_window)
    hi_e = min(len(seq) - 20, approx_end + search_window)

    if lo_s >= hi_s or lo_e >= hi_e:
        return approx_start, approx_end

    candidates = [
        (s, e)
        for s in range(lo_s, hi_s + 1, step)
        for e in range(lo_e, hi_e + 1, step)
        if e > s
    ]
    if not candidates:
        return approx_start, approx_end

    X = np.stack([_extract_paired_features(seq, s, e, true_len) for s, e in candidates])
    probs = model.predict_proba(X)[:, 1]
    best_idx = int(np.argmax(probs))
    if probs[best_idx] >= min_confidence:
        return candidates[best_idx]
    return approx_start, approx_end


def refine_boundary(
    seq: str,
    approx_pos: int,
    search_window: int = 100,
    min_confidence: float = 0.55,
) -> int:
    """Single-boundary refinement fallback (uses single-position features).

    Prefer refine_boundaries() when both boundaries are known.
    """
    model = _load_model()
    if model is None:
        return approx_pos

    lo = max(20, approx_pos - search_window)
    hi = min(len(seq) - 20, approx_pos + search_window)
    if lo >= hi:
        return approx_pos

    candidates = np.arange(lo, hi + 1)
    X = np.stack([extract_features(seq, int(p)) for p in candidates])
    probs = model.predict_proba(X)[:, 1]
    best_idx = int(np.argmax(probs))
    if probs[best_idx] >= min_confidence:
        return int(candidates[best_idx])
    return approx_pos
