#!/usr/bin/env python3
"""Train the ML boundary refinement model on TnCentral ground-truth boundaries.

Uses PAIRED (start, end) features so TSD signal (seq[start-k:start] vs
seq[end:end+k]) can be used — this is the dominant structural signal at
IS element boundaries and cannot be computed from a single boundary position.

For each IS element in the TnCentral ground truth:
  - True boundaries: (true_start, true_end)
  - Positive examples: offsets (ds, de) within POSITIVE_RADIUS of (0, 0)
  - Negative examples: offsets where |ds| > NEGATIVE_GAP or |de| > NEGATIVE_GAP

Features per (start_candidate, end_candidate) pair:
  - TSD score: quality of seq[start-k:start] == seq[end:end+k] (2 features)
  - IR mismatch: seq[start:start+20] vs RC of seq[end-20:end] (1 feature)
  - GC delta at start and end (2 features)
  - 2-mer composition of 20bp windows at start and end (32 features)
  - Element length ratio: candidate_length / true_length (1 feature)
Total: 38 features.

Trains a LogisticRegression and saves to matryoshka/models/boundary_model.pkl.

Usage:
    pixi run python3 scripts/train_boundary_refiner.py
"""

import csv
import pickle
import sys
from pathlib import Path

import numpy as np
from Bio import SeqIO
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import classification_report, roc_auc_score
from sklearn.model_selection import StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

sys.path.insert(0, str(Path(__file__).parent.parent))
from matryoshka.boundary_refine import _dimer_freq, _gc_delta, _ir_mismatch, _tsd_score_paired

GT_TSV = Path("data/tncentral/ground_truth.tsv")
SEQ_DIR = Path("data/tncentral/sequences")
MODEL_DIR = Path("matryoshka/models")
MODEL_PATH = MODEL_DIR / "boundary_model.pkl"

POSITIVE_RADIUS = 5     # offset within this from true boundary = positive
NEGATIVE_GAP = 20       # offset farther than this from true = negative
SEARCH_WINDOW = 80      # scan window ± around true boundary
NEG_PER_POS = 6         # negative samples per positive example


def extract_paired_features(seq: str, start: int, end: int, true_len: int) -> np.ndarray:
    """38-element feature vector for a (start, end) boundary pair."""
    tsd_exact, tsd_ham = _tsd_score_paired(seq, start, end)      # 2
    ir_mm = _ir_mismatch(seq, start, 20)                          # 1 (uses seq[-20:] proxy)
    gc_s = _gc_delta(seq, start)                                   # 1
    gc_e = _gc_delta(seq, end)                                     # 1
    win_s = seq[max(0, start - 10): start + 10]
    win_e = seq[max(0, end - 10): end + 10]
    dim_s = _dimer_freq(win_s)                                     # 16
    dim_e = _dimer_freq(win_e)                                     # 16
    cand_len = end - start
    len_ratio = cand_len / true_len if true_len > 0 else 1.0      # 1
    return np.concatenate([[tsd_exact, tsd_ham, ir_mm, gc_s, gc_e], dim_s, dim_e, [len_ratio]])


def load_ground_truth() -> dict[str, list[dict]]:
    by_seq: dict[str, list[dict]] = {}
    with open(GT_TSV) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if row["element_type"] == "IS":
                by_seq.setdefault(row["tn_name"], []).append(row)
    return by_seq


def load_sequence(tn_name: str) -> str | None:
    fa = SEQ_DIR / (tn_name + ".fasta")
    if not fa.exists():
        return None
    recs = list(SeqIO.parse(fa, "fasta"))
    return str(recs[0].seq) if recs else None


def build_dataset(by_seq: dict) -> tuple[np.ndarray, np.ndarray]:
    X_parts, y_parts = [], []
    rng = np.random.default_rng(42)
    skipped = used = 0

    for tn_name, feats in by_seq.items():
        seq = load_sequence(tn_name)
        if seq is None or len(seq) < 200:
            skipped += 1
            continue

        for feat in feats:
            ts = int(feat["start"])
            te = int(feat["end"])
            true_len = te - ts
            if ts < 20 or te > len(seq) - 20 or true_len < 50:
                continue

            # Positive: paired offsets both within POSITIVE_RADIUS
            pos_pairs = [
                (ts + ds, te + de)
                for ds in range(-POSITIVE_RADIUS, POSITIVE_RADIUS + 1)
                for de in range(-POSITIVE_RADIUS, POSITIVE_RADIUS + 1)
                if 20 <= ts + ds < te + de <= len(seq) - 20
            ]
            for s, e in pos_pairs:
                X_parts.append(extract_paired_features(seq, s, e, true_len))
                y_parts.append(1)

            # Negative: at least one boundary off by > NEGATIVE_GAP
            neg_pool = [
                (ts + ds, te + de)
                for ds in range(-SEARCH_WINDOW, SEARCH_WINDOW + 1, 3)
                for de in range(-SEARCH_WINDOW, SEARCH_WINDOW + 1, 3)
                if (abs(ds) > NEGATIVE_GAP or abs(de) > NEGATIVE_GAP)
                and 20 <= ts + ds < te + de <= len(seq) - 20
            ]
            if not neg_pool:
                continue
            n_neg = min(len(neg_pool), len(pos_pairs) * NEG_PER_POS)
            chosen = [neg_pool[i] for i in rng.choice(len(neg_pool), size=n_neg, replace=False)]
            for s, e in chosen:
                X_parts.append(extract_paired_features(seq, s, e, true_len))
                y_parts.append(0)
            used += 1

    print(f"Sequences used: {used}, skipped: {skipped}", file=sys.stderr)
    return np.stack(X_parts), np.array(y_parts)


def main():
    MODEL_DIR.mkdir(exist_ok=True)
    print("Loading ground truth...", file=sys.stderr)
    by_seq = load_ground_truth()
    print(f"Sequences with IS annotations: {len(by_seq)}", file=sys.stderr)

    print("Building dataset...", file=sys.stderr)
    X, y = build_dataset(by_seq)
    print(f"Dataset: {X.shape[0]} examples, {y.sum()} positive ({y.mean():.1%})",
          file=sys.stderr)

    model = Pipeline([
        ("scaler", StandardScaler()),
        ("clf", LogisticRegression(
            class_weight="balanced",
            max_iter=2000,
            C=0.5,
            solver="lbfgs",
        )),
    ])

    # Cross-validate to check generalisation
    print("\nCross-validation (5-fold):", file=sys.stderr)
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    aucs = []
    for fold, (tr_idx, val_idx) in enumerate(skf.split(X, y)):
        model.fit(X[tr_idx], y[tr_idx])
        probs = model.predict_proba(X[val_idx])[:, 1]
        auc = roc_auc_score(y[val_idx], probs)
        aucs.append(auc)
        print(f"  Fold {fold+1}: AUC={auc:.3f}", file=sys.stderr)
    print(f"  Mean AUC: {np.mean(aucs):.3f} ± {np.std(aucs):.3f}", file=sys.stderr)

    # Train final model on all data
    print("\nTraining final model on full dataset...", file=sys.stderr)
    model.fit(X, y)
    probs = model.predict_proba(X)[:, 1]
    preds = model.predict(X)
    print(classification_report(y, preds, target_names=["non-boundary", "boundary"]),
          file=sys.stderr)

    with open(MODEL_PATH, "wb") as f:
        pickle.dump(model, f)
    print(f"\nModel saved to {MODEL_PATH}", file=sys.stderr)

    # Feature importance (logistic regression coefficients)
    feature_names = (
        ["tsd_exact", "tsd_hamming", "ir_mismatch", "gc_delta_start", "gc_delta_end"]
        + [f"2mer_start_{a+b}" for a in "ACGT" for b in "ACGT"]
        + [f"2mer_end_{a+b}" for a in "ACGT" for b in "ACGT"]
        + ["len_ratio"]
    )
    coef = model.named_steps["clf"].coef_[0]
    print("\nTop positive features (boundary indicators):", file=sys.stderr)
    for i in np.argsort(coef)[::-1][:8]:
        print(f"  {feature_names[i]:28s}: {coef[i]:+.3f}", file=sys.stderr)
    print("Top negative features:", file=sys.stderr)
    for i in np.argsort(coef)[:5]:
        print(f"  {feature_names[i]:28s}: {coef[i]:+.3f}", file=sys.stderr)


if __name__ == "__main__":
    main()
