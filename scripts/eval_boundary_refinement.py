#!/usr/bin/env python3
"""
eval_boundary_refinement.py — Evaluate the ML boundary refinement model.

For each TnCentral benchmark element that Matryoshka correctly identifies,
compares boundaries with and without ML refinement against ground truth.

Measures:
  - Mean absolute boundary error (MAE) before/after ML refinement
  - Number of elements where ML moved boundary closer/further
  - Net improvement (MAE_before - MAE_after)

Usage:
    pixi run python3 scripts/eval_boundary_refinement.py
"""

from __future__ import annotations

import copy
import json
import sys
import tempfile
import warnings
from collections import defaultdict
from pathlib import Path

from Bio import SeqIO

PROJECT_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from matryoshka.boundaries import confirm_boundaries
from matryoshka.detect import MGEFeature

DATA_DIR = PROJECT_ROOT / "data" / "tncentral"
RESULTS_DIR = DATA_DIR / "benchmark_results"
PARSED_DIR = DATA_DIR / "parsed"
GB_DIR = DATA_DIR / "gb"

# Benchmark representative set — same as run_benchmark.py
REPRESENTATIVE_SET: list[str] = [
    "Tn1412-L36547",
    "Tn2000-AF205943",
    "Tn6023-GU562437.2",
    "Tn6010-EU370913",
    "Tn6309-KX710094",
    "TnAph1-L36457",
    "TnPMLUA4-KC964607.1",
    "Tn2012-EU523120",
    "TnEcp1.1-MF062700",
    "Tn2.1-CP028717",
    "Tn125-JN872328",
    "Tn2006-EF127491.1",
    "Tn3-V00613",
    "Tn2-KT002541",
    "Tn21-AF071413",
    "Tn1696-U12338.3",
    "Tn1721-X61367.1",
    "Tn1331-KC354802.1",
    "Tn1546-M97297.1",
    "Tn5393-M95402.1",
    "Tn4401a-KT378596",
    "Tn4401b-KT378597",
    "In2-AF071413",
    "In0-U49101",
    "In1-AY046276",
    "In104-AY463797",
    "In336-KP873172",
    "Tn10-AF162223",
    "Tn5-U00004.1",
    "Tn4001-GU235985",
]


def _load_sequence(stem: str) -> str | None:
    """Load the genomic sequence from a GenBank file."""
    gb_path = GB_DIR / f"{stem}.gb"
    if not gb_path.exists():
        candidates = list(GB_DIR.glob(f"{stem}*.gb"))
        if candidates:
            gb_path = candidates[0]
        else:
            return None
    try:
        warnings.filterwarnings("ignore", category=BiopythonWarning)
    except NameError:
        pass
    warnings.filterwarnings("ignore")
    rec = next(SeqIO.parse(str(gb_path), "genbank"))
    return str(rec.seq)


def _load_benchmark_result(stem: str) -> dict | None:
    """Load a benchmark result JSON."""
    result_path = RESULTS_DIR / f"{stem}.json"
    if not result_path.exists():
        return None
    with open(result_path) as fh:
        return json.load(fh)


def _load_ground_truth(stem: str) -> dict | None:
    """Load parsed TnCentral ground truth."""
    json_path = PARSED_DIR / f"{stem}.json"
    if not json_path.exists():
        return None
    with open(json_path) as fh:
        return json.load(fh)


def _dict_to_feature(d: dict) -> MGEFeature:
    """Reconstruct an MGEFeature from a serialised dict."""
    attrs = d.get("attributes", {})
    # Convert string values back to native types where appropriate
    clean_attrs = {}
    for k, v in attrs.items():
        clean_attrs[k] = v
    f = MGEFeature(
        element_type=d["element_type"],
        family=d["family"],
        name=d["name"],
        start=d["start"],
        end=d["end"],
        strand=d.get("strand", "."),
        tsd_seq=d.get("tsd_seq"),
        ir_left=d.get("ir_left"),
        ir_right=d.get("ir_right"),
        score=d.get("score"),
        attributes=clean_attrs,
    )
    return f


def _get_gt_features(gt: dict) -> list[dict]:
    """Extract ground truth features with boundaries for comparison."""
    features = []
    # IS elements (not spanning full record)
    for ise in gt.get("is_elements", []):
        if not ise.get("spans_full_record"):
            features.append({
                "type": "IS",
                "name": ise["name"],
                "start": ise["start"],
                "end": ise["end"],
            })
    # Transposons (not spanning full record)
    for me in gt.get("mobile_elements", []):
        if not me.get("spans_full_record"):
            features.append({
                "type": "transposon",
                "name": me["name"],
                "start": me["start"],
                "end": me["end"],
            })
    # Integrons (not spanning full record)
    for integ in gt.get("integrons", []):
        if not integ.get("spans_full_record"):
            features.append({
                "type": "integron",
                "name": integ["name"],
                "start": integ["start"],
                "end": integ["end"],
            })
    return features


def _match_to_gt(
    pred_feature: MGEFeature,
    gt_features: list[dict],
    tolerance: int = 200,
) -> dict | None:
    """Find the closest ground-truth feature matching a prediction.

    Uses a generous tolerance to find matches even when boundaries are
    significantly off — we want to measure the error, not filter by it.
    """
    best_match = None
    best_distance = float("inf")

    for gf in gt_features:
        # Type compatibility
        if pred_feature.element_type == "IS" and gf["type"] != "IS":
            continue
        if pred_feature.element_type == "transposon" and gf["type"] not in ("transposon",):
            continue

        # Distance: sum of absolute boundary differences
        start_diff = abs(pred_feature.start - gf["start"])
        end_diff = abs(pred_feature.end - gf["end"])
        total_diff = start_diff + end_diff

        if total_diff < best_distance and start_diff <= tolerance and end_diff <= tolerance:
            best_distance = total_diff
            best_match = gf

    return best_match


def evaluate_element(
    stem: str,
) -> list[dict] | None:
    """Evaluate ML boundary refinement for a single benchmark element.

    Returns a list of per-feature evaluation records, or None if the
    element cannot be evaluated.
    """
    result = _load_benchmark_result(stem)
    if result is None or not result.get("detected"):
        return None

    gt = _load_ground_truth(stem)
    if gt is None:
        return None

    seq = _load_sequence(stem)
    if seq is None:
        return None

    gt_features = _get_gt_features(gt)
    if not gt_features:
        return None

    # Extract BLAST-sourced features from the benchmark result
    all_features = result.get("all_features", [])
    blast_features = [
        f for f in all_features
        if f.get("attributes", {}).get("source") == "reference_scan"
        and f.get("element_type") in ("IS", "transposon")
        and f.get("start", 0) > 0
        and f.get("end", 0) > 0
        and (f.get("end", 0) - f.get("start", 0)) > 100
    ]

    if not blast_features:
        return None

    records = []
    for bf_dict in blast_features:
        # Create two copies: one for without ML, one for with ML
        feat_no_ml = _dict_to_feature(bf_dict)
        feat_ml = _dict_to_feature(bf_dict)

        # Record original boundaries (before any confirm_boundaries call)
        orig_start = feat_no_ml.start
        orig_end = feat_no_ml.end

        # Run confirm_boundaries WITHOUT ML refinement
        confirm_boundaries(seq, [feat_no_ml], ml_refine=False)
        no_ml_start = feat_no_ml.start
        no_ml_end = feat_no_ml.end

        # Run confirm_boundaries WITH ML refinement
        confirm_boundaries(seq, [feat_ml], ml_refine=True)
        ml_start = feat_ml.start
        ml_end = feat_ml.end

        # Find the ground truth match
        gt_match = _match_to_gt(feat_no_ml, gt_features)
        if gt_match is None:
            # Try matching with the ML-refined version
            gt_match = _match_to_gt(feat_ml, gt_features)
        if gt_match is None:
            continue

        gt_start = gt_match["start"]
        gt_end = gt_match["end"]

        # Compute boundary errors
        no_ml_start_err = abs(no_ml_start - gt_start)
        no_ml_end_err = abs(no_ml_end - gt_end)
        ml_start_err = abs(ml_start - gt_start)
        ml_end_err = abs(ml_end - gt_end)

        no_ml_mae = (no_ml_start_err + no_ml_end_err) / 2.0
        ml_mae = (ml_start_err + ml_end_err) / 2.0

        # Determine if ML was applied
        ml_changed = (ml_start != no_ml_start or ml_end != no_ml_end)

        records.append({
            "element": stem.split("-")[0],
            "feature_type": bf_dict["element_type"],
            "feature_name": bf_dict["name"],
            "feature_family": bf_dict["family"],
            "gt_start": gt_start,
            "gt_end": gt_end,
            "orig_start": orig_start,
            "orig_end": orig_end,
            "no_ml_start": no_ml_start,
            "no_ml_end": no_ml_end,
            "ml_start": ml_start,
            "ml_end": ml_end,
            "no_ml_start_err": no_ml_start_err,
            "no_ml_end_err": no_ml_end_err,
            "no_ml_mae": no_ml_mae,
            "ml_start_err": ml_start_err,
            "ml_end_err": ml_end_err,
            "ml_mae": ml_mae,
            "ml_changed": ml_changed,
            "ml_improved": ml_mae < no_ml_mae if ml_changed else False,
            "ml_worsened": ml_mae > no_ml_mae if ml_changed else False,
        })

    return records if records else None


def main() -> None:
    print("ML Boundary Refinement Evaluation")
    print("=" * 60)
    print()

    all_records: list[dict] = []
    elements_evaluated = 0
    elements_skipped = 0

    for stem in REPRESENTATIVE_SET:
        records = evaluate_element(stem)
        if records is None:
            elements_skipped += 1
            continue
        elements_evaluated += 1
        all_records.extend(records)

    if not all_records:
        print("No BLAST-sourced features found for evaluation.")
        return

    # Compute summary statistics
    n_total = len(all_records)
    n_changed = sum(1 for r in all_records if r["ml_changed"])
    n_improved = sum(1 for r in all_records if r["ml_improved"])
    n_worsened = sum(1 for r in all_records if r["ml_worsened"])
    n_unchanged = sum(1 for r in all_records if not r["ml_changed"])

    mae_before_values = [r["no_ml_mae"] for r in all_records]
    mae_after_values = [r["ml_mae"] for r in all_records]
    mae_before = sum(mae_before_values) / n_total
    mae_after = sum(mae_after_values) / n_total
    net_improvement = mae_before - mae_after

    # Compute median
    mae_before_sorted = sorted(mae_before_values)
    mae_after_sorted = sorted(mae_after_values)
    median_before = mae_before_sorted[n_total // 2]
    median_after = mae_after_sorted[n_total // 2]

    # Print summary
    print(f"Elements evaluated: {elements_evaluated}")
    print(f"Elements skipped (not detected or no BLAST features): {elements_skipped}")
    print(f"Total BLAST-sourced features evaluated: {n_total}")
    print()
    print("Summary")
    print("-" * 60)
    print(f"  MAE before ML refinement:  {mae_before:.1f} bp (median {median_before:.1f} bp)")
    print(f"  MAE after ML refinement:   {mae_after:.1f} bp (median {median_after:.1f} bp)")
    print(f"  Net improvement:           {net_improvement:+.1f} bp")
    print()
    print(f"  Features where ML changed boundary: {n_changed}/{n_total}")
    print(f"  Features where ML moved CLOSER:     {n_improved}/{n_total}")
    print(f"  Features where ML moved FURTHER:    {n_worsened}/{n_total}")
    print(f"  Features where ML did NOT change:   {n_unchanged}/{n_total}")
    print()

    # Per-feature detail table
    print("Per-feature details")
    print("-" * 60)
    header = (
        f"{'Element':<12} {'Feature':<15} {'GT span':<15} "
        f"{'MAE before':>10} {'MAE after':>10} {'Delta':>8} {'ML?':<5}"
    )
    print(header)
    print("-" * len(header))

    for r in sorted(all_records, key=lambda x: (x["element"], x["feature_name"])):
        gt_span = f"{r['gt_start']}-{r['gt_end']}"
        delta = r["no_ml_mae"] - r["ml_mae"]
        ml_flag = ""
        if r["ml_improved"]:
            ml_flag = "better"
        elif r["ml_worsened"]:
            ml_flag = "worse"
        elif r["ml_changed"]:
            ml_flag = "same"
        else:
            ml_flag = "n/a"
        print(
            f"{r['element']:<12} {r['feature_name']:<15} {gt_span:<15} "
            f"{r['no_ml_mae']:>10.1f} {r['ml_mae']:>10.1f} {delta:>+8.1f} {ml_flag:<5}"
        )
    print()

    # Breakdown by element type
    by_type: dict[str, list[dict]] = defaultdict(list)
    for r in all_records:
        by_type[r["feature_type"]].append(r)

    print("Breakdown by feature type")
    print("-" * 60)
    for ftype, recs in sorted(by_type.items()):
        n = len(recs)
        mae_b = sum(r["no_ml_mae"] for r in recs) / n
        mae_a = sum(r["ml_mae"] for r in recs) / n
        improved = sum(1 for r in recs if r["ml_improved"])
        worsened = sum(1 for r in recs if r["ml_worsened"])
        print(
            f"  {ftype}: {n} features, "
            f"MAE {mae_b:.1f} -> {mae_a:.1f} bp, "
            f"{improved} improved, {worsened} worsened"
        )
    print()

    # Save results to markdown
    out_path = DATA_DIR / "boundary_refinement_eval.md"
    lines = []
    lines.append("# ML Boundary Refinement Evaluation")
    lines.append("")
    lines.append(f"Elements evaluated: {elements_evaluated}")
    lines.append(f"Total BLAST-sourced features: {n_total}")
    lines.append("")
    lines.append("## Summary")
    lines.append("")
    lines.append("| Metric | Value |")
    lines.append("|---|---|")
    lines.append(f"| MAE before ML refinement | {mae_before:.1f} bp |")
    lines.append(f"| MAE after ML refinement | {mae_after:.1f} bp |")
    lines.append(f"| Median MAE before | {median_before:.1f} bp |")
    lines.append(f"| Median MAE after | {median_after:.1f} bp |")
    lines.append(f"| Net improvement | {net_improvement:+.1f} bp |")
    lines.append(f"| Features where ML changed boundary | {n_changed}/{n_total} |")
    lines.append(f"| Features where ML moved closer | {n_improved}/{n_total} |")
    lines.append(f"| Features where ML moved further | {n_worsened}/{n_total} |")
    lines.append(f"| Features unchanged | {n_unchanged}/{n_total} |")
    lines.append("")
    lines.append("## Per-feature results")
    lines.append("")
    lines.append("| Element | Feature | GT span | MAE before | MAE after | Delta | ML effect |")
    lines.append("|---|---|---|---|---|---|---|")
    for r in sorted(all_records, key=lambda x: (x["element"], x["feature_name"])):
        gt_span = f"{r['gt_start']}-{r['gt_end']}"
        delta = r["no_ml_mae"] - r["ml_mae"]
        ml_flag = ""
        if r["ml_improved"]:
            ml_flag = "improved"
        elif r["ml_worsened"]:
            ml_flag = "worsened"
        elif r["ml_changed"]:
            ml_flag = "neutral"
        else:
            ml_flag = "unchanged"
        lines.append(
            f"| {r['element']} | {r['feature_name']} | {gt_span} | "
            f"{r['no_ml_mae']:.1f} | {r['ml_mae']:.1f} | {delta:+.1f} | {ml_flag} |"
        )
    lines.append("")
    lines.append("## Breakdown by feature type")
    lines.append("")
    lines.append("| Type | Count | MAE before | MAE after | Improved | Worsened |")
    lines.append("|---|---|---|---|---|---|")
    for ftype, recs in sorted(by_type.items()):
        n = len(recs)
        mae_b = sum(r["no_ml_mae"] for r in recs) / n
        mae_a = sum(r["ml_mae"] for r in recs) / n
        improved = sum(1 for r in recs if r["ml_improved"])
        worsened = sum(1 for r in recs if r["ml_worsened"])
        lines.append(f"| {ftype} | {n} | {mae_b:.1f} | {mae_a:.1f} | {improved} | {worsened} |")
    lines.append("")

    with open(out_path, "w") as fh:
        fh.write("\n".join(lines))
    print(f"Results saved to: {out_path}")


if __name__ == "__main__":
    main()
