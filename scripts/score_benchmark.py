#!/usr/bin/env python3
"""
score_benchmark.py — Compare Matryoshka predictions against TnCentral ground truth.

Computes:
  1. Transposon-level recall: did Matryoshka detect the composite/unit transposon structure?
  2. Feature-level recall: for each annotated feature, was it detected within +-50bp?
  3. Feature-level precision: for each predicted feature, is there a ground truth match?
  4. Failure mode classification: upstream tool miss / hierarchy failure / boundary failure
  5. Summary report (Markdown)

Usage:
    pixi run python3 scripts/score_benchmark.py
    pixi run python3 scripts/score_benchmark.py --tolerance=100
"""

from __future__ import annotations

import argparse
import json
from collections import Counter, defaultdict
from pathlib import Path

PROJECT_ROOT = Path(__file__).parent.parent
DATA_DIR = PROJECT_ROOT / "data" / "tncentral"
PARSED_DIR = DATA_DIR / "parsed"
RESULTS_DIR = DATA_DIR / "benchmark_results"
REPORT_PATH = DATA_DIR / "benchmark_report.md"

# Boundary tolerance in bp for feature matching
DEFAULT_TOLERANCE = 50


def _flatten_features(features: list[dict], depth: int = 0) -> list[dict]:
    """Recursively flatten a nested feature tree."""
    flat: list[dict] = []
    for f in features:
        entry = dict(f)
        entry.pop("children", None)
        entry["_depth"] = depth
        flat.append(entry)
        for child in f.get("children", []):
            flat.extend(_flatten_features([child], depth + 1))
    return flat


def _features_overlap(
    a_start: int, a_end: int,
    b_start: int, b_end: int,
    tolerance: int,
) -> bool:
    """Check if two features match within a boundary tolerance."""
    return (
        abs(a_start - b_start) <= tolerance
        and abs(a_end - b_end) <= tolerance
    )


def _features_overlap_any(
    a_start: int, a_end: int,
    b_start: int, b_end: int,
    tolerance: int,
) -> bool:
    """Check if two features overlap at all, accounting for tolerance."""
    return (
        a_start - tolerance <= b_end
        and b_start - tolerance <= a_end
    )


def _normalise_element_type(etype: str) -> str:
    """Map TnCentral and Matryoshka element types to comparable categories."""
    etype_lower = etype.lower()
    if etype_lower in ("is", "insertion sequence"):
        return "IS"
    if etype_lower in ("transposon", "composite_transposon", "unit_transposon"):
        return "transposon"
    if etype_lower in ("integron",):
        return "integron"
    if etype_lower in ("amr",):
        return "AMR"
    if etype_lower in ("genomic_island",):
        return "genomic_island"
    if etype_lower in ("replicon",):
        return "replicon"
    if etype_lower in ("mobile_element",):
        return "mobile_element"
    return etype_lower


def score_transposon_level(
    ground_truth: dict, predictions: list[dict]
) -> dict:
    """Score whether Matryoshka detected the top-level transposon structure."""
    gt_struct = ground_truth.get("structure", {})
    gt_type = gt_struct.get("type", "")
    gt_name = gt_struct.get("name", "")
    gt_family = gt_struct.get("family", "")
    element_class = ground_truth.get("element_class", "")
    seq_len = ground_truth.get("length", 0)

    result = {
        "element_name": ground_truth["element_name"],
        "element_class": element_class,
        "gt_type": gt_type,
        "gt_name": gt_name,
        "gt_family": gt_family,
        "gt_length": seq_len,
        "detected": False,
        "detected_type": None,
        "detected_family": None,
        "detected_name": None,
        "boundary_match": False,
        "notes": [],
    }

    if not predictions:
        result["notes"].append("no_predictions")
        return result

    # Flatten all predictions
    flat_pred = _flatten_features(predictions)

    # Look for any transposon/integron/genomic_island prediction
    # that spans a substantial portion of the sequence
    for p in flat_pred:
        p_type = p.get("element_type", "")
        p_start = p.get("start", 0)
        p_end = p.get("end", 0)
        p_span = p_end - p_start

        if p_type in ("transposon", "genomic_island", "integron"):
            # Check if prediction covers >=50% of the sequence
            overlap_start = max(1, p_start)
            overlap_end = min(seq_len, p_end)
            overlap = max(0, overlap_end - overlap_start)

            if overlap >= seq_len * 0.5 or p_span >= seq_len * 0.5:
                result["detected"] = True
                result["detected_type"] = p_type
                result["detected_family"] = p.get("family", "")
                result["detected_name"] = p.get("name", "")
                # Check boundary accuracy
                if abs(p_start - 1) <= 50 and abs(p_end - seq_len) <= 50:
                    result["boundary_match"] = True
                break

    # Also check for family-name matches even without span coverage
    if not result["detected"]:
        for p in flat_pred:
            p_name = (p.get("name", "") or "").lower()
            p_family = (p.get("family", "") or "").lower()
            gt_name_lower = gt_name.lower()

            if gt_name_lower and (gt_name_lower in p_name or gt_name_lower in p_family):
                result["detected"] = True
                result["detected_type"] = p.get("element_type", "")
                result["detected_family"] = p.get("family", "")
                result["detected_name"] = p.get("name", "")
                result["notes"].append("name_match_only")
                break

    if not result["detected"]:
        # Classify why it was missed
        if element_class == "IS_element":
            result["notes"].append("IS_element_not_in_scope")
        elif element_class == "composite_transposon":
            # Check if any IS elements were detected
            has_is = any(p["element_type"] == "IS" for p in flat_pred)
            has_amr = any(p["element_type"] == "AMR" for p in flat_pred)
            if not has_is and not has_amr:
                result["notes"].append("no_component_features_detected")
            elif has_is and not has_amr:
                result["notes"].append("IS_detected_but_no_AMR")
            elif has_amr and not has_is:
                result["notes"].append("AMR_detected_but_no_IS")
            else:
                result["notes"].append("components_found_but_no_transposon_inferred")
        elif element_class == "unit_transposon":
            result["notes"].append("unit_transposon_no_blast_match")
        elif element_class == "integron":
            result["notes"].append("integron_not_detected")

    return result


def score_feature_level(
    ground_truth: dict,
    predictions: list[dict],
    tolerance: int = DEFAULT_TOLERANCE,
) -> dict:
    """Score feature-level recall and precision."""
    flat_pred = _flatten_features(predictions)

    # Ground truth features: IS elements + AMR genes + integrons
    gt_features = []
    for is_elem in ground_truth.get("is_elements", []):
        if is_elem.get("spans_full_record"):
            continue  # skip the full-record-spanning element
        gt_features.append({
            "type": "IS",
            "name": is_elem["name"],
            "start": is_elem["start"],
            "end": is_elem["end"],
            "strand": is_elem.get("strand", "."),
        })
    for amr in ground_truth.get("amr_genes", []):
        gt_features.append({
            "type": "AMR",
            "name": amr["name"],
            "start": amr["start"],
            "end": amr["end"],
            "strand": amr.get("strand", "."),
        })
    for integ in ground_truth.get("integrons", []):
        if integ.get("spans_full_record"):
            continue
        gt_features.append({
            "type": "integron",
            "name": integ["name"],
            "start": integ["start"],
            "end": integ["end"],
            "strand": integ.get("strand", "."),
        })

    # Predicted features (exclude res_site for matching)
    pred_features = []
    for p in flat_pred:
        p_type = _normalise_element_type(p.get("element_type", ""))
        if p_type in ("IS", "AMR", "integron", "transposon", "genomic_island",
                       "cassette", "integrase", "replicon"):
            pred_features.append({
                "type": p_type,
                "name": p.get("name", ""),
                "family": p.get("family", ""),
                "start": p.get("start", 0),
                "end": p.get("end", 0),
                "strand": p.get("strand", "."),
            })

    # Feature-level recall: for each GT feature, find best match in predictions
    gt_matched = [False] * len(gt_features)
    gt_match_details = []
    for i, gt in enumerate(gt_features):
        best_match = None
        best_score = float("inf")
        for p in pred_features:
            # Type must be compatible
            gt_norm = _normalise_element_type(gt["type"])
            p_norm = _normalise_element_type(p["type"])
            if gt_norm != p_norm:
                # Allow transposon to match mobile_element
                if not (gt_norm == "IS" and p_norm == "mobile_element"):
                    continue
            if _features_overlap(gt["start"], gt["end"],
                                 p["start"], p["end"], tolerance):
                score = abs(gt["start"] - p["start"]) + abs(gt["end"] - p["end"])
                if score < best_score:
                    best_score = score
                    best_match = p
        if best_match:
            gt_matched[i] = True
            gt_match_details.append({
                "gt": gt, "pred": best_match,
                "start_offset": best_match["start"] - gt["start"],
                "end_offset": best_match["end"] - gt["end"],
            })
        else:
            gt_match_details.append({"gt": gt, "pred": None})

    # Feature-level precision: for each prediction, find GT match
    pred_matched = [False] * len(pred_features)
    for j, p in enumerate(pred_features):
        for gt in gt_features:
            gt_norm = _normalise_element_type(gt["type"])
            p_norm = _normalise_element_type(p["type"])
            if gt_norm != p_norm:
                if not (gt_norm == "IS" and p_norm == "mobile_element"):
                    continue
            if _features_overlap_any(gt["start"], gt["end"],
                                      p["start"], p["end"], tolerance):
                pred_matched[j] = True
                break

    # Compute metrics
    gt_total = len(gt_features)
    gt_hit = sum(gt_matched)
    pred_total = len(pred_features)
    pred_hit = sum(pred_matched)

    return {
        "gt_total": gt_total,
        "gt_matched": gt_hit,
        "recall": gt_hit / gt_total if gt_total > 0 else None,
        "pred_total": pred_total,
        "pred_matched": pred_hit,
        "precision": pred_hit / pred_total if pred_total > 0 else None,
        "match_details": gt_match_details,
        "gt_by_type": dict(Counter(f["type"] for f in gt_features)),
        "gt_matched_by_type": dict(Counter(
            gt_features[i]["type"] for i in range(gt_total) if gt_matched[i]
        )),
        "pred_by_type": dict(Counter(p["type"] for p in pred_features)),
    }


def classify_failure(
    ground_truth: dict,
    tn_score: dict,
    feat_score: dict,
) -> list[str]:
    """Classify why detection failed for this element.

    Returns a list of failure mode labels.
    """
    modes = []
    element_class = ground_truth.get("element_class", "")

    if tn_score["detected"]:
        if not tn_score["boundary_match"]:
            modes.append("boundary_inaccuracy")
        return modes if modes else ["success"]

    # Not detected at all — classify why
    if element_class == "IS_element":
        modes.append("IS_only_element_out_of_scope")
        return modes

    # Check if BLAST should have found it
    gt_struct = ground_truth.get("structure", {})
    gt_family = gt_struct.get("family", "")
    gt_name = gt_struct.get("name", "")

    # Known Matryoshka BLAST references
    blast_covered_families = {
        "Tn3", "Tn21", "Tn4401", "Tn1546", "Tn1331", "Tn5393",
        "Tn7", "Tn552",
    }
    blast_covered_names = {
        "Tn3", "Tn2", "Tn21", "Tn1696", "Tn1721", "Tn6452",
        "Tn4401", "Tn1546", "Tn1331", "Tn5393", "Tn7", "Tn552",
        "Tn6022", "Tn6019", "Tn6172",
    }

    # Transposon rule coverage
    flanked_rule_families = {"Tn4401", "Tn1999", "Tn6330", "Tn2006", "Tn125"}
    signature_rule_families = {"Tn1546", "Tn1331"}

    if gt_name in blast_covered_names or gt_family in blast_covered_families:
        modes.append("blast_reference_miss_despite_coverage")
    elif gt_family in ("Tn3",) and gt_name not in blast_covered_names:
        modes.append("tn3_family_variant_not_in_reference_library")
    elif gt_family in ("Tn402", "Tn7"):
        modes.append(f"{gt_family}_family_no_matching_reference")
    elif element_class == "composite_transposon":
        # Check sub-features
        has_is = any(not i.get("spans_full_record") for i in ground_truth.get("is_elements", []))
        has_amr = bool(ground_truth.get("amr_genes"))
        if not has_is:
            modes.append("composite_no_IS_annotated_in_ground_truth")
        elif not has_amr:
            modes.append("composite_no_AMR_cargo")
        else:
            # IS and AMR present but no rule fired — probably not a covered family
            is_families = set(i.get("family", "") for i in ground_truth.get("is_elements", []))
            covered_is = {"IS6", "IS4", "IS1380", "IS30", "ISCR", "IS91"}
            if is_families & covered_is:
                modes.append("flanked_rule_exists_but_did_not_fire")
            else:
                modes.append(f"no_flanked_rule_for_IS_families_{','.join(is_families)}")
    elif element_class == "unit_transposon":
        modes.append("unit_transposon_not_in_reference_library")
    elif element_class == "integron":
        modes.append("integron_detection_requires_integron_finder")
    else:
        modes.append("unknown_element_type")

    return modes if modes else ["unclassified_failure"]


def generate_report(
    tn_scores: list[dict],
    feat_scores: list[dict],
    failure_modes: list[tuple[str, list[str]]],
    ground_truths: list[dict],
    tolerance: int,
) -> str:
    """Generate the Markdown benchmark report."""
    lines = []
    lines.append("# Matryoshka TnCentral Benchmark Report")
    lines.append("")
    lines.append(f"Benchmark run against {len(tn_scores)} TnCentral sequences.")
    lines.append(f"Feature matching tolerance: +/-{tolerance} bp.")
    lines.append("")

    # --- Overall transposon-level metrics ---
    lines.append("## Transposon-level detection")
    lines.append("")

    by_class: dict[str, list[dict]] = defaultdict(list)
    for ts in tn_scores:
        by_class[ts["element_class"]].append(ts)

    lines.append("| Element class | Total | Detected | Recall | Boundary match |")
    lines.append("|---|---|---|---|---|")
    total_all = 0
    detected_all = 0
    boundary_all = 0
    for cls in ["composite_transposon", "unit_transposon", "integron", "IS_element", "other"]:
        items = by_class.get(cls, [])
        if not items:
            continue
        total = len(items)
        detected = sum(1 for t in items if t["detected"])
        boundary = sum(1 for t in items if t["boundary_match"])
        recall = detected / total if total > 0 else 0
        lines.append(
            f"| {cls} | {total} | {detected} | {recall:.1%} | {boundary} |"
        )
        total_all += total
        detected_all += detected
        boundary_all += boundary

    overall_recall = detected_all / total_all if total_all > 0 else 0
    lines.append(
        f"| **Total** | **{total_all}** | **{detected_all}** | "
        f"**{overall_recall:.1%}** | **{boundary_all}** |"
    )
    lines.append("")

    # --- Feature-level metrics ---
    lines.append("## Feature-level metrics")
    lines.append("")

    # Aggregate by feature type
    type_gt_total: Counter = Counter()
    type_gt_matched: Counter = Counter()
    type_pred_total: Counter = Counter()
    type_pred_matched: Counter = Counter()

    for fs in feat_scores:
        for t, c in fs.get("gt_by_type", {}).items():
            type_gt_total[t] += c
        for t, c in fs.get("gt_matched_by_type", {}).items():
            type_gt_matched[t] += c
        for t, c in fs.get("pred_by_type", {}).items():
            type_pred_total[t] += c

    lines.append("### Recall by feature type")
    lines.append("")
    lines.append("| Feature type | Ground truth | Matched | Recall |")
    lines.append("|---|---|---|---|")
    for ftype in sorted(type_gt_total.keys()):
        gt = type_gt_total[ftype]
        matched = type_gt_matched.get(ftype, 0)
        recall = matched / gt if gt > 0 else 0
        lines.append(f"| {ftype} | {gt} | {matched} | {recall:.1%} |")
    total_gt = sum(type_gt_total.values())
    total_matched = sum(type_gt_matched.values())
    total_recall = total_matched / total_gt if total_gt > 0 else 0
    lines.append(
        f"| **Total** | **{total_gt}** | **{total_matched}** | **{total_recall:.1%}** |"
    )
    lines.append("")

    # Overall precision
    total_pred = sum(fs.get("pred_total", 0) for fs in feat_scores)
    total_pred_matched = sum(fs.get("pred_matched", 0) for fs in feat_scores)
    overall_precision = total_pred_matched / total_pred if total_pred > 0 else 0
    lines.append(f"**Overall feature-level precision**: {total_pred_matched}/{total_pred} "
                 f"({overall_precision:.1%})")
    lines.append("")

    # --- Failure mode breakdown ---
    lines.append("## Failure mode breakdown")
    lines.append("")

    mode_counts: Counter = Counter()
    for _, modes in failure_modes:
        for m in modes:
            mode_counts[m] += 1

    lines.append("| Failure mode | Count |")
    lines.append("|---|---|")
    for mode, count in mode_counts.most_common():
        lines.append(f"| {mode} | {count} |")
    lines.append("")

    # --- Top 10 failure cases ---
    lines.append("## Top 10 failure cases")
    lines.append("")
    lines.append("Transposons that were not detected, ordered by biological importance.")
    lines.append("")

    failures = [
        (ts, fm)
        for ts, (_, fm) in zip(tn_scores, failure_modes)
        if not ts["detected"] and ts["element_class"] != "IS_element"
    ]
    # Sort by: composites first, then unit transposons, then others
    priority = {"composite_transposon": 0, "unit_transposon": 1, "integron": 2, "other": 3}
    failures.sort(key=lambda x: (priority.get(x[0]["element_class"], 9), x[0]["element_name"]))

    lines.append("| # | Element | Class | Family | Length | Failure mode |")
    lines.append("|---|---|---|---|---|---|")
    for i, (ts, fm) in enumerate(failures[:10], 1):
        modes_str = "; ".join(fm)
        lines.append(
            f"| {i} | {ts['gt_name']} | {ts['element_class']} | "
            f"{ts['gt_family']} | {ts['gt_length']:,} bp | {modes_str} |"
        )
    lines.append("")

    # --- Cross-reference with GAPS.md ---
    lines.append("## Cross-reference with GAPS.md known limitations")
    lines.append("")
    lines.append("The following patterns in the failure modes align with documented gaps:")
    lines.append("")

    gap_cross_refs = [
        ("tn3_family_variant_not_in_reference_library",
         "GAPS.md: Tn3-family coverage limited to shipped references"),
        ("unit_transposon_not_in_reference_library",
         "GAPS.md: many unit transposons (Tn916, Tn554, etc.) not covered"),
        ("integron_detection_requires_integron_finder",
         "GAPS.md: integrons require IntegronFinder (not in BLAST-only mode)"),
        ("Tn402_family_no_matching_reference",
         "GAPS.md: Tn402/Tn5053 tni module detection explicitly listed as priority #1"),
        ("no_flanked_rule_for_IS_families",
         "GAPS.md: IS3, IS110, IS200, IS1, IS5, IS66 have no family-specific rules"),
        ("IS_only_element_out_of_scope",
         "Matryoshka IS detection depends on ISEScan — standalone IS detection "
         "is not tested in BLAST-only mode"),
    ]
    for mode, note in gap_cross_refs:
        if mode_counts.get(mode, 0) > 0 or any(
            mode in m for _, modes in failure_modes for m in modes
        ):
            lines.append(f"- **{mode}** ({mode_counts.get(mode, 0)} cases): {note}")
    lines.append("")

    # --- Detailed precision analysis ---
    lines.append("## Prediction analysis")
    lines.append("")
    lines.append(f"Total predictions across all sequences: {total_pred}")
    lines.append(f"Predictions with ground-truth match: {total_pred_matched} "
                 f"({overall_precision:.1%})")
    lines.append("")

    # Prediction type breakdown
    lines.append("### Predictions by type")
    lines.append("")
    lines.append("| Type | Count |")
    lines.append("|---|---|")
    pred_type_counts: Counter = Counter()
    for fs in feat_scores:
        for t, c in fs.get("pred_by_type", {}).items():
            pred_type_counts[t] += c
    for ptype, count in pred_type_counts.most_common():
        lines.append(f"| {ptype} | {count} |")
    lines.append("")

    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser(description="Score Matryoshka TnCentral benchmark")
    parser.add_argument("--tolerance", type=int, default=DEFAULT_TOLERANCE,
                        help=f"Boundary tolerance in bp (default: {DEFAULT_TOLERANCE})")
    args = parser.parse_args()

    # Load parsed ground truth
    parsed_files = sorted(PARSED_DIR.glob("*.json"))
    ground_truth_by_file: dict[str, dict] = {}
    for pf in parsed_files:
        with open(pf) as fh:
            gt = json.load(fh)
        ground_truth_by_file[gt["filename"]] = gt

    # Load benchmark results
    result_files = sorted(RESULTS_DIR.glob("*.json"))
    if not result_files:
        print("No benchmark results found. Run benchmark_tncentral.py first.")
        return

    print(f"Scoring {len(result_files)} benchmark results...")

    tn_scores = []
    feat_scores = []
    failure_modes = []
    ground_truths = []

    for rf in result_files:
        with open(rf) as fh:
            result = json.load(fh)

        filename = result.get("filename", "")
        gt = ground_truth_by_file.get(filename)
        if gt is None:
            continue

        predictions = result.get("predictions", [])
        if result.get("error"):
            predictions = []

        # Transposon-level score
        ts = score_transposon_level(gt, predictions)
        tn_scores.append(ts)

        # Feature-level score
        fs = score_feature_level(gt, predictions, tolerance=args.tolerance)
        feat_scores.append(fs)

        # Failure mode
        fm = classify_failure(gt, ts, fs)
        failure_modes.append((filename, fm))
        ground_truths.append(gt)

    # Generate report
    report = generate_report(
        tn_scores, feat_scores, failure_modes, ground_truths, args.tolerance
    )
    with open(REPORT_PATH, "w") as fh:
        fh.write(report)
    print(f"Report written to {REPORT_PATH}")

    # Print summary to stdout
    detected = sum(1 for ts in tn_scores if ts["detected"])
    total = len(tn_scores)
    print(f"\nTransposon-level recall: {detected}/{total} ({detected/total:.1%})")

    total_gt = sum(fs["gt_total"] for fs in feat_scores)
    total_matched = sum(fs["gt_matched"] for fs in feat_scores)
    if total_gt > 0:
        print(f"Feature-level recall: {total_matched}/{total_gt} ({total_matched/total_gt:.1%})")

    total_pred = sum(fs["pred_total"] for fs in feat_scores)
    total_pred_m = sum(fs["pred_matched"] for fs in feat_scores)
    if total_pred > 0:
        print(f"Feature-level precision: {total_pred_m}/{total_pred} "
              f"({total_pred_m/total_pred:.1%})")

    # Top failure modes
    mode_counts: Counter = Counter()
    for _, modes in failure_modes:
        for m in modes:
            mode_counts[m] += 1
    print("\nTop failure modes:")
    for mode, count in mode_counts.most_common(5):
        print(f"  {mode}: {count}")


if __name__ == "__main__":
    main()
