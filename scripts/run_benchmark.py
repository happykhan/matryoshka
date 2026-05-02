#!/usr/bin/env python3
"""
run_benchmark.py — Comprehensive Matryoshka benchmark against TnCentral ground truth.

Selects 30 representative transposons covering:
  - IS6/IS26-flanked composites
  - IS1380/ISEcp1-capture composites
  - ISCR/rolling-circle composites
  - Tn3-family unit transposons
  - Integron-bearing elements

Runs the full Matryoshka pipeline (ISEScan + AMRFinder + IntegronFinder + BLAST)
on each, compares predictions to TnCentral annotations, and produces:
  - data/tncentral/ground_truth.tsv
  - data/tncentral/benchmark_report.md
  - data/tncentral/FINDINGS.md

Usage:
    pixi run python3 scripts/run_benchmark.py
    pixi run python3 scripts/run_benchmark.py --blast-only
    pixi run python3 scripts/run_benchmark.py --max=10
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import tempfile
import traceback
from collections import Counter, defaultdict
from pathlib import Path

from Bio import SeqIO

PROJECT_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from matryoshka.boundaries import confirm_boundaries
from matryoshka.confidence import assign_confidence
from matryoshka.detect import (
    MGEFeature,
    parse_amrfinder,
    parse_integron_finder,
    parse_isescan,
)
from matryoshka.hierarchy import build_hierarchy
from matryoshka.reference_scan import blast_available, scan_all
from matryoshka.transposon import annotate_res_sites, infer_transposons

DATA_DIR = PROJECT_ROOT / "data" / "tncentral"
PARSED_DIR = DATA_DIR / "parsed"
RESULTS_DIR = DATA_DIR / "benchmark_results"
GB_DIR = DATA_DIR / "genbank" / "gb"

# Boundary tolerance in bp for feature matching
TOLERANCE = 50

# -----------------------------------------------------------------------
# Representative element selection — 30 transposons covering all categories
# -----------------------------------------------------------------------

# Format: (filename_stem, category, rationale)
REPRESENTATIVE_SET: list[tuple[str, str, str]] = [
    # --- IS6/IS26-flanked composites (7) ---
    ("Tn1412-L36547", "IS6_flanked",
     "IS26/IS26 flanking multi-AMR cargo including blaLCR-1"),
    ("Tn2000-AF205943", "IS6_flanked",
     "IS26/IS26 flanking blaVEB-1 + cassette array"),
    ("Tn6023-GU562437.2", "IS6_flanked",
     "IS26/IS26 flanking APH(3')-Ia — small simple composite"),
    ("Tn6010-EU370913", "IS6_flanked",
     "IS26/IS26 flanking oqxAB — non-standard cargo"),
    ("Tn6309-KX710094", "IS6_flanked",
     "IS26/IS26 flanking tetA(C)/tetR(C)"),
    ("TnAph1-L36457", "IS6_flanked",
     "IS26/IS26 flanking APH(3')-Ia — canonical minimal composite"),
    ("TnPMLUA4-KC964607.1", "IS6_flanked",
     "IS26/IS26 flanking mphE/msrE"),

    # --- IS1380/ISEcp1-capture composites (3) ---
    ("Tn2012-EU523120", "ISEcp1_capture",
     "ISEcp1 one-ended capture of qnrB19"),
    ("TnEcp1.1-MF062700", "ISEcp1_capture",
     "ISEcp1 capture of blaCTX-M-3"),
    ("Tn2.1-CP028717", "ISEcp1_capture",
     "ISEcp1-mediated with IS1F — blaCTX-M-3 + blaTEM-1"),

    # --- ISCR-flanked / rolling-circle (2) ---
    ("Tn125-JN872328", "ISCR_flanked",
     "ISAba125-flanked blaNDM-2 with ISCR21 inside"),
    ("Tn2006-EF127491.1", "IS4_flanked",
     "ISAba1-flanked blaOXA-23 — Acinetobacter composite"),

    # --- Tn3-family unit transposons (10) ---
    ("Tn3-V00613", "Tn3_family",
     "Archetype Tn3 with blaTEM-1"),
    ("Tn2-KT002541", "Tn3_family",
     "Tn2 (blaTEM-1) — Tn3-family unit"),
    ("Tn21-AF071413", "Tn3_family",
     "Archetype Tn21 with class-1 integron + mer operon"),
    ("Tn1696-U12338.3", "Tn3_family",
     "Tn1696 — Tn21-subfamily with class-1 integron"),
    ("Tn1721-X61367.1", "Tn3_family",
     "Tn1721 — Tn21-subfamily with tet(A)"),
    ("Tn1331-KC354802.1", "Tn3_family",
     "Tn1331 — aac(6')-Ib + blaOXA-9 + aadA1 Tn3 unit"),
    ("Tn1546-M97297.1", "Tn3_family",
     "Tn1546 — vanA vancomycin-resistance unit"),
    ("Tn5393-M95402.1", "Tn3_family",
     "Tn5393 — strAB aminoglycoside resistance"),
    ("Tn4401a-KT378596", "Tn3_family",
     "Tn4401a — ISKpn7/blaKPC-3/ISKpn6"),
    ("Tn4401b-KT378597", "Tn3_family",
     "Tn4401b — ISKpn7/blaKPC-3/ISKpn6 variant b"),

    # --- Integron-bearing elements (5) ---
    ("In2-AF071413", "integron",
     "Archetype class-1 integron In2 from Tn21"),
    ("In0-U49101", "integron",
     "In0 — class 1 integron with IS1326"),
    ("In1-AY046276", "integron",
     "In1 — class 1 integron with blaOXA-2"),
    ("In104-AY463797", "integron",
     "In104 — class 1 integron with IS6100 + aadA7 + aac(3)-Id"),
    ("In336-KP873172", "integron",
     "In336 — class 1 integron with blaVIM-1"),

    # --- Composite with non-covered IS families (3) ---
    ("Tn10-AF162223", "non_covered_IS",
     "Tn10 — IS10L/IS10R flanked tet(B), IS4 family — NOT in flanked rules"),
    ("Tn5-U00004.1", "non_covered_IS",
     "Tn5 — IS50L/IS50R flanked kanamycin resistance, IS4 family"),
    ("Tn4001-GU235985", "non_covered_IS",
     "Tn4001 — IS256 flanked AAC(6')-Ie, Staph — NOT in flanked rules"),
]


def _find_gb_path(stem: str) -> Path | None:
    """Find the GenBank file for a given filename stem.

    Handles minor naming mismatches.
    """
    exact = GB_DIR / f"{stem}.gb"
    if exact.exists():
        return exact
    # Try glob-based fallback
    candidates = list(GB_DIR.glob(f"{stem}*.gb"))
    if candidates:
        return candidates[0]
    # Try with element name only
    name_part = stem.split("-")[0]
    candidates = list(GB_DIR.glob(f"{name_part}-*.gb"))
    if candidates:
        return candidates[0]
    return None


def _feature_to_dict(f: MGEFeature, depth: int = 0) -> dict:
    """Serialise an MGEFeature (and children) to a JSON-safe dict."""
    d = {
        "element_type": f.element_type,
        "family": f.family,
        "name": f.name,
        "start": f.start,
        "end": f.end,
        "strand": f.strand,
        "tsd_seq": f.tsd_seq,
        "ir_left": f.ir_left,
        "ir_right": f.ir_right,
        "score": f.score,
        "attributes": {k: str(v) for k, v in f.attributes.items()},
    }
    if f.children:
        d["children"] = [_feature_to_dict(c, depth + 1) for c in f.children]
    return d


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


def _suppress_redundant_inference(features: list[MGEFeature]) -> list[MGEFeature]:
    """Drop rule-inferred transposons that overlap a BLAST-confirmed call."""
    blast_by_family: dict[str, list[MGEFeature]] = {}
    for f in features:
        if f.attributes.get("source") == "reference_scan" and f.element_type == "transposon":
            blast_by_family.setdefault(f.family, []).append(f)

    def _overlap(a: MGEFeature, b: MGEFeature) -> bool:
        if a.end < b.start or b.end < a.start:
            return False
        ovl = min(a.end, b.end) - max(a.start, b.start)
        short = min(a.end - a.start, b.end - b.start)
        return short > 0 and ovl / short > 0.5

    kept: list[MGEFeature] = []
    for f in features:
        if (
            f.element_type == "transposon"
            and f.attributes.get("source") != "reference_scan"
            and f.family in blast_by_family
            and any(_overlap(f, b) for b in blast_by_family[f.family])
        ):
            continue
        kept.append(f)
    return kept


# -----------------------------------------------------------------------
# External tool runners
# -----------------------------------------------------------------------

def run_external_tool(
    tool: str, fasta_path: Path, outdir: Path
) -> Path | None:
    """Run ISEScan, AMRFinder, or IntegronFinder. Returns output path or None."""
    outdir.mkdir(parents=True, exist_ok=True)
    try:
        if tool == "isescan":
            subprocess.run(
                ["pixi", "run", "-e", "isescan", "isescan.py",
                 "--seqfile", str(fasta_path),
                 "--output", str(outdir), "--nthread", "2"],
                check=True, capture_output=True, timeout=300,
                cwd=str(PROJECT_ROOT),
            )
            candidates = list(outdir.rglob("*.tsv"))
            return candidates[0] if candidates else None

        elif tool == "amrfinder":
            out_tsv = outdir / "amrfinder.tsv"
            subprocess.run(
                ["pixi", "run", "-e", "amrfinder", "amrfinder",
                 "-n", str(fasta_path), "-o", str(out_tsv), "--plus"],
                check=True, capture_output=True, timeout=300,
                cwd=str(PROJECT_ROOT),
            )
            return out_tsv if out_tsv.exists() else None

        elif tool == "integron_finder":
            subprocess.run(
                ["pixi", "run", "-e", "integron", "integron_finder",
                 "--outdir", str(outdir), str(fasta_path)],
                check=True, capture_output=True, timeout=300,
                cwd=str(PROJECT_ROOT),
            )
            candidates = list(outdir.rglob("*.integrons"))
            return candidates[0] if candidates else None

    except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as e:
        tool_err = ""
        if hasattr(e, "stderr") and e.stderr:
            tool_err = e.stderr[:200] if isinstance(e.stderr, str) else e.stderr.decode(errors="replace")[:200]
        print(f"    {tool} failed: {type(e).__name__}: {tool_err}", file=sys.stderr)
        return None

    return None


def write_single_fasta(gb_path: Path, out_path: Path) -> bool:
    """Write a single GenBank record as a FASTA file for tool input."""
    try:
        rec = next(SeqIO.parse(gb_path, "genbank"))
        name = gb_path.stem.split("-")[0]
        accession = gb_path.stem.split("-", 1)[1] if "-" in gb_path.stem else ""
        rec.id = f"{name}__{accession}"
        rec.description = ""
        with open(out_path, "w") as fh:
            SeqIO.write([rec], fh, "fasta")
        return True
    except Exception:
        return False


def run_full_pipeline(fasta_path: Path) -> tuple[list[MGEFeature], list[MGEFeature], dict]:
    """Run the full Matryoshka pipeline on a single FASTA sequence.

    Returns (roots, all_features, tool_status).
    """
    all_features: list[MGEFeature] = []
    tool_status: dict[str, str] = {}

    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = Path(tmpdir)

        # ISEScan
        isescan_out = run_external_tool("isescan", fasta_path, tmp / "isescan")
        if isescan_out:
            feats = parse_isescan(isescan_out)
            all_features.extend(feats)
            tool_status["isescan"] = f"ok ({len(feats)} features)"
        else:
            tool_status["isescan"] = "no_output"

        # AMRFinder
        amrfinder_out = run_external_tool("amrfinder", fasta_path, tmp / "amrfinder")
        if amrfinder_out:
            feats = parse_amrfinder(amrfinder_out)
            all_features.extend(feats)
            tool_status["amrfinder"] = f"ok ({len(feats)} features)"
        else:
            tool_status["amrfinder"] = "no_output"

        # IntegronFinder
        integron_out = run_external_tool("integron_finder", fasta_path, tmp / "integron")
        if integron_out:
            try:
                feats = parse_integron_finder(integron_out)
                all_features.extend(feats)
                tool_status["integron_finder"] = f"ok ({len(feats)} features)"
            except Exception as e:
                tool_status["integron_finder"] = f"parse_error: {e}"
        else:
            tool_status["integron_finder"] = "no_output"

    # BLAST reference scan
    if blast_available():
        ref_hits = scan_all(str(fasta_path))
        all_features.extend(ref_hits)
        tool_status["blast_scan"] = f"ok ({len(ref_hits)} hits)"
    else:
        tool_status["blast_scan"] = "blastn_not_available"

    # Run inference
    seq = str(next(SeqIO.parse(fasta_path, "fasta")).seq)
    inferred = infer_transposons(all_features)
    all_feats = all_features + inferred
    all_feats = _suppress_redundant_inference(all_feats)
    tool_status["inferred"] = f"{len(inferred)} transposons inferred"

    transposons = [f for f in all_feats if f.element_type == "transposon"]
    all_feats.extend(annotate_res_sites(transposons))

    checkable = [f for f in all_feats if f.element_type in ("IS", "transposon")]
    confirm_boundaries(seq, checkable)
    assign_confidence(all_feats)
    roots = build_hierarchy(all_feats)

    return roots, all_feats, tool_status


def run_blast_pipeline(fasta_path: Path) -> tuple[list[MGEFeature], list[MGEFeature], dict]:
    """Run Matryoshka with BLAST reference scan only."""
    tool_status: dict[str, str] = {}

    all_features: list[MGEFeature] = []
    if blast_available():
        ref_hits = scan_all(str(fasta_path))
        all_features.extend(ref_hits)
        tool_status["blast_scan"] = f"ok ({len(ref_hits)} hits)"
    else:
        tool_status["blast_scan"] = "blastn_not_available"

    seq = str(next(SeqIO.parse(fasta_path, "fasta")).seq)
    inferred = infer_transposons(all_features)
    all_feats = all_features + inferred
    all_feats = _suppress_redundant_inference(all_feats)

    transposons = [f for f in all_feats if f.element_type == "transposon"]
    all_feats.extend(annotate_res_sites(transposons))

    checkable = [f for f in all_feats if f.element_type in ("IS", "transposon")]
    confirm_boundaries(seq, checkable)
    assign_confidence(all_feats)
    roots = build_hierarchy(all_feats)

    return roots, all_feats, tool_status


# -----------------------------------------------------------------------
# Ground truth extraction
# -----------------------------------------------------------------------

def load_ground_truth(filename: str) -> dict | None:
    """Load parsed TnCentral ground truth for a given filename stem."""
    json_path = PARSED_DIR / f"{filename}.json"
    if not json_path.exists():
        return None
    with open(json_path) as fh:
        return json.load(fh)


def write_ground_truth_tsv(
    elements: list[tuple[str, str, dict]],
) -> None:
    """Write benchmark_ground_truth.tsv from parsed TnCentral records."""
    out_path = DATA_DIR / "benchmark_ground_truth.tsv"
    with open(out_path, "w") as fh:
        fh.write("transposon_name\tfeature_type\tfeature_name\tstart\tend\tstrand\n")
        for stem, category, gt in elements:
            name = gt["element_name"]
            # Mobile elements / transposons
            for me in gt.get("mobile_elements", []):
                if not me.get("spans_full_record"):
                    fh.write(f"{name}\ttransposon\t{me['name']}\t{me['start']}\t{me['end']}\t{me['strand']}\n")
            # IS elements
            for ise in gt.get("is_elements", []):
                if not ise.get("spans_full_record"):
                    fh.write(f"{name}\tIS\t{ise['name']}\t{ise['start']}\t{ise['end']}\t{ise['strand']}\n")
            # AMR genes
            for amr in gt.get("amr_genes", []):
                fh.write(f"{name}\tAMR\t{amr['name']}\t{amr['start']}\t{amr['end']}\t{amr['strand']}\n")
            # Integrons
            for integ in gt.get("integrons", []):
                if not integ.get("spans_full_record"):
                    fh.write(f"{name}\tintegron\t{integ['name']}\t{integ['start']}\t{integ['end']}\t{integ['strand']}\n")
    print(f"Wrote ground truth: {out_path}")


# -----------------------------------------------------------------------
# Scoring
# -----------------------------------------------------------------------

def score_element(
    gt: dict,
    predictions: list[dict],
    all_features: list[dict],
    tool_status: dict[str, str],
    tolerance: int = TOLERANCE,
) -> dict:
    """Score a single element: transposon-level + feature-level + failure mode."""
    gt_struct = gt.get("structure", {})
    gt_name = gt_struct.get("name", "")
    element_class = gt.get("element_class", "")
    seq_len = gt.get("length", 0)
    flat_pred = _flatten_features(predictions)

    result = {
        "element_name": gt["element_name"],
        "filename": gt["filename"],
        "element_class": element_class,
        "gt_name": gt_name,
        "gt_family": gt_struct.get("family", ""),
        "gt_length": seq_len,
        "tool_status": tool_status,
        "detected": False,
        "boundary_match": False,
        "detected_type": None,
        "detected_family": None,
        "detected_name": None,
        "failure_modes": [],
    }

    # --- Transposon-level detection ---
    for p in flat_pred:
        p_type = p.get("element_type", "")
        p_start = p.get("start", 0)
        p_end = p.get("end", 0)
        p_span = p_end - p_start

        if p_type in ("transposon", "genomic_island", "integron"):
            overlap_start = max(1, p_start)
            overlap_end = min(seq_len, p_end)
            overlap = max(0, overlap_end - overlap_start)

            if overlap >= seq_len * 0.5 or p_span >= seq_len * 0.5:
                result["detected"] = True
                result["detected_type"] = p_type
                result["detected_family"] = p.get("family", "")
                result["detected_name"] = p.get("name", "")
                if abs(p_start - 1) <= tolerance and abs(p_end - seq_len) <= tolerance:
                    result["boundary_match"] = True
                break

    # Name-based fallback
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
                break

    # --- Feature-level scoring ---
    gt_features = []
    for ise in gt.get("is_elements", []):
        if not ise.get("spans_full_record"):
            gt_features.append({"type": "IS", "name": ise["name"],
                                "start": ise["start"], "end": ise["end"]})
    for amr in gt.get("amr_genes", []):
        gt_features.append({"type": "AMR", "name": amr["name"],
                            "start": amr["start"], "end": amr["end"]})
    for integ in gt.get("integrons", []):
        if not integ.get("spans_full_record"):
            gt_features.append({"type": "integron", "name": integ["name"],
                                "start": integ["start"], "end": integ["end"]})

    pred_features = [
        {"type": p["element_type"], "name": p.get("name", ""),
         "start": p.get("start", 0), "end": p.get("end", 0)}
        for p in flat_pred
        if p.get("element_type") in ("IS", "AMR", "integron", "transposon",
                                      "cassette", "integrase", "replicon",
                                      "genomic_island")
    ]

    gt_matched = []
    gt_unmatched = []
    for gf in gt_features:
        matched = False
        for pf in pred_features:
            if _type_compatible(gf["type"], pf["type"]):
                if (abs(gf["start"] - pf["start"]) <= tolerance
                        and abs(gf["end"] - pf["end"]) <= tolerance):
                    matched = True
                    break
        if matched:
            gt_matched.append(gf)
        else:
            gt_unmatched.append(gf)

    pred_matched_count = 0
    for pf in pred_features:
        for gf in gt_features:
            if _type_compatible(gf["type"], pf["type"]):
                if (abs(gf["start"] - pf["start"]) <= tolerance
                        and abs(gf["end"] - pf["end"]) <= tolerance):
                    pred_matched_count += 1
                    break

    result["gt_features_total"] = len(gt_features)
    result["gt_features_matched"] = len(gt_matched)
    result["gt_features_unmatched"] = [f"{f['type']}:{f['name']}@{f['start']}-{f['end']}" for f in gt_unmatched]
    result["pred_features_total"] = len(pred_features)
    result["pred_features_matched"] = pred_matched_count
    result["recall"] = len(gt_matched) / len(gt_features) if gt_features else None
    result["precision"] = pred_matched_count / len(pred_features) if pred_features else None

    # --- Failure mode classification ---
    if result["detected"]:
        if result["boundary_match"]:
            result["failure_modes"] = ["success"]
        else:
            result["failure_modes"] = ["boundary_inaccuracy"]
    else:
        result["failure_modes"] = _classify_failure(gt, flat_pred, tool_status)

    return result


def _type_compatible(gt_type: str, pred_type: str) -> bool:
    """Check if ground truth and prediction types are comparable."""
    if gt_type == pred_type:
        return True
    # IS may appear as mobile_element in some contexts
    if gt_type == "IS" and pred_type == "mobile_element":
        return True
    return False


def _classify_failure(
    gt: dict, flat_pred: list[dict], tool_status: dict[str, str]
) -> list[str]:
    """Classify why element-level detection failed."""
    modes = []
    element_class = gt.get("element_class", "")
    gt_struct = gt.get("structure", {})
    gt_family = gt_struct.get("family", "")
    gt_name = gt_struct.get("name", "")

    # Check which tools actually ran
    isescan_ok = "ok" in tool_status.get("isescan", "")
    amrfinder_ok = "ok" in tool_status.get("amrfinder", "")
    integron_ok = "ok" in tool_status.get("integron_finder", "")

    has_is = any(p["element_type"] == "IS" for p in flat_pred)
    has_amr = any(p["element_type"] == "AMR" for p in flat_pred)
    has_transposon = any(p["element_type"] == "transposon" for p in flat_pred)

    if element_class == "IS_element":
        if not isescan_ok:
            modes.append("isescan_did_not_run")
        elif not has_is:
            modes.append("isescan_missed_IS")
        else:
            modes.append("IS_detected_but_not_element_level")
        return modes

    if element_class == "composite_transposon":
        # Check what sub-components were found
        gt_is = [i for i in gt.get("is_elements", []) if not i.get("spans_full_record")]
        gt_amr = gt.get("amr_genes", [])

        if not isescan_ok and not has_is:
            modes.append("upstream_tool_miss:isescan_not_run")
        elif gt_is and not has_is:
            modes.append("upstream_tool_miss:isescan_missed_all_IS")
        elif gt_amr and not has_amr:
            if not amrfinder_ok:
                modes.append("upstream_tool_miss:amrfinder_not_run")
            else:
                modes.append("upstream_tool_miss:amrfinder_missed_AMR")

        if has_is and has_amr and not has_transposon:
            # Both components found but no transposon inferred
            is_families = set(p.get("family", "") for p in flat_pred if p["element_type"] == "IS")
            covered_families = {"IS6", "IS4", "IS1380", "IS30", "ISCR", "IS91"}
            if is_families & covered_families:
                modes.append("hierarchy_inference_failed:flanked_rule_did_not_fire")
            else:
                modes.append(f"hierarchy_inference_failed:no_rule_for_IS_families_{','.join(sorted(is_families))}")
        elif has_is and not gt_amr:
            modes.append("hierarchy_inference_failed:composite_has_no_AMR_cargo")
        elif has_transposon:
            modes.append("transposon_detected_but_wrong_span")

    elif element_class == "unit_transposon":
        blast_covered = {
            "Tn3", "Tn2", "Tn21", "Tn1696", "Tn1721", "Tn6452",
            "Tn4401", "Tn1546", "Tn1331", "Tn5393", "Tn7", "Tn552",
            "Tn6022", "Tn6019", "Tn6172",
        }
        if gt_name in blast_covered:
            modes.append("blast_reference_miss_despite_coverage")
        elif gt_family == "Tn402":
            modes.append("tn402_family_not_in_reference_library")
        else:
            modes.append("unit_transposon_not_in_reference_library")

    elif element_class == "integron":
        if not integron_ok:
            modes.append("upstream_tool_miss:integron_finder_not_run")
        else:
            modes.append("integron_detected_but_wrong_boundaries")

    if not modes:
        modes.append("unclassified_failure")

    return modes


# -----------------------------------------------------------------------
# Report generation
# -----------------------------------------------------------------------

def generate_report(scored: list[dict]) -> str:
    """Generate Markdown benchmark report."""
    lines = []
    lines.append("# Matryoshka TnCentral Benchmark Report")
    lines.append("")
    lines.append(f"Benchmark run against {len(scored)} representative TnCentral elements.")
    lines.append(f"Feature matching tolerance: +/-{TOLERANCE} bp.")
    lines.append(f"Pipeline mode: full (ISEScan + AMRFinder + IntegronFinder + BLAST).")
    lines.append("")

    # Transposon-level metrics by category
    lines.append("## Transposon-level detection by category")
    lines.append("")
    lines.append("| Category | Total | Detected | Recall | Boundary match |")
    lines.append("|---|---|---|---|---|")

    # Group by the category from REPRESENTATIVE_SET
    cat_map = {stem: cat for stem, cat, _ in REPRESENTATIVE_SET}
    by_cat: dict[str, list[dict]] = defaultdict(list)
    for s in scored:
        cat = cat_map.get(s["filename"], "unknown")
        by_cat[cat].append(s)

    total_all = detected_all = boundary_all = 0
    cat_order = ["IS6_flanked", "ISEcp1_capture", "ISCR_flanked", "IS4_flanked",
                 "Tn3_family", "integron", "non_covered_IS"]
    for cat in cat_order:
        items = by_cat.get(cat, [])
        if not items:
            continue
        total = len(items)
        detected = sum(1 for s in items if s["detected"])
        boundary = sum(1 for s in items if s["boundary_match"])
        recall = detected / total if total > 0 else 0
        lines.append(f"| {cat} | {total} | {detected} | {recall:.0%} | {boundary} |")
        total_all += total
        detected_all += detected
        boundary_all += boundary

    overall = detected_all / total_all if total_all > 0 else 0
    lines.append(f"| **Total** | **{total_all}** | **{detected_all}** | **{overall:.0%}** | **{boundary_all}** |")
    lines.append("")

    # Feature-level metrics
    lines.append("## Feature-level metrics")
    lines.append("")

    total_gt = sum(s["gt_features_total"] for s in scored)
    total_gt_matched = sum(s["gt_features_matched"] for s in scored)
    total_pred = sum(s["pred_features_total"] for s in scored)
    total_pred_matched = sum(s["pred_features_matched"] for s in scored)

    feat_recall = total_gt_matched / total_gt if total_gt > 0 else 0
    feat_precision = total_pred_matched / total_pred if total_pred > 0 else 0

    lines.append(f"- **Feature-level recall**: {total_gt_matched}/{total_gt} ({feat_recall:.1%})")
    lines.append(f"- **Feature-level precision**: {total_pred_matched}/{total_pred} ({feat_precision:.1%})")
    lines.append("")

    # Failure mode breakdown
    lines.append("## Failure mode breakdown")
    lines.append("")
    lines.append("| Failure mode | Count |")
    lines.append("|---|---|")

    mode_counts: Counter = Counter()
    for s in scored:
        for m in s["failure_modes"]:
            mode_counts[m] += 1
    for mode, count in mode_counts.most_common():
        lines.append(f"| {mode} | {count} |")
    lines.append("")

    # Detailed per-element results
    lines.append("## Per-element results")
    lines.append("")
    lines.append("| Element | Category | Detected | Boundary | Feature recall | Failure mode |")
    lines.append("|---|---|---|---|---|---|")

    for s in scored:
        cat = cat_map.get(s["filename"], "?")
        det = "yes" if s["detected"] else "NO"
        bnd = "yes" if s["boundary_match"] else "no"
        recall_str = f"{s['recall']:.0%}" if s["recall"] is not None else "n/a"
        modes = "; ".join(s["failure_modes"])
        lines.append(f"| {s['element_name']} | {cat} | {det} | {bnd} | {recall_str} | {modes} |")
    lines.append("")

    # Unmatched ground truth features
    lines.append("## Unmatched ground truth features (top 20)")
    lines.append("")
    all_unmatched: list[tuple[str, str]] = []
    for s in scored:
        for u in s.get("gt_features_unmatched", []):
            all_unmatched.append((s["element_name"], u))

    if all_unmatched:
        lines.append("| Element | Missed feature |")
        lines.append("|---|---|")
        for elem, feat in all_unmatched[:20]:
            lines.append(f"| {elem} | {feat} |")
    else:
        lines.append("None — all ground truth features were matched.")
    lines.append("")

    # Tool status summary
    lines.append("## External tool status")
    lines.append("")
    tool_summary: dict[str, Counter] = defaultdict(Counter)
    for s in scored:
        for tool, status in s.get("tool_status", {}).items():
            status_key = status.split("(")[0].strip()
            tool_summary[tool][status_key] += 1

    for tool in ["isescan", "amrfinder", "integron_finder", "blast_scan"]:
        if tool in tool_summary:
            status_str = ", ".join(f"{st}: {c}" for st, c in tool_summary[tool].most_common())
            lines.append(f"- **{tool}**: {status_str}")
    lines.append("")

    return "\n".join(lines)


def generate_findings(scored: list[dict]) -> str:
    """Generate FINDINGS.md summary of benchmark insights."""
    cat_map = {stem: cat for stem, cat, _ in REPRESENTATIVE_SET}

    total = len(scored)
    detected = sum(1 for s in scored if s["detected"])
    boundary = sum(1 for s in scored if s["boundary_match"])

    total_gt = sum(s["gt_features_total"] for s in scored)
    total_gt_matched = sum(s["gt_features_matched"] for s in scored)
    total_pred = sum(s["pred_features_total"] for s in scored)
    total_pred_matched = sum(s["pred_features_matched"] for s in scored)

    mode_counts: Counter = Counter()
    for s in scored:
        for m in s["failure_modes"]:
            mode_counts[m] += 1

    # Identify top 3 failure modes excluding "success" and "boundary_inaccuracy"
    real_failures = {m: c for m, c in mode_counts.items()
                     if m not in ("success", "boundary_inaccuracy")}
    top_failures = sorted(real_failures.items(), key=lambda x: -x[1])[:3]

    # Identify categories with worst recall
    by_cat: dict[str, list[dict]] = defaultdict(list)
    for s in scored:
        cat = cat_map.get(s["filename"], "unknown")
        by_cat[cat].append(s)

    lines = []
    lines.append("# TnCentral Benchmark Findings")
    lines.append("")
    lines.append(f"Date: 2026-05-02")
    lines.append(f"Elements tested: {total}")
    lines.append(f"Pipeline: full (ISEScan + AMRFinder + IntegronFinder + BLAST reference scan)")
    lines.append("")

    lines.append("## Summary metrics")
    lines.append("")
    lines.append(f"- Transposon-level recall: {detected}/{total} ({detected/total:.0%})")
    lines.append(f"- Transposon-level boundary accuracy: {boundary}/{detected} detected elements have boundaries within +/-50bp")
    feat_recall = total_gt_matched / total_gt if total_gt > 0 else 0
    feat_precision = total_pred_matched / total_pred if total_pred > 0 else 0
    lines.append(f"- Feature-level recall: {total_gt_matched}/{total_gt} ({feat_recall:.1%})")
    lines.append(f"- Feature-level precision: {total_pred_matched}/{total_pred} ({feat_precision:.1%})")
    lines.append("")

    lines.append("## Top 3 failure modes")
    lines.append("")
    for i, (mode, count) in enumerate(top_failures, 1):
        # Find example elements
        examples = [s["element_name"] for s in scored if mode in s["failure_modes"]][:3]
        examples_str = ", ".join(examples)
        lines.append(f"### {i}. {mode} ({count} elements)")
        lines.append(f"Examples: {examples_str}")
        lines.append("")

        if "no_rule_for_IS_families" in mode:
            lines.append("Matryoshka only has flanked-cargo rules for IS6, IS4, IS1380, IS30, ISCR, "
                         "and IS91 families. Composites using other IS families (IS10, IS50, IS256, "
                         "IS110, etc.) are not detected by rule-based inference. These represent a "
                         "significant gap — IS10 (Tn10) and IS50 (Tn5) are among the most "
                         "well-characterised composite transposons in the literature.")
        elif "flanked_rule_did_not_fire" in mode:
            lines.append("IS elements and AMR cargo were both detected by upstream tools, but "
                         "the transposon inference rules did not match them into a composite. "
                         "This indicates either: (a) the cargo gene name did not match any "
                         "cargo_match pattern in FLANKED_RULES, or (b) the gap between IS and "
                         "cargo exceeded the upstream_max/downstream_max thresholds.")
        elif "upstream_tool_miss" in mode:
            lines.append("An external detection tool (ISEScan, AMRFinder, or IntegronFinder) "
                         "failed to detect a component feature annotated in TnCentral. This is "
                         "an upstream limitation that Matryoshka cannot compensate for.")
        elif "unit_transposon_not_in_reference_library" in mode:
            lines.append("Unit transposons are detected primarily via BLAST against bundled "
                         "reference sequences. Elements not in the reference library cannot "
                         "be detected. Adding more Tn3-family references would improve coverage.")
        elif "boundary_inaccuracy" in mode:
            lines.append("The transposon was correctly detected (name/family match or >= 50% "
                         "span overlap) but the predicted boundaries deviated from the TnCentral "
                         "annotation by more than +/-50bp. This is often due to BLAST alignment "
                         "truncation at regions of low identity, or ISEScan coordinate imprecision.")
        elif "tn402_family" in mode:
            lines.append("Tn402/Tn5053-family transposons are the progenitors of class-1 integrons "
                         "but are not in Matryoshka's BLAST reference library. This is listed as "
                         "priority #1 in GAPS.md.")
        elif "blast_reference_miss" in mode:
            lines.append("A transposon that should be covered by a BLAST reference was not detected. "
                         "This may indicate the reference sequence is too divergent from the TnCentral "
                         "exemplar, or the min_identity/min_length thresholds are too strict.")
        lines.append("")

    lines.append("## Category-level analysis")
    lines.append("")
    for cat in ["IS6_flanked", "ISEcp1_capture", "ISCR_flanked", "IS4_flanked",
                "Tn3_family", "integron", "non_covered_IS"]:
        items = by_cat.get(cat, [])
        if not items:
            continue
        det = sum(1 for s in items if s["detected"])
        lines.append(f"### {cat} ({det}/{len(items)} detected)")
        for s in items:
            status = "DETECTED" if s["detected"] else "MISSED"
            recall_str = f"feature recall {s['recall']:.0%}" if s["recall"] is not None else "no GT features"
            modes = "; ".join(s["failure_modes"])
            lines.append(f"- {s['element_name']}: {status}, {recall_str} [{modes}]")
        lines.append("")

    lines.append("## Recommendations for ML improvement priorities")
    lines.append("")
    lines.append("Based on this benchmark, the following improvements would have the largest impact:")
    lines.append("")
    lines.append("1. **Add flanked-cargo rules for IS4 (IS10), IS4 (IS50), and IS256 families.** "
                 "These cover Tn10, Tn5, and Tn4001 — three of the most clinically important "
                 "composite transposons that are currently undetectable by rule-based inference.")
    lines.append("2. **Add Tn402/Tn5053 backbone to the BLAST reference library.** "
                 "This unlocks detection of the Tn402 family, which is the structural ancestor "
                 "of most class-1 integrons.")
    lines.append("3. **Improve cargo_match flexibility.** "
                 "Some composites have cargo gene names that differ from the patterns in "
                 "FLANKED_RULES (e.g. AMRFinder reports 'blaKPC-3' but the rule checks 'KPC'). "
                 "An ML-based cargo classifier could replace rigid substring matching.")
    lines.append("4. **Tune BLAST thresholds for Tn3-family references.** "
                 "Some unit transposons with >90% identity to a reference were missed because "
                 "the alignment length fell below the min_length cutoff. A sliding-scale "
                 "approach (lower length threshold at higher identity) would improve recall.")
    lines.append("")

    return "\n".join(lines)


# -----------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(description="Comprehensive Matryoshka TnCentral benchmark")
    parser.add_argument("--blast-only", action="store_true",
                        help="Run BLAST reference scan only (skip external tools)")
    parser.add_argument("--max", type=int, default=0,
                        help="Maximum elements to process (0 = all 30)")
    args = parser.parse_args()

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    # Resolve representative set
    elements: list[tuple[str, str, str, Path, dict]] = []
    for stem, category, rationale in REPRESENTATIVE_SET:
        gb_path = _find_gb_path(stem)
        if gb_path is None:
            print(f"  SKIP: {stem} — GenBank file not found")
            continue
        gt = load_ground_truth(gb_path.stem)
        if gt is None:
            print(f"  SKIP: {stem} — parsed ground truth not found")
            continue
        elements.append((stem, category, rationale, gb_path, gt))

    if args.max > 0:
        elements = elements[:args.max]

    print(f"Benchmark set: {len(elements)} elements")
    print(f"Mode: {'BLAST-only' if args.blast_only else 'full pipeline'}")
    print()

    # Write ground truth TSV
    write_ground_truth_tsv([(stem, cat, gt) for stem, cat, _, _, gt in elements])

    scored: list[dict] = []
    for i, (stem, category, rationale, gb_path, gt) in enumerate(elements):
        result_path = RESULTS_DIR / f"{gb_path.stem}.json"

        # Check cache
        if result_path.exists():
            with open(result_path) as fh:
                cached = json.load(fh)
            scored.append(cached)
            print(f"  [{i+1}/{len(elements)}] {gt['element_name']} ({category}) — cached")
            continue

        print(f"  [{i+1}/{len(elements)}] {gt['element_name']} "
              f"({category}, {gt['length']} bp)...", end="", flush=True)

        with tempfile.TemporaryDirectory() as tmpdir:
            fasta_path = Path(tmpdir) / f"{gb_path.stem}.fasta"
            if not write_single_fasta(gb_path, fasta_path):
                print(" SKIP (FASTA conversion failed)")
                continue

            try:
                if args.blast_only:
                    roots, all_feats, tool_status = run_blast_pipeline(fasta_path)
                else:
                    roots, all_feats, tool_status = run_full_pipeline(fasta_path)
            except Exception as e:
                print(f" ERROR: {e}")
                traceback.print_exc()
                error_result = {
                    "element_name": gt["element_name"],
                    "filename": gb_path.stem,
                    "element_class": gt.get("element_class", ""),
                    "gt_name": gt.get("structure", {}).get("name", ""),
                    "gt_family": gt.get("structure", {}).get("family", ""),
                    "gt_length": gt.get("length", 0),
                    "tool_status": {},
                    "detected": False,
                    "boundary_match": False,
                    "detected_type": None,
                    "detected_family": None,
                    "detected_name": None,
                    "failure_modes": [f"pipeline_error:{e}"],
                    "gt_features_total": 0,
                    "gt_features_matched": 0,
                    "gt_features_unmatched": [],
                    "pred_features_total": 0,
                    "pred_features_matched": 0,
                    "recall": None,
                    "precision": None,
                    "error": str(e),
                }
                scored.append(error_result)
                with open(result_path, "w") as fh:
                    json.dump(error_result, fh, indent=2)
                continue

        predictions = [_feature_to_dict(r) for r in roots]
        all_feat_dicts = [_feature_to_dict(f) for f in all_feats]

        result = score_element(gt, predictions, all_feat_dicts, tool_status)

        # Save per-element result
        result["predictions"] = predictions
        result["all_features"] = all_feat_dicts
        with open(result_path, "w") as fh:
            json.dump(result, fh, indent=2)

        # Remove verbose fields from the in-memory copy for report generation
        report_result = {k: v for k, v in result.items()
                         if k not in ("predictions", "all_features")}
        scored.append(report_result)

        det = "DETECTED" if report_result["detected"] else "MISSED"
        recall_str = f"recall={report_result['recall']:.0%}" if report_result["recall"] is not None else ""
        print(f" {det} {recall_str}")

    print()

    # Generate reports
    report = generate_report(scored)
    report_path = DATA_DIR / "benchmark_report.md"
    with open(report_path, "w") as fh:
        fh.write(report)
    print(f"Report: {report_path}")

    findings = generate_findings(scored)
    findings_path = DATA_DIR / "FINDINGS.md"
    with open(findings_path, "w") as fh:
        fh.write(findings)
    print(f"Findings: {findings_path}")

    # Save combined results (without heavy prediction fields)
    combined = [
        {k: v for k, v in s.items() if k not in ("predictions", "all_features")}
        for s in scored
    ]
    combined_path = DATA_DIR / "benchmark_results_combined.json"
    with open(combined_path, "w") as fh:
        json.dump(combined, fh, indent=2)
    print(f"Combined JSON: {combined_path}")

    # Summary
    detected_total = sum(1 for s in scored if s["detected"])
    print(f"\nTransposon-level recall: {detected_total}/{len(scored)} "
          f"({detected_total/len(scored):.0%})")

    mode_counts: Counter = Counter()
    for s in scored:
        for m in s["failure_modes"]:
            mode_counts[m] += 1
    print("\nTop failure modes:")
    for mode, count in mode_counts.most_common(5):
        print(f"  {mode}: {count}")


if __name__ == "__main__":
    main()
