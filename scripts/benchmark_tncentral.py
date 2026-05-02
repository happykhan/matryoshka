#!/usr/bin/env python3
"""
benchmark_tncentral.py — Run Matryoshka on TnCentral sequences and compare
predictions against ground-truth annotations.

Modes:
  --mode=blast     Reference-scan only (fast, ~1 min for all 514 sequences)
  --mode=full      Full pipeline: ISEScan + AMRFinder + IntegronFinder + BLAST
                   (slow, ~2-5 min per sequence x 514 = hours)
  --mode=subset    Full pipeline on a filtered subset (default: composites + integrons
                   with AMR genes, ~80 records)

Checkpointing: results are saved per-sequence so interrupted runs can resume.

Usage:
    pixi run python3 scripts/benchmark_tncentral.py --mode=blast
    pixi run python3 scripts/benchmark_tncentral.py --mode=subset --max=50
    pixi run python3 scripts/benchmark_tncentral.py --mode=full --max=20
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
import tempfile
import traceback
from pathlib import Path

from Bio import SeqIO

# Add project root to path
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
FASTA_DIR = DATA_DIR / "fasta" / "fa"
GB_DIR = DATA_DIR / "genbank" / "gb"


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


def run_blast_only(fasta_path: Path) -> list[MGEFeature]:
    """Run only the BLAST reference scan against bundled Matryoshka references."""
    if not blast_available():
        return []
    return scan_all(str(fasta_path))


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
        return None

    return None


def run_full_pipeline(fasta_path: Path) -> tuple[list[MGEFeature], list[MGEFeature]]:
    """Run the full Matryoshka pipeline on a single FASTA sequence.

    Returns (roots, all_features).
    """
    all_features: list[MGEFeature] = []

    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = Path(tmpdir)

        # ISEScan
        isescan_out = run_external_tool("isescan", fasta_path, tmp / "isescan")
        if isescan_out:
            all_features.extend(parse_isescan(isescan_out))

        # AMRFinder
        amrfinder_out = run_external_tool("amrfinder", fasta_path, tmp / "amrfinder")
        if amrfinder_out:
            all_features.extend(parse_amrfinder(amrfinder_out))

        # IntegronFinder
        integron_out = run_external_tool("integron_finder", fasta_path, tmp / "integron")
        if integron_out:
            try:
                all_features.extend(parse_integron_finder(integron_out))
            except Exception:
                pass

    # BLAST reference scan
    if blast_available():
        ref_hits = scan_all(str(fasta_path))
        all_features.extend(ref_hits)

    # Run inference
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

    return roots, all_feats


def run_blast_pipeline(fasta_path: Path) -> tuple[list[MGEFeature], list[MGEFeature]]:
    """Run Matryoshka with BLAST reference scan only (no external tools).

    Returns (roots, all_features).
    """
    all_features: list[MGEFeature] = run_blast_only(fasta_path)

    seq = str(next(SeqIO.parse(fasta_path, "fasta")).seq)

    # No IS/AMR features means no rule-based inference,
    # but we still run the pipeline stages for consistency
    inferred = infer_transposons(all_features)
    all_feats = all_features + inferred
    all_feats = _suppress_redundant_inference(all_feats)

    transposons = [f for f in all_feats if f.element_type == "transposon"]
    all_feats.extend(annotate_res_sites(transposons))

    checkable = [f for f in all_feats if f.element_type in ("IS", "transposon")]
    confirm_boundaries(seq, checkable)
    assign_confidence(all_feats)
    roots = build_hierarchy(all_feats)

    return roots, all_feats


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


def process_one(
    element_name: str, gb_path: Path, mode: str
) -> dict | None:
    """Process a single TnCentral record. Returns results dict or None."""
    result_path = RESULTS_DIR / f"{gb_path.stem}.json"
    if result_path.exists():
        with open(result_path) as fh:
            return json.load(fh)

    with tempfile.TemporaryDirectory() as tmpdir:
        fasta_path = Path(tmpdir) / f"{gb_path.stem}.fasta"
        if not write_single_fasta(gb_path, fasta_path):
            return None

        try:
            if mode == "blast":
                roots, all_feats = run_blast_pipeline(fasta_path)
            else:
                roots, all_feats = run_full_pipeline(fasta_path)
        except Exception as e:
            result = {
                "element_name": element_name,
                "filename": gb_path.stem,
                "error": str(e),
                "traceback": traceback.format_exc(),
                "predictions": [],
                "all_features": [],
            }
            with open(result_path, "w") as fh:
                json.dump(result, fh, indent=2)
            return result

    result = {
        "element_name": element_name,
        "filename": gb_path.stem,
        "error": None,
        "predictions": [_feature_to_dict(r) for r in roots],
        "all_features": [_feature_to_dict(f) for f in all_feats],
    }
    with open(result_path, "w") as fh:
        json.dump(result, fh, indent=2)
    return result


def select_subset(parsed_records: list[dict]) -> list[dict]:
    """Select a useful benchmark subset: composites and integrons with AMR."""
    subset = []
    for rec in parsed_records:
        cls = rec.get("element_class", "")
        if cls == "composite_transposon":
            subset.append(rec)
        elif cls == "integron" and rec.get("amr_genes"):
            subset.append(rec)
        elif cls == "unit_transposon" and rec.get("amr_genes"):
            # Include Tn3-family units with AMR — these test reference scan
            subset.append(rec)
    return subset


def main() -> None:
    parser = argparse.ArgumentParser(description="Benchmark Matryoshka against TnCentral")
    parser.add_argument("--mode", choices=["blast", "full", "subset"], default="blast",
                        help="Pipeline mode (default: blast)")
    parser.add_argument("--max", type=int, default=0,
                        help="Maximum sequences to process (0 = all)")
    parser.add_argument("--clean", action="store_true",
                        help="Remove cached results and start fresh")
    args = parser.parse_args()

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    if args.clean:
        for f in RESULTS_DIR.glob("*.json"):
            f.unlink()

    # Load parsed records
    parsed_files = sorted(PARSED_DIR.glob("*.json"))
    if not parsed_files:
        print("No parsed TnCentral records found. Run parse_tncentral.py first.")
        return

    parsed_records = []
    for pf in parsed_files:
        with open(pf) as fh:
            parsed_records.append(json.load(fh))

    # Filter to subset if requested
    if args.mode == "subset":
        records_to_process = select_subset(parsed_records)
        print(f"Subset mode: {len(records_to_process)} records selected "
              f"(composites + integrons + unit Tn with AMR)")
    else:
        records_to_process = parsed_records

    if args.max > 0:
        records_to_process = records_to_process[:args.max]

    pipeline_mode = "blast" if args.mode == "blast" else "full"
    print(f"Mode: {pipeline_mode}")
    print(f"Processing {len(records_to_process)} sequences...")

    results: list[dict] = []
    errors = 0
    cached = 0

    for i, rec in enumerate(records_to_process):
        filename = rec["filename"]
        gb_path = GB_DIR / f"{filename}.gb"
        if not gb_path.exists():
            continue

        # Check cache
        result_path = RESULTS_DIR / f"{filename}.json"
        if result_path.exists():
            with open(result_path) as fh:
                result = json.load(fh)
            results.append(result)
            cached += 1
            continue

        print(f"  [{i+1}/{len(records_to_process)}] {rec['element_name']} "
              f"({rec['length']} bp, {rec.get('element_class', '?')})...",
              end="", flush=True)

        result = process_one(rec["element_name"], gb_path, pipeline_mode)
        if result is None:
            print(" SKIP")
            continue
        if result.get("error"):
            print(f" ERROR: {result['error'][:60]}")
            errors += 1
        else:
            n_pred = len(result.get("predictions", []))
            n_all = len(result.get("all_features", []))
            print(f" OK ({n_pred} roots, {n_all} features)")
        results.append(result)

    print(f"\nProcessed: {len(results)} (cached: {cached}, errors: {errors})")

    # Save combined results
    combined_path = DATA_DIR / "benchmark_results_combined.json"
    with open(combined_path, "w") as fh:
        json.dump(results, fh, indent=2)
    print(f"Combined results: {combined_path}")


if __name__ == "__main__":
    main()
