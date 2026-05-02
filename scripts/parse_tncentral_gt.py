#!/usr/bin/env python3
"""Parse TnCentral GenBank files into ground-truth annotations.

Outputs:
  data/tncentral/ground_truth.tsv   — one row per annotated feature
  data/tncentral/sequences/         — one FASTA per transposon
  data/tncentral/summary.tsv        — one row per transposon with counts
"""

import csv
import re
import sys
from pathlib import Path

from Bio import SeqIO

GB_DIR = Path(__file__).parent.parent / "data" / "tncentral" / "gb"
OUT_DIR = Path(__file__).parent.parent / "data" / "tncentral"
SEQ_DIR = OUT_DIR / "sequences"


def extract_dr_length(note: str) -> int | None:
    m = re.search(r'DR\s+(\d+)\s*bp', note, re.IGNORECASE)
    if m:
        return int(m.group(1))
    return None


def classify_feature(feat):
    """Return (element_type, name) from a GenBank feature."""
    ftype = feat.type
    qualifiers = feat.qualifiers

    if ftype == "mobile_element":
        me_type = qualifiers.get("mobile_element_type", [""])[0]
        label = qualifiers.get("label", [""])[0]
        if "transposon" in me_type.lower():
            name = me_type.replace("transposon:", "").strip() or label
            return "transposon", name
        elif "insertion sequence" in me_type.lower():
            name = me_type.replace("insertion sequence:", "").strip() or label
            return "IS", name
        else:
            return "mobile_element", label

    elif ftype == "CDS":
        gene = qualifiers.get("gene", qualifiers.get("label", [""]))[0]
        product = qualifiers.get("product", [""])[0].lower()
        # Classify by product
        if any(k in product for k in ("transposase", "integrase", "resolvase", "recombinase")):
            return "transposase", gene
        elif any(k in product for k in ("resistance", "beta-lactam", "aminoglycoside",
                                         "tetracycline", "chloramphenicol", "sulfonamide",
                                         "trimethoprim", "quinolone", "macrolide", "colistin")):
            return "AMR", gene
        else:
            return "CDS", gene

    elif ftype == "repeat_region":
        label = qualifiers.get("label", ["repeat_region"])[0]
        note = qualifiers.get("note", [""])[0]
        if "IR" in label or "inverted" in note.lower():
            return "IR", label
        elif "DR" in label or "direct" in note.lower():
            return "DR", label
        return "repeat", label

    elif ftype == "integron":
        return "integron", qualifiers.get("label", ["integron"])[0]

    return None, None


def parse_gb(gb_path: Path) -> tuple[dict, list[dict]]:
    """Parse a single GenBank file. Returns (meta, features)."""
    record = SeqIO.read(gb_path, "genbank")
    tn_name = gb_path.stem  # e.g. Tn1999.1-CP019080

    # Top-level mobile_element (the transposon itself)
    top_note = ""
    for feat in record.features:
        if feat.type == "mobile_element":
            me_type = feat.qualifiers.get("mobile_element_type", [""])[0]
            if "transposon" in me_type.lower():
                top_note = feat.qualifiers.get("note", [""])[0]
                break

    dr_length = extract_dr_length(top_note) if top_note else None
    family = ""
    group = ""
    if top_note:
        m = re.search(r'Family\s*=\s*([^;]+)', top_note)
        if m:
            family = m.group(1).strip()
        m = re.search(r'Group\s*=\s*([^;]+)', top_note)
        if m:
            group = m.group(1).strip()

    meta = {
        "tn_name": tn_name,
        "length": len(record.seq),
        "family": family,
        "group": group,
        "dr_length": dr_length,
        "accession": record.id,
    }

    features = []
    for feat in record.features:
        if feat.type == "source":
            continue
        etype, ename = classify_feature(feat)
        if etype is None:
            continue
        loc = feat.location
        start = int(loc.start) + 1  # 1-based
        end = int(loc.end)
        strand = "+" if loc.strand >= 0 else "-"
        features.append({
            "tn_name": tn_name,
            "element_type": etype,
            "element_name": ename,
            "start": start,
            "end": end,
            "strand": strand,
        })

    return meta, features


def main():
    SEQ_DIR.mkdir(parents=True, exist_ok=True)
    gb_files = sorted(GB_DIR.glob("*.gb"))
    print(f"Parsing {len(gb_files)} GenBank files...", file=sys.stderr)

    all_meta = []
    all_features = []
    errors = []

    for gb_path in gb_files:
        try:
            record = SeqIO.read(gb_path, "genbank")
            meta, features = parse_gb(gb_path)
            all_meta.append(meta)
            all_features.extend(features)

            # Write FASTA
            fa_path = SEQ_DIR / (gb_path.stem + ".fasta")
            with open(fa_path, "w") as fh:
                fh.write(f">{meta['tn_name']}\n{record.seq}\n")

        except Exception as e:
            errors.append((gb_path.name, str(e)))

    # Write ground truth TSV
    gt_path = OUT_DIR / "ground_truth.tsv"
    with open(gt_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=[
            "tn_name", "element_type", "element_name", "start", "end", "strand"
        ], delimiter="\t")
        writer.writeheader()
        writer.writerows(all_features)

    # Write summary TSV
    summary_path = OUT_DIR / "summary.tsv"
    with open(summary_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=[
            "tn_name", "length", "family", "group", "dr_length", "accession",
            "n_IS", "n_AMR", "n_integron", "n_IR", "n_transposase"
        ], delimiter="\t")
        writer.writeheader()
        for meta in all_meta:
            tn = meta["tn_name"]
            feats = [f for f in all_features if f["tn_name"] == tn]
            meta["n_IS"] = sum(1 for f in feats if f["element_type"] == "IS")
            meta["n_AMR"] = sum(1 for f in feats if f["element_type"] == "AMR")
            meta["n_integron"] = sum(1 for f in feats if f["element_type"] == "integron")
            meta["n_IR"] = sum(1 for f in feats if f["element_type"] == "IR")
            meta["n_transposase"] = sum(1 for f in feats if f["element_type"] == "transposase")
            writer.writerow(meta)

    print(f"Parsed {len(all_meta)} transposons, {len(all_features)} features", file=sys.stderr)
    print(f"Errors: {len(errors)}", file=sys.stderr)
    for name, err in errors[:10]:
        print(f"  {name}: {err}", file=sys.stderr)
    print(f"Ground truth: {gt_path}", file=sys.stderr)
    print(f"Summary: {summary_path}", file=sys.stderr)
    print(f"Sequences: {SEQ_DIR}/", file=sys.stderr)


if __name__ == "__main__":
    main()
