#!/usr/bin/env python3
"""
parse_tncentral.py — Extract ground-truth annotations from TnCentral GenBank files.

Parses TnCentral's annotated transposon GenBank records to produce:
  1. A structured JSON file per transposon (in data/tncentral/parsed/)
  2. A summary TSV listing all transposons and their feature counts
  3. A combined FASTA of all sequences (for batch Matryoshka input)

TnCentral GenBank annotations use:
  - mobile_element features with mobile_element_type qualifiers
    (e.g. "transposon:Tn21", "insertion sequence:IS26")
  - CDS features with ARO: tags in gene qualifiers for AMR genes
  - repeat_region for IRs and TSDs (labelled via note/label qualifiers)
  - misc_feature for attC sites, res sites, etc.
"""

from __future__ import annotations

import json
import re
import warnings
from pathlib import Path

from Bio import SeqIO

warnings.filterwarnings("ignore", category=UserWarning)

DATA_DIR = Path(__file__).parent.parent / "data" / "tncentral"
GB_DIR = DATA_DIR / "genbank" / "gb"
PARSED_DIR = DATA_DIR / "parsed"
FASTA_OUT = DATA_DIR / "tncentral_sequences.fasta"
SUMMARY_TSV = DATA_DIR / "tncentral_summary.tsv"


def _parse_mobile_element_type(qual: str) -> tuple[str, str]:
    """Parse mobile_element_type qualifier into (type, name).

    Examples:
        "transposon:Tn21" -> ("transposon", "Tn21")
        "insertion sequence:IS26" -> ("IS", "IS26")
        "integron:In2" -> ("integron", "In2")
    """
    if ":" in qual:
        raw_type, name = qual.split(":", 1)
    else:
        raw_type, name = qual, qual
    raw_type = raw_type.strip().lower()
    if raw_type == "insertion sequence":
        raw_type = "IS"
    elif raw_type == "other":
        raw_type = "mobile_element"
    return raw_type, name.strip()


def _parse_note_metadata(note: str) -> dict[str, str]:
    """Parse TnCentral-style note field into key-value pairs.

    Notes follow the pattern: "Key1 = Value1; Key2 = Value2; ..."
    """
    meta: dict[str, str] = {}
    for pair in note.split(";"):
        pair = pair.strip()
        if "=" in pair:
            k, v = pair.split("=", 1)
            meta[k.strip()] = v.strip()
    return meta


def _is_amr_gene(feat) -> bool:
    """Check whether a CDS feature represents an AMR gene via ARO tag."""
    gene = feat.qualifiers.get("gene", [""])[0]
    return "ARO:" in gene


def _extract_amr_name(gene_str: str) -> str:
    """Extract clean AMR gene name, stripping the ARO accession.

    Example: "sul1 (ARO:3000410)" -> "sul1"
    """
    return re.sub(r"\s*\(ARO:\d+\)", "", gene_str).strip()


def _extract_ir_info(feat) -> dict:
    """Extract IR information from a repeat_region feature."""
    label = feat.qualifiers.get("label", [""])[0]
    note = feat.qualifiers.get("note", [""])[0]
    ir_type = ""
    associated = ""
    if label.startswith("IRL"):
        ir_type = "IRL"
    elif label.startswith("IRR"):
        ir_type = "IRR"
    elif label.startswith("IRt"):
        ir_type = "IRt"
    elif label.startswith("IRi"):
        ir_type = "IRi"
    elif label.startswith("IR "):
        ir_type = "IR"
    elif "DR" in label:
        ir_type = "TSD"

    meta = _parse_note_metadata(note)
    associated = meta.get("AssociatedElement", "")

    return {
        "type": ir_type,
        "label": label,
        "associated_element": associated,
        "start": int(feat.location.start) + 1,  # convert to 1-based
        "end": int(feat.location.end),
        "strand": "+" if feat.location.strand == 1 else "-",
    }


def parse_genbank_record(gb_path: Path) -> dict | None:
    """Parse a single TnCentral GenBank record into a structured dict."""
    try:
        rec = next(SeqIO.parse(gb_path, "genbank"))
    except Exception:
        return None

    filename = gb_path.stem
    # Extract element name and accession from filename
    # Format: ElementName-Accession.gb (e.g. "Tn21-AF071413.gb")
    parts = filename.split("-", 1)
    element_name = parts[0]
    accession = parts[1] if len(parts) > 1 else ""

    result = {
        "filename": filename,
        "element_name": element_name,
        "accession": accession,
        "description": rec.description,
        "length": len(rec.seq),
        "sequence_id": rec.id,
        "mobile_elements": [],
        "is_elements": [],
        "amr_genes": [],
        "integrons": [],
        "repeat_regions": [],
        "res_sites": [],
        "other_cds": [],
        "structure": {},  # top-level element metadata
    }

    for feat in rec.features:
        loc_start = int(feat.location.start) + 1  # 1-based
        loc_end = int(feat.location.end)
        strand = "+" if feat.location.strand == 1 else ("-" if feat.location.strand == -1 else ".")

        if feat.type == "mobile_element":
            me_type_raw = feat.qualifiers.get("mobile_element_type", [""])[0]
            label = feat.qualifiers.get("label", [""])[0]
            note = feat.qualifiers.get("note", [""])[0]
            me_type, me_name = _parse_mobile_element_type(me_type_raw)
            meta = _parse_note_metadata(note)

            entry = {
                "type": me_type,
                "name": me_name,
                "label": label,
                "family": meta.get("Family", ""),
                "group": meta.get("Group", ""),
                "accession": meta.get("Accession", ""),
                "partial": meta.get("Partial", "").strip(),
                "start": loc_start,
                "end": loc_end,
                "strand": strand,
                "spans_full_record": (loc_start == 1 and loc_end == len(rec.seq)),
            }

            if me_type == "IS":
                result["is_elements"].append(entry)
            elif me_type == "integron":
                result["integrons"].append(entry)
            elif me_type == "transposon":
                result["mobile_elements"].append(entry)
            else:
                result["mobile_elements"].append(entry)

            # If this spans the full record, it's the top-level structure
            if entry["spans_full_record"]:
                result["structure"] = {
                    "type": me_type,
                    "name": me_name,
                    "family": entry["family"],
                    "group": entry["group"],
                }

        elif feat.type == "CDS":
            if _is_amr_gene(feat):
                gene = feat.qualifiers.get("gene", [""])[0]
                product = feat.qualifiers.get("product", [""])[0]
                note = feat.qualifiers.get("note", [""])[0]
                meta = _parse_note_metadata(note)
                result["amr_genes"].append({
                    "name": _extract_amr_name(gene),
                    "product": product,
                    "start": loc_start,
                    "end": loc_end,
                    "strand": strand,
                    "class": meta.get("Subclass", ""),
                    "associated_element": meta.get("AssociatedElement", ""),
                })
            else:
                gene = feat.qualifiers.get("gene", [""])[0]
                product = feat.qualifiers.get("product", [""])[0]
                note = feat.qualifiers.get("note", [""])[0]
                meta = _parse_note_metadata(note)
                func = feat.qualifiers.get("function", [""])[0]

                entry = {
                    "gene": gene,
                    "product": product,
                    "function": func,
                    "start": loc_start,
                    "end": loc_end,
                    "strand": strand,
                }

                # Classify transposase genes
                if "transposase" in func.lower() or "transposase" in product.lower():
                    entry["role"] = "transposase"
                elif "resolvase" in func.lower() or gene.lower().startswith("tnpr"):
                    entry["role"] = "resolvase"
                elif "integrase" in func.lower() or gene.lower().startswith("inti"):
                    entry["role"] = "integrase"
                else:
                    entry["role"] = "other"

                result["other_cds"].append(entry)

        elif feat.type == "repeat_region":
            ir_info = _extract_ir_info(feat)
            result["repeat_regions"].append(ir_info)

        elif feat.type == "misc_feature":
            label = feat.qualifiers.get("label", [""])[0]
            note = feat.qualifiers.get("note", [""])[0]
            if "res" in label.lower() and ("site" in label.lower() or "res_site" in note.lower()):
                result["res_sites"].append({
                    "label": label,
                    "start": loc_start,
                    "end": loc_end,
                    "strand": strand,
                })

    return result


def classify_element(parsed: dict) -> str:
    """Classify the element for benchmark grouping.

    Returns one of: composite_transposon, unit_transposon, IS_element,
    integron, compound, other
    """
    struct = parsed.get("structure", {})
    family = struct.get("family", "")
    etype = struct.get("type", "")

    # IS elements (standalone)
    if etype == "IS":
        return "IS_element"
    if (parsed["element_name"].startswith("IS")
            and not parsed["element_name"].startswith("ISEcp")
            and not any(m["type"] == "transposon" for m in parsed["mobile_elements"])):
        return "IS_element"

    # Integrons
    if etype == "integron" or parsed["element_name"].startswith("In"):
        return "integron"

    # Compound / composite transposons
    if family in ("Compound Transposon", "Compound", "compound"):
        return "composite_transposon"
    # Composite if named Tn and has >=2 IS elements inside
    internal_is = [
        i for i in parsed["is_elements"]
        if not i.get("spans_full_record", False)
    ]
    if etype == "transposon" and len(internal_is) >= 2:
        return "composite_transposon"

    # Unit transposons (Tn3-family etc)
    if etype == "transposon" or parsed["element_name"].startswith("Tn"):
        return "unit_transposon"

    return "other"


def main() -> None:
    PARSED_DIR.mkdir(parents=True, exist_ok=True)

    if not GB_DIR.exists():
        print(f"GenBank directory not found: {GB_DIR}")
        print("Run the download step first.")
        return

    gb_files = sorted(GB_DIR.glob("*.gb"))
    print(f"Found {len(gb_files)} GenBank files to parse")

    all_records: list[dict] = []
    fasta_records = []

    for gb_file in gb_files:
        parsed = parse_genbank_record(gb_file)
        if parsed is None:
            print(f"  SKIP: {gb_file.name} (parse error)")
            continue

        parsed["element_class"] = classify_element(parsed)
        all_records.append(parsed)

        # Write individual JSON
        json_path = PARSED_DIR / f"{parsed['filename']}.json"
        with open(json_path, "w") as fh:
            json.dump(parsed, fh, indent=2)

        # Collect FASTA
        try:
            rec = next(SeqIO.parse(gb_file, "genbank"))
            # Use a clean ID: element_name__accession
            clean_id = f"{parsed['element_name']}__{parsed['accession']}"
            rec.id = clean_id
            rec.description = f"{parsed['element_name']} {parsed['accession']} len={len(rec.seq)}"
            fasta_records.append(rec)
        except Exception:
            pass

    # Write combined FASTA
    with open(FASTA_OUT, "w") as fh:
        SeqIO.write(fasta_records, fh, "fasta")
    print(f"Wrote {len(fasta_records)} sequences to {FASTA_OUT}")

    # Write summary TSV
    with open(SUMMARY_TSV, "w") as fh:
        headers = [
            "element_name", "accession", "element_class", "length",
            "num_mobile_elements", "num_is_elements", "num_amr_genes",
            "num_integrons", "num_repeat_regions", "num_res_sites",
            "family", "group",
            "is_names", "amr_names",
        ]
        fh.write("\t".join(headers) + "\n")
        for rec in all_records:
            struct = rec.get("structure", {})
            is_names = ",".join(i["name"] for i in rec["is_elements"])
            amr_names = ",".join(a["name"] for a in rec["amr_genes"])
            row = [
                rec["element_name"],
                rec["accession"],
                rec["element_class"],
                str(rec["length"]),
                str(len(rec["mobile_elements"])),
                str(len(rec["is_elements"])),
                str(len(rec["amr_genes"])),
                str(len(rec["integrons"])),
                str(len(rec["repeat_regions"])),
                str(len(rec["res_sites"])),
                struct.get("family", ""),
                struct.get("group", ""),
                is_names,
                amr_names,
            ]
            fh.write("\t".join(row) + "\n")
    print(f"Wrote summary to {SUMMARY_TSV}")

    # Print statistics
    from collections import Counter
    classes = Counter(r["element_class"] for r in all_records)
    print("\nElement class distribution:")
    for cls, count in classes.most_common():
        print(f"  {cls}: {count}")

    amr_count = sum(1 for r in all_records if r["amr_genes"])
    is_count = sum(1 for r in all_records if r["is_elements"])
    composite_with_is = sum(
        1 for r in all_records
        if r["element_class"] == "composite_transposon" and len(r["is_elements"]) >= 2
    )
    print(f"\nRecords with AMR genes: {amr_count}")
    print(f"Records with IS elements: {is_count}")
    print(f"Composite transposons with >=2 IS: {composite_with_is}")


if __name__ == "__main__":
    main()
