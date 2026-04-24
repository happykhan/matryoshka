"""
output.py — Write Matryoshka results in multiple formats.

Formats:
  - Wolvercote (cell format) — compact nested notation
  - SVG — schematic circular diagram via wolvercote renderer
  - GFF3 with IR / TSD / attribute annotations
  - JSON hierarchy (full feature payload)
  - GenBank (one record per contig with nested features as sub-features)
"""

from __future__ import annotations

import json
from io import StringIO
from typing import Any

from .detect import MGEFeature

# ---------------------------------------------------------------------------
# Wolvercote output
# ---------------------------------------------------------------------------

_WOL_UNSAFE = str.maketrans({c: " " for c in "(){}[],;=\"'"})


def _wol_name(s: str) -> str:
    """Sanitise a feature name for use as a Wolvercote label."""
    return s.translate(_WOL_UNSAFE).strip()


def _feature_to_wolvercote(f: MGEFeature) -> str:
    inner = ", ".join(_feature_to_wolvercote(c) for c in f.children)
    attrs = {}
    if f.family:
        attrs["family"] = f.family
    if f.tsd_seq:
        attrs["tsd"] = f.tsd_seq
    if f.element_type == "AMR":
        attrs["type"] = "AMR"
    attr_str = ""
    if attrs:
        pairs = ", ".join(f'{k}="{v}"' for k, v in attrs.items())
        attr_str = f"[{pairs}]"
    return f"{{{inner}}}{_wol_name(f.name)}{attr_str}"


def to_wolvercote(
    chromosome_features: list[MGEFeature],
    plasmid_features: list[list[MGEFeature]],
    sample_name: str = "",
) -> str:
    chr_inner = ", ".join(_feature_to_wolvercote(f) for f in chromosome_features)
    chr_str = f"({chr_inner}){sample_name}"

    plasmid_strs = []
    for i, pf in enumerate(plasmid_features):
        inner = ", ".join(_feature_to_wolvercote(f) for f in pf)
        label = f"plasmid_{i + 1}"
        plasmid_strs.append(f"{{{inner}}}{label}")

    parts = [chr_str] + plasmid_strs
    return ", ".join(parts)


# ---------------------------------------------------------------------------
# GFF3 output
# ---------------------------------------------------------------------------

def _gff_escape(s: str) -> str:
    return (
        s.replace("%", "%25")
         .replace(";", "%3B")
         .replace("=", "%3D")
         .replace("&", "%26")
         .replace(",", "%2C")
    )


def _feature_to_gff3_rows(
    f: MGEFeature,
    seqid: str,
    parent_id: str | None = None,
) -> list[str]:
    feat_id = f"{f.element_type}_{f.start}_{f.end}"
    attrs: list[str] = [
        f"ID={_gff_escape(feat_id)}",
        f"Name={_gff_escape(f.name)}",
        f"family={_gff_escape(f.family)}",
    ]
    if parent_id:
        attrs.append(f"Parent={_gff_escape(parent_id)}")
    if f.tsd_seq:
        attrs.append(f"tsd={f.tsd_seq}")
    if f.ir_left:
        attrs.append(f"ir_left={f.ir_left}")
    if f.ir_right:
        attrs.append(f"ir_right={f.ir_right}")
    if f.score is not None:
        attrs.append(f"score={f.score}")
    # Domain-specific attributes (seqid is meta — skip in output)
    for k, v in (f.attributes or {}).items():
        if k == "seqid" or v in (None, ""):
            continue
        attrs.append(f"{_gff_escape(k)}={_gff_escape(str(v))}")

    score_col = f"{f.score}" if f.score is not None else "."
    row = (
        f"{seqid}\tMatryoshka\t{f.element_type}\t{f.start}\t{f.end}\t"
        f"{score_col}\t{f.strand}\t.\t{';'.join(attrs)}"
    )
    rows = [row]
    for child in f.children:
        rows.extend(_feature_to_gff3_rows(child, seqid, feat_id))
    return rows


def to_gff3(features: list[MGEFeature], seqid: str = "sequence") -> str:
    rows = ["##gff-version 3"]
    for f in features:
        rows.extend(_feature_to_gff3_rows(f, seqid))
    return "\n".join(rows)


# ---------------------------------------------------------------------------
# JSON output
# ---------------------------------------------------------------------------

def _feature_to_dict(f: MGEFeature) -> dict[str, Any]:
    return {
        "element_type": f.element_type,
        "family": f.family,
        "name": f.name,
        "start": f.start,
        "end": f.end,
        "strand": f.strand,
        "tsd_length": f.tsd_length,
        "tsd_seq": f.tsd_seq,
        "ir_left": f.ir_left,
        "ir_right": f.ir_right,
        "score": f.score,
        "attributes": dict(f.attributes) if f.attributes else {},
        "children": [_feature_to_dict(c) for c in f.children],
    }


def to_json(features: list[MGEFeature], indent: int = 2) -> str:
    return json.dumps([_feature_to_dict(f) for f in features], indent=indent)


# ---------------------------------------------------------------------------
# SVG / PNG output via wolvercote renderer
# ---------------------------------------------------------------------------

def to_svg(features: list[MGEFeature], sample_name: str = "") -> str:
    from wolvercote import parse
    from wolvercote.renderer import render_svg
    wol = to_wolvercote(features, [], sample_name)
    cell_set = parse(wol)
    return render_svg(cell_set)


def to_png(features: list[MGEFeature], sample_name: str = "", dpi: int = 200) -> bytes:
    from wolvercote import parse
    from wolvercote.render_png import render_png
    wol = to_wolvercote(features, [], sample_name)
    cell_set = parse(wol)
    return render_png(cell_set, dpi=dpi)


# ---------------------------------------------------------------------------
# GenBank output (via Biopython)
# ---------------------------------------------------------------------------

def _feature_to_seqfeatures(f: MGEFeature) -> list:
    """Flatten a feature tree into a list of SeqFeature objects for GenBank."""
    from Bio.SeqFeature import SeqFeature, SimpleLocation

    strand = {"+": 1, "-": -1}.get(f.strand)
    loc = SimpleLocation(f.start - 1, f.end, strand=strand)
    qualifiers: dict[str, list[str]] = {
        "label": [f.name],
        "family": [f.family],
    }
    if f.tsd_seq:
        qualifiers["tsd"] = [f.tsd_seq]
    if f.ir_left:
        qualifiers["ir_left"] = [f.ir_left]
    if f.ir_right:
        qualifiers["ir_right"] = [f.ir_right]
    for k, v in (f.attributes or {}).items():
        if k == "seqid" or v in (None, ""):
            continue
        qualifiers[k] = [str(v)]

    feature_type_map = {
        "IS": "mobile_element",
        "transposon": "mobile_element",
        "integron": "mobile_element",
        "cassette": "CDS",
        "attC": "misc_feature",
        "integrase": "CDS",
        "AMR": "CDS",
        "replicon": "rep_origin",
    }
    ftype = feature_type_map.get(f.element_type, "misc_feature")
    out = [SeqFeature(location=loc, type=ftype, qualifiers=qualifiers)]
    for c in f.children:
        out.extend(_feature_to_seqfeatures(c))
    return out


def to_genbank(features: list[MGEFeature], seq: str, sample_name: str = "sequence") -> str:
    """Produce a GenBank flat-file string for one contig."""
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    rec = SeqRecord(Seq(seq), id=sample_name[:16] or "unknown",
                    name=sample_name[:16] or "unknown",
                    description=f"Matryoshka annotation of {sample_name}")
    rec.annotations["molecule_type"] = "DNA"
    for f in features:
        rec.features.extend(_feature_to_seqfeatures(f))

    buf = StringIO()
    from Bio import SeqIO as _SeqIO
    _SeqIO.write(rec, buf, "genbank")
    return buf.getvalue()
