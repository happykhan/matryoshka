"""
output.py — Write Matryoshka results in multiple formats.

Formats:
  - Wolvercote (cell format) — compact nested notation
  - GFF3
  - JSON hierarchy
"""

from __future__ import annotations
import json
import sys
from pathlib import Path

from .detect import MGEFeature


# ---------------------------------------------------------------------------
# Wolvercote output
# ---------------------------------------------------------------------------

def _feature_to_wolvercote(f: MGEFeature) -> str:
    """Recursively convert an MGEFeature to a Wolvercote { } block."""
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
    return f"{{{inner}}}{f.name}{attr_str}"


def to_wolvercote(
    chromosome_features: list[MGEFeature],
    plasmid_features: list[list[MGEFeature]],
    sample_name: str = "",
) -> str:
    """
    Build a Wolvercote string from the hierarchy.

    chromosome_features: root-level features on the chromosome
    plasmid_features: list of root-level feature lists, one per plasmid
    """
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

def _feature_to_gff3_rows(f: MGEFeature, seqid: str, parent_id: str | None = None) -> list[str]:
    feat_id = f"{f.element_type}_{f.start}_{f.end}"
    attrs = [f"ID={feat_id}", f"Name={f.name}", f"family={f.family}"]
    if parent_id:
        attrs.append(f"Parent={parent_id}")
    if f.tsd_seq:
        attrs.append(f"tsd={f.tsd_seq}")
    attr_str = ";".join(attrs)
    row = f"{seqid}\tMatryoshka\t{f.element_type}\t{f.start}\t{f.end}\t.\t{f.strand}\t.\t{attr_str}"
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

def _feature_to_dict(f: MGEFeature) -> dict:
    return {
        "element_type": f.element_type,
        "family": f.family,
        "name": f.name,
        "start": f.start,
        "end": f.end,
        "strand": f.strand,
        "tsd_seq": f.tsd_seq,
        "ir_left": f.ir_left,
        "ir_right": f.ir_right,
        "children": [_feature_to_dict(c) for c in f.children],
    }


def to_json(features: list[MGEFeature], indent: int = 2) -> str:
    return json.dumps([_feature_to_dict(f) for f in features], indent=indent)
