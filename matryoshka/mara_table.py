"""MARA-style tabular SVG output generated from an MGE hierarchy."""

from __future__ import annotations

import html
import textwrap
from dataclasses import dataclass

from .detect import MGEFeature

WIDTH = 1280
MARGIN = 18
TITLE_H = 34
HEADER_H = 30
BASE_ROW_H = 28
KEY_H = 38
FONT = "Arial,Helvetica,sans-serif"
COLUMN_X = (MARGIN, 170, 520, 650, 760, WIDTH - MARGIN)
COLUMN_LABELS = ("Position", "Name*", "FID", "Type", "Notes")


@dataclass(frozen=True)
class TableRow:
    feature: MGEFeature
    depth: int
    position: str
    name: str
    fid: str
    feature_type: str
    notes: str


def _flatten(features: list[MGEFeature], depth: int = 0) -> list[tuple[MGEFeature, int]]:
    out: list[tuple[MGEFeature, int]] = []
    for feature in sorted(features, key=lambda item: (item.start, item.end, item.name)):
        out.append((feature, depth))
        out.extend(_flatten(feature.children, depth + 1))
    return out


def _feature_type(feature: MGEFeature) -> str:
    return {
        "transposon": "Tn",
        "transposition_unit": "TPU",
        "IS": "IS",
        "AMR": "R_gene",
        "gene": "gene",
        "res_site": "res",
        "IR": "IR",
        "integron": "region",
        "cassette": "cassette",
        "cassette_remnant": "cassette fragment",
        "cassette_array": "cassette array",
        "attC": "attC",
        "attI": "attI",
        "promoter": "promoter",
        "Pc_promoter": "promoter",
        "integrase": "gene",
        "transposition_gene": "gene",
        "integron_segment": "CS",
        "group_II_intron": "intron",
        "intron": "intron",
        "inserted_sequence": "inserted",
        "transposon_component": "region",
        "ISCR": "ISCR",
        "oriIS": "site",
        "terIS": "site",
        "captured_segment": "captured region",
        "DR": "DR",
        "direct_repeat": "DR",
        "replicon": "replicon",
        "oriV": "site",
        "oriT": "site",
        "ncRNA": "ncRNA",
        "plasmid_backbone": "backbone",
        "accessory_region": "accessory region",
        "multiresistance_region": "MRR",
        "annotation_gap": "gap",
        "unknown_fragment": "fragment",
    }.get(feature.element_type, feature.element_type)


def _number(value: object) -> str:
    try:
        number = float(value)  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return ""
    return f"{number:.2f}".rstrip("0").rstrip(".")


def _notes(feature: MGEFeature) -> str:
    attrs = feature.attributes
    if feature.element_type in {"transposon", "transposition_unit"}:
        parts: list[str] = []
        identity = _number(attrs.get("blast_identity"))
        coverage = _number(
            attrs.get("blast_subject_coverage") or attrs.get("blast_coverage")
        )
        if identity:
            parts.append(f"BLAST identity={identity}%")
        if coverage:
            parts.append(f"reference coverage={coverage}%")
        definition = attrs.get("closest_definition")
        subtype = attrs.get("defined_subtype")
        if definition:
            label = f"closest definition={definition}"
            if subtype and subtype != "canonical":
                label += f" (subtype {subtype})"
            parts.append(label)
        variant_status = str(attrs.get("variant_status", ""))
        structural_status = str(attrs.get("structural_status", ""))
        if variant_status and variant_status != "exact_reference":
            parts.append(variant_status.replace("_", " "))
        if structural_status and structural_status != "intact":
            parts.append(structural_status.replace("_", " "))
        inserted_bases = int(attrs.get("inserted_bases", 0) or 0)
        deleted_bases = int(attrs.get("deleted_bases", 0) or 0)
        mismatch_bases = int(attrs.get("mismatch_bases", 0) or 0)
        if mismatch_bases:
            parts.append(f"substitutions={mismatch_bases}")
        if inserted_bases:
            parts.append(f"inserted sequence≈{inserted_bases} bp")
        if deleted_bases:
            parts.append(f"deleted reference sequence≈{deleted_bases} bp")
        known_differences = attrs.get("known_differences_from_parent", [])
        if isinstance(known_differences, list) and known_differences:
            parts.append(
                "reviewed subtype differences: "
                + ", ".join(str(item) for item in known_differences)
            )
        assembly_status = attrs.get("component_assembly_status")
        if assembly_status:
            count = attrs.get("detected_component_count", 0)
            parts.append(
                f"component grammar={assembly_status}; "
                f"{count} independently sequence-detected components"
            )
        if feature.tsd_seq:
            strength = attrs.get("tsd_evidence_strength", "supporting")
            left = f"{attrs.get('tsd_left_start', '?')}..{attrs.get('tsd_left_end', '?')}"
            right = f"{attrs.get('tsd_right_start', '?')}..{attrs.get('tsd_right_end', '?')}"
            parts.append(
                f"boundary-adjacent DR={feature.tsd_seq} ({strength} evidence; "
                f"left={left}; right={right})"
            )
        elif feature.tsd_length:
            evidence = attrs.get("tsd_evidence")
            if evidence == "short_repeat_candidate":
                candidate = attrs.get("tsd_candidate_seq", "")
                parts.append(
                    f"expected DR={feature.tsd_length} bp; short-repeat "
                    f"candidate={candidate}; weak evidence"
                )
            elif evidence == "searched_not_found":
                parts.append(f"expected DR={feature.tsd_length} bp; searched but not found")
            elif evidence == "untestable_missing_flank":
                parts.append(
                    f"expected DR={feature.tsd_length} bp; flanking sequence unavailable"
                )
            else:
                parts.append(f"expected DR={feature.tsd_length} bp; not sequence-confirmed")
        if attrs.get("fragment"):
            parts.append("exact Tn1/Tn2/Tn3 identity unresolved")
        return "; ".join(parts)
    if attrs.get("source") == "tn123_component_scan":
        identity = attrs.get("blast_identity", "")
        coverage = attrs.get("blast_coverage", "")
        status = attrs.get("structural_status", "")
        detail = (
            f"independently sequence-detected {attrs.get('component_role', 'component')}; "
            f"identity={identity}%; coverage={coverage}%"
        )
        if status and status != "intact":
            detail += f"; {str(status).replace('_', ' ')}"
        return detail
    if feature.element_type == "IR":
        return f"terminal inverted repeat; IR={feature.end - feature.start + 1} bp"
    if feature.element_type == "res_site":
        return f"curated {feature.family} resolution site"
    if feature.element_type == "inserted_sequence":
        inserted = int(feature.attributes.get("inserted_bases", 0) or 0)
        return f"unclassified insertion; estimated inserted sequence={inserted} bp"
    if feature.element_type in {"DR", "direct_repeat"}:
        sequence = attrs.get("sequence") or feature.tsd_seq
        return f"confirmed direct repeat={sequence}" if sequence else "direct repeat; sequence unresolved"
    if feature.element_type == "cassette_remnant":
        return "incomplete cassette; not promoted to a complete gene cassette"
    if feature.element_type in {"multiresistance_region", "accessory_region"}:
        return "context region; independent mobility of the whole region is not implied"
    if attrs.get("source") == "curated_reference":
        accession = attrs.get("source_accession", "")
        return f"curated internal feature from {accession}".rstrip()
    return str(attrs.get("note", ""))


def _row(feature: MGEFeature, depth: int) -> TableRow:
    attrs = feature.attributes
    fid = str(
        attrs.get("fid")
        or attrs.get("source_accession")
        or attrs.get("reference_id")
        or ""
    )
    start = int(attrs.get("absolute_start", feature.start))
    end = int(attrs.get("absolute_end", feature.end))
    return TableRow(
        feature=feature,
        depth=depth,
        position=f"{start}..{end}",
        name=feature.name,
        fid=fid,
        feature_type=_feature_type(feature),
        notes=_notes(feature),
    )


def _text(
    x: float,
    y: float,
    value: str,
    *,
    size: int = 11,
    weight: str = "normal",
    anchor: str = "start",
    fill: str = "#171717",
) -> str:
    return (
        f'<text x="{x:.1f}" y="{y:.1f}" font-size="{size}" '
        f'font-weight="{weight}" text-anchor="{anchor}" fill="{fill}">'
        f"{html.escape(value)}</text>"
    )


def _arrow(x: float, y: float, strand: str, *, projected: bool = False) -> str:
    if strand not in {"+", "-"}:
        return ""
    length = 34
    if strand == "+":
        x1, x2 = x, x + length
        head = f"{x2:.1f},{y:.1f} {x2 - 8:.1f},{y - 5:.1f} {x2 - 8:.1f},{y + 5:.1f}"
    else:
        x1, x2 = x + length, x
        head = f"{x2:.1f},{y:.1f} {x2 + 8:.1f},{y - 5:.1f} {x2 + 8:.1f},{y + 5:.1f}"
    colour = "#ffffff" if projected else "#000000"
    dash = ' stroke-dasharray="4,2"' if projected else ""
    return (
        f'<line x1="{x1:.1f}" y1="{y:.1f}" x2="{x2:.1f}" y2="{y:.1f}" '
        f'stroke="#000" stroke-width="6"{dash}/>'
        f'<polygon points="{head}" fill="{colour}" stroke="#000"{dash}/>'
    )


def _ir_triangle(x: float, y: float, name: str) -> str:
    direction = 1 if name == "IRL" else -1
    return (
        f'<polygon points="{x:.1f},{y - 6:.1f} {x + direction * 10:.1f},{y:.1f} '
        f'{x:.1f},{y + 6:.1f}" fill="#fff" stroke="#111" stroke-width="1.2"/>'
    )


def _wrapped_notes(notes: str) -> list[str]:
    if not notes:
        return [""]
    return textwrap.wrap(notes, width=72, break_long_words=False, break_on_hyphens=False)


def to_mara_table_svg(roots: list[MGEFeature], sample_name: str = "") -> str:
    """Render a MARA-style annotation table as selectable-text SVG."""
    rows = [_row(feature, depth) for feature, depth in _flatten(roots)]
    note_lines = [_wrapped_notes(row.notes) for row in rows]
    row_heights = [max(BASE_ROW_H, 12 + len(lines) * 14) for lines in note_lines]
    table_bottom = TITLE_H + HEADER_H + sum(row_heights)
    height = table_bottom + KEY_H + MARGIN

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{WIDTH}" height="{height}" '
        f'style="background:#fff;font-family:{FONT}">',
        f'<rect x="0" y="0" width="{WIDTH}" height="{height}" fill="#fff"/>',
        _text(MARGIN, 23, sample_name or "MARA annotation table", size=14, weight="bold"),
        f'<rect x="{MARGIN}" y="{TITLE_H}" width="{WIDTH - 2 * MARGIN}" '
        f'height="{HEADER_H}" fill="#eeeeec" stroke="#c9c9c6"/>',
    ]

    for index, label in enumerate(COLUMN_LABELS):
        parts.append(_text(COLUMN_X[index] + 8, TITLE_H + 20, label, size=11, weight="bold"))

    y = TITLE_H + HEADER_H
    for index, (row, lines, row_height) in enumerate(zip(rows, note_lines, row_heights, strict=True)):
        fill = "#f2f1ef" if index % 2 == 0 else "#ffffff"
        if row.feature.element_type in {"transposon", "transposition_unit"}:
            fill = "#e8e8e5"
        parts.append(
            f'<rect x="{MARGIN}" y="{y}" width="{WIDTH - 2 * MARGIN}" '
            f'height="{row_height}" fill="{fill}" stroke="#ddddda" stroke-width="0.7"/>'
        )
        baseline = y + 18
        weight = (
            "bold"
            if row.feature.element_type in {"transposon", "transposition_unit"}
            else "normal"
        )
        parts.append(_text(COLUMN_X[0] + 8, baseline, row.position, size=10))
        name_x = COLUMN_X[1] + 8 + row.depth * 20
        parts.append(_text(name_x, baseline, row.name, size=10, weight=weight))
        glyph_x = min(name_x + len(row.name) * 6 + 12, COLUMN_X[2] - 48)
        if row.feature.element_type == "IR":
            parts.append(_ir_triangle(glyph_x, baseline - 3, row.name))
        else:
            parts.append(_arrow(
                glyph_x,
                baseline - 3,
                row.feature.strand,
                projected=row.feature.attributes.get("evidence_class")
                == "reference_projected",
            ))
        parts.append(_text(COLUMN_X[2] + 8, baseline, row.fid, size=10))
        parts.append(_text(COLUMN_X[3] + 8, baseline, row.feature_type, size=10))
        for line_index, line in enumerate(lines):
            parts.append(_text(COLUMN_X[4] + 8, baseline + line_index * 14, line, size=10))
        y += row_height

    for x in COLUMN_X:
        parts.append(
            f'<line x1="{x}" y1="{TITLE_H}" x2="{x}" y2="{table_bottom}" '
            'stroke="#d2d2cf" stroke-width="0.8"/>'
        )
    key_y = table_bottom + 23
    parts.extend([
        _text(MARGIN, key_y, "Key:", size=10, weight="bold"),
        _arrow(62, key_y - 3, "+"),
        _text(105, key_y, "strand", size=10),
        _ir_triangle(180, key_y - 3, "IRL"),
        _text(198, key_y, "terminal inverted repeat (IR)", size=10),
        f'<circle cx="430" cy="{key_y - 5}" r="5" fill="#9ca3aa" stroke="#111"/>',
        _text(443, key_y, "sequence-matched DR/TSD", size=10),
        f'<circle cx="590" cy="{key_y - 5}" r="5" fill="#fff" stroke="#737b83"/>',
        _text(603, key_y, "expected, unconfirmed DR/TSD", size=10),
        _text(
            850,
            key_y,
            "Solid components are sequence-detected; dashed/outlined components are reference-projected.",
            size=10,
            fill="#555555",
        ),
    ])
    parts.append("</svg>")
    return "\n".join(parts)
