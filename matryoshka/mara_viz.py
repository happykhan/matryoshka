"""Render annotations using MARA's compact single-line visual grammar.

This renderer is intentionally separate from :mod:`matryoshka.viz`. The
existing nested linear visualisation remains available as ``--format linear``.
"""

from __future__ import annotations

import html

from .detect import MGEFeature

WIDTH = 1280
HEIGHT = 272
LEFT = 72
RIGHT = 72
TRACK_Y = 66
BLOCK_TOP = 44
BLOCK_BOTTOM = 84
FONT = "Arial,Helvetica,sans-serif"
TN123_COLOUR = "#009b55"
OUTLINE = "#111111"
GAP_COLOUR = "#e3262e"
CASSETTE_COLOUR = "#68c5e3"
INTEGRON_COLOUR = "#ed873b"
INTRON_COLOUR = "#b8b8b8"
ISCR_COLOUR = "#fff0ad"
REPLICON_COLOUR = "#50a6a6"
NCRNA_COLOUR = "#7aa6d8"
REGION_COLOUR = "#777777"
CAPTURE_COLOUR = "#e4e4df"
SUPPORTED_RENDER_SYMBOLS = frozenset({
    "att_site",
    "captured_region",
    "cassette_box",
    "cassette_fragment",
    "conserved_segment",
    "dashed_gap",
    "dr_lollipop",
    "fragment_block",
    "gene_arrow",
    "intron_block_arrow",
    "ir_flag",
    "is_block_arrow",
    "iscr_block_arrow",
    "jagged_insertion",
    "ncrna_block",
    "promoter_flag",
    "region_outline",
    "replicon_block",
    "res_site",
    "site_marker",
})
TRANSPOSON_COLOURS = {
    "Tn21": "#9f1d2d",
    "Tn1696": "#758895",
    "Tn1721": "#d74747",
    "Tn1722": "#d74747",
    "Tn402": "#9b4b9d",
    "Tn4401": "#ef92dc",
    "Tn5393": "#efce53",
    "Tn5403": "#eceb83",
}


def _walk(features: list[MGEFeature]) -> list[MGEFeature]:
    out: list[MGEFeature] = []
    for feature in features:
        out.append(feature)
        out.extend(_walk(feature.children))
    return out


def _x(position: int, seq_len: int) -> float:
    drawing_width = WIDTH - LEFT - RIGHT
    return LEFT + drawing_width * (max(1, position) - 1) / max(seq_len - 1, 1)


def _text(
    x: float,
    y: float,
    value: str,
    *,
    size: int = 12,
    anchor: str = "middle",
    weight: str = "normal",
    style: str = "normal",
    fill: str = OUTLINE,
) -> str:
    return (
        f'<text x="{x:.1f}" y="{y:.1f}" text-anchor="{anchor}" '
        f'font-size="{size}" font-weight="{weight}" font-style="{style}" '
        f'fill="{fill}">{html.escape(value)}</text>'
    )


def _transposon_colour(feature: MGEFeature) -> str:
    return TRANSPOSON_COLOURS.get(feature.family, TN123_COLOUR)


def _solid_block(feature: MGEFeature, seq_len: int) -> str:
    x1, x2 = _x(feature.start, seq_len), _x(feature.end, seq_len)
    return (
        f'<rect x="{x1:.1f}" y="{BLOCK_TOP}" width="{max(x2 - x1, 2):.1f}" '
        f'height="{BLOCK_BOTTOM - BLOCK_TOP}" fill="{_transposon_colour(feature)}" '
        f'stroke="{OUTLINE}" stroke-width="2"/>'
    )


def _fragment_block(feature: MGEFeature, seq_len: int) -> str:
    """MARA-style side-specific jagged ends for a partial feature."""
    x1, x2 = _x(feature.start, seq_len), _x(feature.end, seq_len)
    tooth = min(8.0, max((x2 - x1) / 8, 3.0))
    jagged_left = not bool(feature.attributes.get("left_end_covered", False))
    jagged_right = not bool(feature.attributes.get("right_end_covered", False))
    points: list[tuple[float, float]] = [(x1 + tooth if jagged_left else x1, BLOCK_TOP)]
    points.append((x2 - tooth if jagged_right else x2, BLOCK_TOP))
    if jagged_right:
        points.extend([
            (x2, BLOCK_TOP + 8),
            (x2 - tooth, BLOCK_TOP + 16),
            (x2, BLOCK_TOP + 25),
            (x2 - tooth, BLOCK_TOP + 34),
            (x2, BLOCK_BOTTOM - 8),
        ])
    points.append((x2 - tooth if jagged_right else x2, BLOCK_BOTTOM))
    points.append((x1 + tooth if jagged_left else x1, BLOCK_BOTTOM))
    if jagged_left:
        points.extend([
            (x1, BLOCK_BOTTOM - 8),
            (x1 + tooth, BLOCK_BOTTOM - 16),
            (x1, BLOCK_TOP + 25),
            (x1 + tooth, BLOCK_TOP + 16),
            (x1, BLOCK_TOP + 8),
        ])
    encoded = " ".join(f"{x:.1f},{y:.1f}" for x, y in points)
    return (
        f'<polygon points="{encoded}" fill="{_transposon_colour(feature)}" '
        f'stroke="{OUTLINE}" stroke-width="2"/>'
    )


def _block_arrow(
    feature: MGEFeature,
    seq_len: int,
    fill: str,
    *,
    stroke_width: float = 1.5,
) -> str:
    """MARA block arrow; for IS the pointed end denotes IRR."""
    x1, x2 = _x(feature.start, seq_len), _x(feature.end, seq_len)
    width = max(x2 - x1, 4.0)
    tip = min(16.0, width * 0.28)
    top, bottom = BLOCK_TOP + 4, BLOCK_BOTTOM - 4
    if feature.strand == "-":
        points = (
            f"{x1:.1f},{(top + bottom) / 2:.1f} {x1 + tip:.1f},{top:.1f} "
            f"{x2:.1f},{top:.1f} {x2:.1f},{bottom:.1f} {x1 + tip:.1f},{bottom:.1f}"
        )
    elif feature.strand == "+":
        points = (
            f"{x1:.1f},{top:.1f} {x2 - tip:.1f},{top:.1f} "
            f"{x2:.1f},{(top + bottom) / 2:.1f} {x2 - tip:.1f},{bottom:.1f} "
            f"{x1:.1f},{bottom:.1f}"
        )
    else:
        points = f"{x1:.1f},{top:.1f} {x2:.1f},{top:.1f} {x2:.1f},{bottom:.1f} {x1:.1f},{bottom:.1f}"
    return (
        f'<polygon points="{points}" fill="{fill}" stroke="{OUTLINE}" '
        f'stroke-width="{stroke_width}"/>'
    )


def _flat_feature_block(feature: MGEFeature, seq_len: int, fill: str) -> str:
    x1, x2 = _x(feature.start, seq_len), _x(feature.end, seq_len)
    return (
        f'<rect x="{x1:.1f}" y="{BLOCK_TOP + 4}" width="{max(x2 - x1, 2):.1f}" '
        f'height="{BLOCK_BOTTOM - BLOCK_TOP - 8}" fill="{fill}" '
        f'stroke="{OUTLINE}" stroke-width="1.2"/>'
    )


def _attc_marker(feature: MGEFeature, seq_len: int) -> str:
    x = (_x(feature.start, seq_len) + _x(feature.end, seq_len)) / 2
    return (
        f'<line x1="{x:.1f}" y1="{BLOCK_TOP + 2}" x2="{x:.1f}" '
        f'y2="{BLOCK_BOTTOM - 2}" stroke="#5f666d" stroke-width="2"/>'
    )


def _site_marker(feature: MGEFeature, seq_len: int, colour: str) -> str:
    """Compact site marker for attI, oriIS, terIS, oriV and oriT."""
    x = (_x(feature.start, seq_len) + _x(feature.end, seq_len)) / 2
    return (
        f'<line x1="{x:.1f}" y1="{BLOCK_TOP + 1}" x2="{x:.1f}" '
        f'y2="{BLOCK_BOTTOM - 1}" stroke="{colour}" stroke-width="2.5"/>'
        f'<polygon points="{x:.1f},{BLOCK_TOP - 2} {x + 5:.1f},{BLOCK_TOP + 5} '
        f'{x - 5:.1f},{BLOCK_TOP + 5}" fill="{colour}" stroke="{OUTLINE}" '
        'stroke-width="0.7"/>'
    )


def _promoter_flag(feature: MGEFeature, seq_len: int) -> str:
    x = (_x(feature.start, seq_len) + _x(feature.end, seq_len)) / 2
    direction = -1 if feature.strand == "-" else 1
    return (
        f'<line x1="{x:.1f}" y1="{BLOCK_BOTTOM - 2}" x2="{x:.1f}" '
        f'y2="{BLOCK_TOP + 2}" stroke="{INTEGRON_COLOUR}" stroke-width="2"/>'
        f'<line x1="{x:.1f}" y1="{BLOCK_TOP + 2}" '
        f'x2="{x + direction * 11:.1f}" y2="{BLOCK_TOP + 2}" '
        f'stroke="{INTEGRON_COLOUR}" stroke-width="2"/>'
    )


def _region_outline(
    feature: MGEFeature,
    seq_len: int,
    colour: str,
    *,
    y_offset: int = 0,
    dashed: bool = False,
) -> str:
    x1, x2 = _x(feature.start, seq_len), _x(feature.end, seq_len)
    dash = ' stroke-dasharray="6,3"' if dashed else ""
    return (
        f'<rect x="{x1:.1f}" y="{BLOCK_TOP + y_offset}" '
        f'width="{max(x2 - x1, 2):.1f}" '
        f'height="{BLOCK_BOTTOM - BLOCK_TOP - 2 * y_offset}" fill="none" '
        f'stroke="{colour}" stroke-width="2"{dash}/>'
    )


def _gene_arrow(feature: MGEFeature, seq_len: int) -> str:
    x1, x2 = _x(feature.start, seq_len), _x(feature.end, seq_len)
    width = max(x2 - x1, 4.0)
    tip = min(15.0, width * 0.35)
    top, middle, bottom = TRACK_Y - 8, TRACK_Y, TRACK_Y + 8
    if feature.strand == "-":
        points = (
            f"{x1:.1f},{middle} {x1 + tip:.1f},{top} {x2:.1f},{top} "
            f"{x2:.1f},{bottom} {x1 + tip:.1f},{bottom}"
        )
    elif feature.strand == "+":
        points = (
            f"{x1:.1f},{top} {x2 - tip:.1f},{top} {x2:.1f},{middle} "
            f"{x2 - tip:.1f},{bottom} {x1:.1f},{bottom}"
        )
    else:
        points = (
            f"{x1:.1f},{top} {x2:.1f},{top} {x2:.1f},{bottom} {x1:.1f},{bottom}"
        )
    return f'<polygon points="{points}" fill="#000" stroke="#000" stroke-width="1"/>'


def _res_marker(feature: MGEFeature, seq_len: int) -> str:
    x1, x2 = _x(feature.start, seq_len), _x(feature.end, seq_len)
    return (
        f'<line x1="{x1:.1f}" y1="{BLOCK_TOP + 3}" x2="{x1:.1f}" '
        f'y2="{BLOCK_BOTTOM - 3}" stroke="{OUTLINE}" stroke-width="1" '
        'stroke-dasharray="2,2"/>'
        f'<line x1="{x2:.1f}" y1="{BLOCK_TOP + 3}" x2="{x2:.1f}" '
        f'y2="{BLOCK_BOTTOM - 3}" stroke="{OUTLINE}" stroke-width="1" '
        'stroke-dasharray="2,2"/>'
    )


def _insertion_block(feature: MGEFeature, seq_len: int) -> str:
    x1, x2 = _x(feature.start, seq_len), _x(feature.end, seq_len)
    tooth = min(6.0, max((x2 - x1) / 6, 2.0))
    points = [
        (x1 + tooth, BLOCK_TOP),
        (x2 - tooth, BLOCK_TOP),
        (x2, BLOCK_TOP + 7),
        (x2 - tooth, BLOCK_TOP + 14),
        (x2, BLOCK_TOP + 21),
        (x2 - tooth, BLOCK_BOTTOM),
        (x1 + tooth, BLOCK_BOTTOM),
        (x1, BLOCK_BOTTOM - 7),
        (x1 + tooth, BLOCK_BOTTOM - 14),
        (x1, BLOCK_TOP + 19),
        (x1 + tooth, BLOCK_TOP + 10),
    ]
    encoded = " ".join(f"{x:.1f},{y:.1f}" for x, y in points)
    return (
        f'<polygon points="{encoded}" fill="#ffffff" stroke="{OUTLINE}" '
        'stroke-width="1.5"/>'
    )


def _ir_flag(feature: MGEFeature, seq_len: int) -> str:
    left_end = feature.name in {"IRL", "IRi"}
    boundary = feature.start if left_end else feature.end
    x = _x(boundary, seq_len)
    inward = 1 if left_end else -1
    return (
        f'<line x1="{x:.1f}" y1="{BLOCK_TOP}" x2="{x:.1f}" y2="21" '
        f'stroke="{OUTLINE}" stroke-width="1.5"/>'
        f'<polygon points="{x:.1f},23 {x + inward * 10:.1f},29 {x:.1f},35" '
        f'fill="#ffffff" stroke="{OUTLINE}" stroke-width="1.2"/>'
    )


def _lollipop(
    x: float,
    sequence: str | None,
    expected_length: int | None = None,
    colour: str = "#9ca3aa",
    evidence: str = "",
) -> str:
    confirmed = sequence is not None
    fill = colour if confirmed else "#ffffff"
    stroke = OUTLINE if confirmed else "#737b83"
    if confirmed:
        title = f"confirmed TSD: {sequence}"
    elif evidence == "searched_not_found":
        title = f"expected {expected_length or '?'} bp TSD; searched but not found"
    elif evidence == "untestable_missing_flank":
        title = f"expected {expected_length or '?'} bp TSD; flanking sequence unavailable"
    else:
        title = f"expected {expected_length or '?'} bp TSD; not sequence-confirmed"
    label = _text(x, 9, sequence, size=8) if sequence else ""
    return (
        f'<g><title>{html.escape(title)}</title>'
        f'<line x1="{x:.1f}" y1="{BLOCK_TOP}" x2="{x:.1f}" y2="18" '
        f'stroke="{stroke}" stroke-width="1.4"/>'
        f'<circle cx="{x:.1f}" cy="13" r="5" fill="{fill}" '
        f'stroke="{stroke}" stroke-width="1.2"/>{label}</g>'
    )


def _legend() -> str:
    """Self-contained key for readers who have not seen MARA notation."""
    top = 164
    row1 = 194
    row2 = 222
    parts = [
        f'<line x1="{LEFT}" y1="{top}" x2="{WIDTH - RIGHT}" y2="{top}" '
        'stroke="#d4d4d4" stroke-width="1"/>',
        _text(LEFT, top + 18, "Key", size=10, anchor="start", weight="bold"),
        # Directional gene arrow.
        f'<line x1="115" y1="{row1}" x2="145" y2="{row1}" stroke="#000" stroke-width="6"/>',
        f'<polygon points="145,{row1} 137,{row1 - 5} 137,{row1 + 5}" fill="#000"/>',
        _text(154, row1 + 3, "gene (arrow shows strand)", size=9, anchor="start"),
        # Terminal inverted-repeat flag.
        f'<line x1="350" y1="{row1 + 7}" x2="350" y2="{row1 - 9}" stroke="#111" stroke-width="1.3"/>',
        f'<polygon points="350,{row1 - 7} 360,{row1 - 1} 350,{row1 + 5}" fill="#fff" stroke="#111"/>',
        _text(370, row1 + 3, "terminal inverted repeat (IR)", size=9, anchor="start"),
        # Confirmed and expected/unconfirmed target-site duplications.
        f'<line x1="590" y1="{row1 + 7}" x2="590" y2="{row1 - 7}" stroke="#111"/>',
        f'<circle cx="590" cy="{row1 - 8}" r="5" fill="#9ca3aa" stroke="#111"/>',
        _text(603, row1 + 3, "confirmed DR/TSD", size=9, anchor="start"),
        f'<line x1="744" y1="{row1 + 7}" x2="744" y2="{row1 - 7}" stroke="#737b83"/>',
        f'<circle cx="744" cy="{row1 - 8}" r="5" fill="#fff" stroke="#737b83"/>',
        _text(757, row1 + 3, "expected, unconfirmed DR/TSD", size=9, anchor="start"),
        # White IS block arrow and coloured mobile-element block.
        f'<polygon points="996,{row1 - 9} 1022,{row1 - 9} 1032,{row1} 1022,{row1 + 9} 996,{row1 + 9}" fill="#fff" stroke="#111"/>',
        _text(1040, row1 + 3, "insertion sequence (IS)", size=9, anchor="start"),
        f'<rect x="94" y="{row2 - 9}" width="42" height="18" fill="{TN123_COLOUR}" stroke="#111"/>',
        _text(145, row2 + 3, "mobile element / transposon", size=9, anchor="start"),
        f'<rect x="350" y="{row2 - 9}" width="42" height="18" fill="{CASSETTE_COLOUR}" stroke="#111"/>',
        _text(401, row2 + 3, "gene cassette", size=9, anchor="start"),
        f'<rect x="520" y="{row2 - 9}" width="42" height="18" fill="{INTEGRON_COLOUR}" stroke="#111"/>',
        _text(571, row2 + 3, "integron conserved segment", size=9, anchor="start"),
        f'<polygon points="770,{row2 - 9} 800,{row2 - 9} 806,{row2 - 3} 800,{row2 + 3} 806,{row2 + 9} 770,{row2 + 9} 764,{row2 + 3} 770,{row2 - 3} 764,{row2 - 9}" fill="#fff" stroke="#111"/>',
        _text(815, row2 + 3, "inserted / unresolved sequence", size=9, anchor="start"),
        f'<line x1="1045" y1="{row2}" x2="1085" y2="{row2}" stroke="{GAP_COLOUR}" stroke-width="2" stroke-dasharray="6,4"/>',
        _text(1094, row2 + 3, "unannotated flank", size=9, anchor="start"),
        _text(
            WIDTH / 2,
            252,
            "Symbols are annotation calls; see the matching table Notes for sequence-detected versus reference-projected evidence.",
            size=9,
            fill="#555555",
        ),
    ]
    return "".join(parts)


def _gap_line(start: int, end: int, seq_len: int) -> list[str]:
    length = end - start + 1
    x1, x2 = _x(start, seq_len), _x(end, seq_len)
    return [
        f'<line x1="{x1:.1f}" y1="{TRACK_Y}" x2="{x2:.1f}" y2="{TRACK_Y}" '
        f'stroke="{GAP_COLOUR}" stroke-width="2" stroke-dasharray="6,4"/>',
        _text((x1 + x2) / 2, TRACK_Y + 25, str(length), size=9, fill=GAP_COLOUR),
    ]


def _long_unannotated_gaps(
    visible: list[MGEFeature],
    transposons: list[MGEFeature],
    seq_len: int,
) -> list[str]:
    """Show only gaps not already represented by a complete transposon block."""
    # Work from the most specific visible spans. Container outlines overlap
    # their children and must not generate false backwards/giant gaps.
    leaves = [
        feature
        for feature in visible
        if not any(
            other is not feature
            and feature.start <= other.start
            and other.end <= feature.end
            and (feature.start, feature.end) != (other.start, other.end)
            for other in visible
        )
    ]
    intervals: list[list[int]] = []
    for feature in sorted(leaves, key=lambda item: (item.start, item.end)):
        if intervals and feature.start <= intervals[-1][1] + 1:
            intervals[-1][1] = max(intervals[-1][1], feature.end)
        else:
            intervals.append([feature.start, feature.end])
    parts: list[str] = []
    for left, right in zip(intervals, intervals[1:], strict=False):
        start, end = left[1] + 1, right[0] - 1
        if end - start + 1 <= 50:
            continue
        contained = any(t.start <= start and end <= t.end for t in transposons)
        if not contained:
            parts.extend(_gap_line(start, end, seq_len))
    return parts


def _edge_gaps(visible: list[MGEFeature], seq_len: int) -> list[str]:
    if not visible:
        return []
    parts: list[str] = []
    first = min(feature.start for feature in visible)
    last = max(feature.end for feature in visible)
    if first > 51:
        parts.extend(_gap_line(1, first - 1, seq_len))
    if seq_len - last > 50:
        parts.extend(_gap_line(last + 1, seq_len, seq_len))
    return parts


def to_mara_svg(
    roots: list[MGEFeature],
    seq_len: int,
    sample_name: str = "",
) -> str:
    """Render features on one line using the grammar of MARA Figures 1–2."""
    all_features = _walk(roots)
    transposons = [feature for feature in all_features if feature.element_type == "transposon"]
    transposition_units = [
        feature
        for feature in all_features
        if feature.element_type
        in {
            "transposition_unit",
            "composite_transposon",
            "IS26_translocatable_unit",
            "ISCR_capture_unit",
        }
    ]
    major_units = transposons + transposition_units
    genes = [
        feature
        for feature in all_features
        if feature.element_type
        in {
            "gene",
            "AMR",
            "integrase",
            "transposition_gene",
            "plasmid_backbone_gene",
            "plasmid_function_gene",
        }
    ]
    insertion_sequences = [
        feature for feature in all_features if feature.element_type == "IS"
    ]
    iscrs = [feature for feature in all_features if feature.element_type == "ISCR"]
    introns = [
        feature
        for feature in all_features
        if feature.element_type in {"group_II_intron", "intron"}
    ]
    integrons = [
        feature for feature in all_features if feature.element_type == "integron"
    ]
    conserved_segments = [
        feature
        for feature in all_features
        if feature.element_type == "integron_segment"
        or feature.name in {"5'-CS", "3'-CS", "5-CS", "3-CS"}
    ]
    cassettes = [
        feature
        for feature in all_features
        if feature.element_type in {"cassette", "cassette_remnant"}
    ]
    cassette_arrays = [
        feature for feature in all_features if feature.element_type == "cassette_array"
    ]
    attc_sites = [feature for feature in all_features if feature.element_type == "attC"]
    atti_sites = [feature for feature in all_features if feature.element_type == "attI"]
    promoters = [
        feature for feature in all_features if feature.element_type in {"promoter", "Pc_promoter"}
    ]
    other_sites = [
        feature
        for feature in all_features
        if feature.element_type in {"oriIS", "terIS", "oriV", "oriT", "recombination_site"}
    ]
    res_sites = [feature for feature in all_features if feature.element_type == "res_site"]
    irs = [feature for feature in all_features if feature.element_type == "IR"]
    direct_repeats = [
        feature for feature in all_features if feature.element_type in {"DR", "direct_repeat"}
    ]
    insertions = [
        feature for feature in all_features if feature.element_type == "inserted_sequence"
    ]
    components = [
        feature for feature in all_features if feature.element_type == "transposon_component"
    ]
    captured_segments = [
        feature for feature in all_features if feature.element_type == "captured_segment"
    ]
    unknown_fragments = [
        feature for feature in all_features if feature.element_type == "unknown_fragment"
    ]
    replicons = [feature for feature in all_features if feature.element_type == "replicon"]
    ncrnas = [feature for feature in all_features if feature.element_type == "ncRNA"]
    context_regions = [
        feature
        for feature in all_features
        if feature.element_type
        in {"plasmid", "plasmid_backbone", "accessory_region", "multiresistance_region"}
    ]
    explicit_gaps = [
        feature for feature in all_features if feature.element_type in {"annotation_gap", "unknown_gap"}
    ]

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{WIDTH}" height="{HEIGHT}" '
        f'style="background:#fff;font-family:{FONT}">',
        f'<rect x="0" y="0" width="{WIDTH}" height="{HEIGHT}" fill="#ffffff"/>',
    ]

    for feature in major_units:
        parts.append(
            _fragment_block(feature, seq_len)
            if feature.attributes.get("fragment")
            else _solid_block(feature, seq_len)
        )
        mid = (_x(feature.start, seq_len) + _x(feature.end, seq_len)) / 2
        label = (
            "Tn1/2/3 fragment"
            if feature.attributes.get("fragment") and feature.family == "Tn3_family"
            else feature.name
        )
        # Keep the unit name above any internal IR flag at the locus midpoint
        # (notably the duplicated IRR in Tn1721).
        parts.append(_text(mid, 20, label, size=12, weight="bold"))
        if feature.attributes.get("fragment"):
            x1, x2 = _x(feature.start, seq_len), _x(feature.end, seq_len)
            parts.extend([
                f'<line x1="{x1 - 20:.1f}" y1="{TRACK_Y}" x2="{x1:.1f}" '
                f'y2="{TRACK_Y}" stroke="{GAP_COLOUR}" stroke-width="2" '
                'stroke-dasharray="5,3"/>',
                f'<line x1="{x2:.1f}" y1="{TRACK_Y}" x2="{x2 + 20:.1f}" '
                f'y2="{TRACK_Y}" stroke="{GAP_COLOUR}" stroke-width="2" '
                'stroke-dasharray="5,3"/>',
            ])

    # The outer plasmid/MRR layer is contextual: it groups independently
    # annotated modules but never implies that the complete region is mobile.
    for feature in context_regions:
        parts.append(_region_outline(feature, seq_len, REGION_COLOUR, y_offset=8, dashed=True))
        mid = (_x(feature.start, seq_len) + _x(feature.end, seq_len)) / 2
        parts.append(_text(mid, 133, feature.name, size=8, fill=REGION_COLOUR))

    # IntegronFinder's parent span is shown as a restrained orange outline;
    # explicit conserved segments and cassette units are then drawn on top.
    for feature in integrons:
        x1, x2 = _x(feature.start, seq_len), _x(feature.end, seq_len)
        parts.append(
            f'<rect x="{x1:.1f}" y="{BLOCK_TOP + 1}" width="{max(x2 - x1, 2):.1f}" '
            f'height="{BLOCK_BOTTOM - BLOCK_TOP - 2}" fill="none" '
            f'stroke="{INTEGRON_COLOUR}" stroke-width="3"/>'
        )

    for feature in cassette_arrays:
        parts.append(_region_outline(feature, seq_len, CASSETTE_COLOUR, y_offset=5))

    for feature in conserved_segments:
        parts.append(_flat_feature_block(feature, seq_len, INTEGRON_COLOUR))
        mid = (_x(feature.start, seq_len) + _x(feature.end, seq_len)) / 2
        suffix = "#" if feature.attributes.get("fragment") else ""
        parts.append(_text(mid, 103, f"{feature.name}{suffix}", size=8, weight="bold"))

    for feature in components:
        parts.append(_flat_feature_block(feature, seq_len, _transposon_colour(feature)))
        mid = (_x(feature.start, seq_len) + _x(feature.end, seq_len)) / 2
        suffix = "#" if feature.attributes.get("fragment") else ""
        parts.append(_text(mid, 116, f"{feature.name}{suffix}", size=8, weight="bold"))

    for feature in cassettes:
        parts.append(_flat_feature_block(feature, seq_len, CASSETTE_COLOUR))
        mid = (_x(feature.start, seq_len) + _x(feature.end, seq_len)) / 2
        suffix = "#" if feature.element_type == "cassette_remnant" or feature.attributes.get("fragment") else ""
        parts.append(_text(mid, 116, f"{feature.name}{suffix}", size=8, weight="bold"))

    for feature in attc_sites:
        parts.append(_attc_marker(feature, seq_len))

    for feature in atti_sites:
        parts.append(_site_marker(feature, seq_len, INTEGRON_COLOUR))
        mid = (_x(feature.start, seq_len) + _x(feature.end, seq_len)) / 2
        parts.append(_text(mid, 103, feature.name, size=7, weight="bold"))

    for feature in promoters:
        parts.append(_promoter_flag(feature, seq_len))
        mid = (_x(feature.start, seq_len) + _x(feature.end, seq_len)) / 2
        parts.append(_text(mid, 116, feature.name, size=7, weight="bold"))

    for feature in insertion_sequences:
        parts.append(_block_arrow(feature, seq_len, "#ffffff"))
        mid = (_x(feature.start, seq_len) + _x(feature.end, seq_len)) / 2
        suffix = "#" if feature.attributes.get("fragment") or feature.attributes.get("type") == "p" else ""
        parts.append(_text(mid, 103, f"{feature.name}{suffix}", size=8, weight="bold"))

    for feature in iscrs:
        parts.append(_block_arrow(feature, seq_len, ISCR_COLOUR))
        mid = (_x(feature.start, seq_len) + _x(feature.end, seq_len)) / 2
        parts.append(_text(mid, 103, feature.name, size=8, weight="bold"))

    for feature in introns:
        parts.append(_block_arrow(feature, seq_len, INTRON_COLOUR))
        mid = (_x(feature.start, seq_len) + _x(feature.end, seq_len)) / 2
        parts.append(_text(mid, 103, feature.name, size=8, weight="bold"))

    for feature in captured_segments:
        parts.append(_flat_feature_block(feature, seq_len, CAPTURE_COLOUR))
        mid = (_x(feature.start, seq_len) + _x(feature.end, seq_len)) / 2
        parts.append(_text(mid, 116, feature.name, size=8))

    for feature in unknown_fragments:
        parts.append(_insertion_block(feature, seq_len))
        mid = (_x(feature.start, seq_len) + _x(feature.end, seq_len)) / 2
        parts.append(_text(mid, 116, f"{feature.name}#", size=8))

    for feature in replicons:
        parts.append(_block_arrow(feature, seq_len, REPLICON_COLOUR))
        mid = (_x(feature.start, seq_len) + _x(feature.end, seq_len)) / 2
        parts.append(_text(mid, 103, feature.name, size=8, weight="bold"))

    for feature in ncrnas:
        parts.append(_flat_feature_block(feature, seq_len, NCRNA_COLOUR))
        mid = (_x(feature.start, seq_len) + _x(feature.end, seq_len)) / 2
        parts.append(_text(mid, 116, feature.name, size=8, style="italic"))

    for feature in other_sites:
        colour = ISCR_COLOUR if feature.element_type in {"oriIS", "terIS"} else REPLICON_COLOUR
        parts.append(_site_marker(feature, seq_len, colour))
        mid = (_x(feature.start, seq_len) + _x(feature.end, seq_len)) / 2
        parts.append(_text(mid, 103, feature.name, size=7, weight="bold"))

    for feature in res_sites:
        parts.append(_res_marker(feature, seq_len))

    for feature in insertions:
        parts.append(_insertion_block(feature, seq_len))
        mid = (_x(feature.start, seq_len) + _x(feature.end, seq_len)) / 2
        inserted = int(feature.attributes.get("inserted_bases", 0) or 0)
        parts.append(_text(mid, 105, f"insertion ≈{inserted} bp", size=8))

    for feature in sorted(genes, key=lambda item: item.start):
        parts.append(_gene_arrow(feature, seq_len))
        mid = (_x(feature.start, seq_len) + _x(feature.end, seq_len)) / 2
        parts.append(_text(mid, 105, feature.name, size=9, weight="bold", style="italic"))

    for feature in res_sites:
        mid = (_x(feature.start, seq_len) + _x(feature.end, seq_len)) / 2
        parts.append(_text(mid, 105, "res", size=8, style="italic"))

    for feature in irs:
        parts.append(_ir_flag(feature, seq_len))

    for feature in direct_repeats:
        mid = (_x(feature.start, seq_len) + _x(feature.end, seq_len)) / 2
        sequence = str(feature.attributes.get("sequence") or feature.tsd_seq or "") or None
        parts.append(_lollipop(mid, sequence, feature.tsd_length))

    for feature in explicit_gaps:
        parts.extend(_gap_line(feature.start, feature.end, seq_len))

    visible = (
        major_units
        + integrons
        + conserved_segments
        + components
        + cassettes
        + cassette_arrays
        + attc_sites
        + atti_sites
        + promoters
        + insertion_sequences
        + iscrs
        + introns
        + genes
        + res_sites
        + insertions
        + captured_segments
        + unknown_fragments
        + replicons
        + ncrnas
        + other_sites
        + direct_repeats
        + context_regions
        + explicit_gaps
    )
    auto_gap_features = [
        feature for feature in visible if feature.element_type not in {"annotation_gap", "unknown_gap"}
    ]
    parts.extend(_long_unannotated_gaps(auto_gap_features, major_units, seq_len))
    parts.extend(_edge_gaps(auto_gap_features, seq_len))

    for feature in major_units + insertion_sequences + integrons:
        if feature.tsd_seq:
            parts.append(_lollipop(
                _x(feature.start, seq_len), feature.tsd_seq, feature.tsd_length,
            ))
            parts.append(_lollipop(
                _x(feature.end, seq_len), feature.tsd_seq, feature.tsd_length,
            ))
        elif feature.tsd_length:
            evidence = str(feature.attributes.get("tsd_evidence", ""))
            if not evidence:
                evidence = (
                    "searched_not_found"
                    if feature.start > feature.tsd_length
                    and feature.end + feature.tsd_length <= seq_len
                    else "untestable_missing_flank"
                )
            parts.append(_lollipop(
                _x(feature.start, seq_len), None, feature.tsd_length,
                evidence=evidence,
            ))
            parts.append(_lollipop(
                _x(feature.end, seq_len), None, feature.tsd_length,
                evidence=evidence,
            ))

    if sample_name:
        parts.append(_text(WIDTH - RIGHT, 151, sample_name, size=9, anchor="end", fill="#555"))

    parts.append(_legend())
    parts.append("</svg>")
    return "\n".join(parts)
