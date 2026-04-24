"""
viz.py — Scale-accurate linear SVG map of the MGE hierarchy.

Elements are drawn at their exact genomic coordinates (bp-scaled). Compound
elements (transposons, integrons) are nested visually: each compound occupies
a rounded box with its children drawn inside. Nesting is recursive so pEK499's
IS26-within-IS26 structure renders correctly.

TSD annotation: when a TSD sequence is confirmed on a compound element, two
small gold rectangles flank its boundaries with the TSD sequence labelled.

Colour coding:
  IS          #4a90d9  blue arrow
  AMR         #e74c3c  red rect
  transposon  #27ae60  green box
  integron    #8e44ad  purple box
  cassette    #f39c12  orange rect
  attC        #7f8c8d  grey tick
  integrase   #7f8c8d  grey rect
"""

from __future__ import annotations

from .detect import MGEFeature

# ── Colours ──────────────────────────────────────────────────────────────────
COLORS: dict[str, str] = {
    "IS":         "#4a90d9",
    "AMR":        "#e74c3c",
    "transposon": "#27ae60",
    "integron":   "#8e44ad",
    "cassette":   "#f39c12",
    "attC":       "#7f8c8d",
    "integrase":  "#7f8c8d",
    "replicon":   "#16a085",
}
_DEFAULT = "#95a5a6"
TSD_COLOR = "#f1c40f"   # gold

# ── Canvas ───────────────────────────────────────────────────────────────────
# Width is recomputed per-render based on seq_len so chromosomes get enough
# px-per-bp to remain legible. Module-level names stay so helper functions
# can reference them without threading a config through every call.
_WIDTH_MIN  = 1400           # px for anything up to ~150 kb
_WIDTH_MAX  = 6000           # cap so huge chromosomes don't blow up SVG size
_PX_PER_KB  = 8              # target resolution once past _WIDTH_MIN

W   = _WIDTH_MIN
ML  = 60
MR  = 40
DW  = W - ML - MR
FONT = "monospace"


def _set_canvas_width(seq_len: int, override: int | None = None) -> None:
    """Recompute module-level W / DW based on sequence length."""
    global W, DW
    if override is not None:
        W = max(400, int(override))
    else:
        scaled = int(seq_len / 1000 * _PX_PER_KB) + ML + MR
        W = max(_WIDTH_MIN, min(_WIDTH_MAX, scaled))
    DW = W - ML - MR

# ── Layout constants ─────────────────────────────────────────────────────────
LEAF_H      = 34    # height of a leaf element (IS, AMR, etc.)
LABEL_TOP   = 20    # px reserved for compound label above child area
BOT_PAD     = 8     # padding below child area inside compound box
INNER_PAD   = 6     # horizontal pad inside compound box
TRACK_GAP   = 14    # gap between compound rows and singleton row
TITLE_H     = 32
SCALEBAR_H  = 26
BACKBONE_BELOW_SCALE = 12   # gap from scale bar bottom to backbone line


def _c(et: str) -> str:
    return COLORS.get(et, _DEFAULT)


def _x(pos: int, seq_len: int) -> float:
    return ML + (pos / seq_len) * DW


def _px(s: int, e: int, seq_len: int) -> tuple[float, float]:
    return _x(s, seq_len), _x(e, seq_len)


def _esc(s: str) -> str:
    return s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def _trunc(s: str, px_w: float, ch: float = 7.0) -> str:
    n = max(1, int(px_w / ch))
    return s if len(s) <= n else s[:n - 1] + "…"


# ── Height calculation (recursive) ───────────────────────────────────────────

def _elem_h(f: MGEFeature) -> float:
    """Height needed to render f and all its descendants."""
    if not f.children:
        return LEAF_H
    max_child_h = max(_elem_h(c) for c in f.children)
    return LABEL_TOP + max_child_h + BOT_PAD


# ── SVG primitives ───────────────────────────────────────────────────────────

def _rect(x: float, y: float, w: float, h: float, fill: str,
          stroke: str = "none", sw: float = 1.5, rx: float = 3,
          opacity: float = 1.0) -> str:
    return (
        f'<rect x="{x:.2f}" y="{y:.2f}" width="{max(w,1):.2f}" height="{h:.2f}" '
        f'fill="{fill}" stroke="{stroke}" stroke-width="{sw}" rx="{rx}" '
        f'opacity="{opacity}"/>'
    )


def _text(x: float, y: float, s: str, size: int = 11,
          anchor: str = "middle", fill: str = "white",
          bold: bool = False) -> str:
    weight = ' font-weight="bold"' if bold else ""
    return (
        f'<text x="{x:.2f}" y="{y:.2f}" font-size="{size}" font-family="{FONT}" '
        f'text-anchor="{anchor}" fill="{fill}"{weight}>{_esc(s)}</text>'
    )


def _arrow(x: float, y: float, w: float, h: float,
           strand: str, fill: str) -> str:
    tip = min(h * 0.55, max(w * 0.22, 5.0))
    if w < 3:
        return _rect(x, y, w, h, fill)
    if strand == "+":
        pts = (f"{x:.2f},{y:.2f} {x+w-tip:.2f},{y:.2f} "
               f"{x+w:.2f},{y+h/2:.2f} {x+w-tip:.2f},{y+h:.2f} {x:.2f},{y+h:.2f}")
    elif strand == "-":
        pts = (f"{x+tip:.2f},{y:.2f} {x+w:.2f},{y:.2f} "
               f"{x+w:.2f},{y+h:.2f} {x+tip:.2f},{y+h:.2f} {x:.2f},{y+h/2:.2f}")
    else:
        return _rect(x, y, w, h, fill)
    return f'<polygon points="{pts}" fill="{fill}" stroke="white" stroke-width="0.8"/>'


# ── TSD annotation ───────────────────────────────────────────────────────────

def _draw_tsd(x_left: float, x_right: float, y: float, h: float,
              tsd: str) -> list[str]:
    """Draw gold TSD rectangles flanking element at x_left / x_right."""
    tw = max(len(tsd) * 5.5, 6.0)
    th = h * 0.55
    ty = y + (h - th) / 2
    label_y = y + h + 10

    parts = [
        # left TSD box
        _rect(x_left - tw - 1, ty, tw, th, TSD_COLOR, rx=1),
        _text(x_left - tw / 2 - 1, label_y, tsd, size=8, fill=TSD_COLOR),
        # right TSD box
        _rect(x_right + 1, ty, tw, th, TSD_COLOR, rx=1),
        _text(x_right + tw / 2 + 1, label_y, tsd, size=8, fill=TSD_COLOR),
    ]
    return parts


# ── Recursive element renderer ────────────────────────────────────────────────

def _draw_any(f: MGEFeature, y: float, h: float, seq_len: int,
              depth: int = 0) -> list[str]:
    """
    Render f at vertical position y with total height h.
    Recursively renders compound elements and their children.
    depth controls label size / padding.
    """
    x1, x2 = _px(f.start, f.end, seq_len)
    w = max(x2 - x1, 1.0)
    fill = _c(f.element_type)
    mid_x = (x1 + x2) / 2

    tooltip_extra = ""
    if f.tsd_seq:
        tooltip_extra += f" TSD={f.tsd_seq}"
    if f.ir_left:
        tooltip_extra += f" IR={f.ir_left[:10]}"
    parts: list[str] = [
        f'<g><title>{f.element_type}: {_esc(f.name)} ({f.start}–{f.end}){tooltip_extra}</title>'
    ]

    if not f.children:
        # ── Leaf element ──────────────────────────────────────────────────────
        if f.element_type == "attC":
            mid = (x1 + x2) / 2
            parts.append(
                f'<line x1="{mid:.2f}" y1="{y:.2f}" x2="{mid:.2f}" y2="{y+h:.2f}" '
                f'stroke="{fill}" stroke-width="2"/>'
            )
        elif f.element_type == "IS":
            parts.append(_arrow(x1, y, w, h, f.strand, fill))
            if w >= 18:
                parts.append(_text(mid_x, y + h / 2 + 4,
                                   _trunc(f.name, w - 4), size=8))
            # IR tick marks at both ends when inverted repeats confirmed
            if f.ir_left:
                ir_h = h * 0.7
                ir_y = y + (h - ir_h) / 2
                parts.append(f'<line x1="{x1+2:.2f}" y1="{ir_y:.2f}" '
                              f'x2="{x1+2:.2f}" y2="{ir_y+ir_h:.2f}" '
                              f'stroke="white" stroke-width="2" opacity="0.9"/>')
                parts.append(f'<line x1="{x2-2:.2f}" y1="{ir_y:.2f}" '
                              f'x2="{x2-2:.2f}" y2="{ir_y+ir_h:.2f}" '
                              f'stroke="white" stroke-width="2" opacity="0.9"/>')
        else:
            parts.append(_rect(x1, y, w, h, fill, stroke="white", sw=1.0))
            if w >= 22:
                parts.append(_text(mid_x, y + h / 2 + 4,
                                   _trunc(f.name, w - 4), size=8,
                                   bold=(f.element_type == "AMR")))
    else:
        # ── Compound element ──────────────────────────────────────────────────
        label_size = max(8, 11 - depth)
        # Outer box
        parts.append(_rect(x1, y, w, h, fill, stroke=fill,
                           sw=max(1.5, 2.5 - depth * 0.5),
                           rx=4, opacity=0.10))
        parts.append(_rect(x1, y, w, h, "none", stroke=fill,
                           sw=max(1.5, 2.5 - depth * 0.5), rx=4))

        # Name label top-left
        span_kb = (f.end - f.start) / 1000
        name_label = f"{f.name}  ({span_kb:.1f} kb)"
        parts.append(_text(x1 + 4, y + label_size + 2, name_label,
                           size=label_size, fill=fill, anchor="start", bold=True))

        # TSD flanking marks (below the box)
        if f.tsd_seq:
            parts.extend(_draw_tsd(x1, x2, y, h, f.tsd_seq))

        # Children — all share the same row inside the box
        child_y = y + LABEL_TOP
        child_h = h - LABEL_TOP - BOT_PAD
        for child in sorted(f.children, key=lambda c: c.start):
            parts.extend(_draw_any(child, child_y, child_h, seq_len, depth + 1))

    parts.append("</g>")
    return parts


# ── Scale bar ─────────────────────────────────────────────────────────────────

def _scale_bar(y: float, seq_len: int) -> list[str]:
    parts = [f'<line x1="{ML}" y1="{y+8}" x2="{ML+DW}" y2="{y+8}" '
             f'stroke="#bbb" stroke-width="1.5"/>']
    for i in range(6):
        pos = int(i * seq_len / 5)
        sx = _x(pos, seq_len)
        parts.append(f'<line x1="{sx:.2f}" y1="{y+4}" x2="{sx:.2f}" y2="{y+12}" '
                     f'stroke="#bbb" stroke-width="1.5"/>')
        label = f"{pos // 1000}k" if pos else "0"
        parts.append(_text(sx, y + 22, label, size=9, fill="#888"))
    return parts


# ── Legend ────────────────────────────────────────────────────────────────────

def _legend(y: float) -> list[str]:
    items = [
        ("IS element",  "IS"),
        ("AMR gene",    "AMR"),
        ("Transposon",  "transposon"),
        ("Integron",    "integron"),
        ("Cassette",    "cassette"),
        ("attC site",   "attC"),
        ("TSD",         "__tsd__"),
    ]
    parts: list[str] = []
    lx = ML
    for label, et in items:
        c = TSD_COLOR if et == "__tsd__" else _c(et)
        parts.append(_rect(lx, y, 14, 14, c, rx=2))
        parts.append(_text(lx + 18, y + 11, label, size=10,
                           fill="#444", anchor="start"))
        lx += 120
    return parts


# ── Main entry point ──────────────────────────────────────────────────────────

def to_linear_svg(roots: list[MGEFeature], seq_len: int,
                  sample_name: str = "", width: int | None = None) -> str:
    """Scale-accurate linear SVG of the MGE hierarchy.

    Canvas width auto-scales with sequence length unless `width` is given
    explicitly. A 60 kb plasmid renders at _WIDTH_MIN; a 5 Mb chromosome
    expands up to _WIDTH_MAX so dense feature clusters remain legible.
    """
    _set_canvas_width(seq_len, width)
    compounds  = sorted([r for r in roots if r.children], key=lambda r: r.start)
    singletons = sorted([r for r in roots if not r.children and r.start > 0],
                        key=lambda r: r.start)

    # ── Vertical layout ───────────────────────────────────────────────────────
    backbone_y = TITLE_H + SCALEBAR_H + BACKBONE_BELOW_SCALE
    cursor = backbone_y + 10

    compound_rows: list[tuple[float, float, MGEFeature]] = []
    for f in compounds:
        h = _elem_h(f)
        compound_rows.append((cursor, h, f))
        tsd_extra = 16 if f.tsd_seq else 0
        cursor += h + tsd_extra + TRACK_GAP

    singleton_y = cursor
    legend_y = singleton_y + LEAF_H + TRACK_GAP + 8
    total_h = legend_y + 30

    parts: list[str] = [
        f'<svg xmlns="http://www.w3.org/2000/svg" '
        f'width="{W}" height="{total_h:.0f}" '
        f'style="background:#fff;font-family:{FONT};">',
        _text(W / 2, TITLE_H - 4, sample_name or "MGE linear map",
              size=15, fill="#222", bold=True),
    ]

    parts.extend(_scale_bar(TITLE_H, seq_len))

    # Backbone line
    parts.append(
        f'<line x1="{ML}" y1="{backbone_y}" x2="{ML+DW}" y2="{backbone_y}" '
        f'stroke="#ccc" stroke-width="3"/>'
    )

    # Dashed connector lines from backbone to compound elements
    for (y, _h, f) in compound_rows:
        x1, x2 = _px(f.start, f.end, seq_len)
        c = _c(f.element_type)
        for sx in (x1, x2):
            parts.append(
                f'<line x1="{sx:.2f}" y1="{backbone_y}" x2="{sx:.2f}" y2="{y:.2f}" '
                f'stroke="{c}" stroke-width="1" stroke-dasharray="3,3" opacity="0.4"/>'
            )

    # Draw compound rows (recursive)
    for (y, h, f) in compound_rows:
        parts.extend(_draw_any(f, y, h, seq_len))

    # Singletons on baseline
    for f in singletons:
        parts.extend(_draw_any(f, singleton_y, LEAF_H, seq_len))
        # Small tick above backbone
        x1, x2 = _px(f.start, f.end, seq_len)
        mid = (x1 + x2) / 2
        c = _c(f.element_type)
        parts.append(
            f'<line x1="{mid:.2f}" y1="{backbone_y - 4}" '
            f'x2="{mid:.2f}" y2="{backbone_y + 4}" '
            f'stroke="{c}" stroke-width="2"/>'
        )

    parts.extend(_legend(legend_y))
    parts.append("</svg>")
    return "\n".join(parts)
