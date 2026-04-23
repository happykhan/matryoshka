"""
viz.py — Scale-accurate linear SVG map of the MGE hierarchy.

Each element is drawn at its exact genomic coordinate. Nesting levels are shown
as vertical tracks: compound elements (transposons, integrons) occupy a taller
outer row; their children are drawn as arrows inside them. Singletons sit on a
thin baseline track.

Element style:
  - IS elements: directional pentagon arrows (strand + / -)
  - AMR genes: filled rectangles, bold label
  - Transposons/integrons: open rounded box with children inside
  - attC sites: thin tick mark
  - Cassettes: plain rectangle

Colour coding matches the wolvercote palette.
"""

from __future__ import annotations
import math
from .detect import MGEFeature

# ── Colours ─────────────────────────────────────────────────────────────────
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

# ── Layout ───────────────────────────────────────────────────────────────────
W        = 1400   # total canvas width
ML       = 80     # left margin (room for labels)
MR       = 40     # right margin
DW       = W - ML - MR

COMPOUND_OUTER_H = 90    # height of transposon/integron outer box
COMPOUND_INNER_H = 36    # height of child arrows inside compound box
SINGLETON_H      = 36    # height of standalone elements
TRACK_GAP        = 18    # vertical gap between compound rows / singleton row
SCALEBAR_H       = 28    # height of scale bar area
TITLE_H          = 34

FONT = "monospace"


def _c(et: str) -> str:
    return COLORS.get(et, _DEFAULT)


def _x(pos: int, seq_len: int) -> float:
    return ML + (pos / seq_len) * DW


def _px(start: int, end: int, seq_len: int) -> tuple[float, float]:
    return _x(start, seq_len), _x(end, seq_len)


def _esc(s: str) -> str:
    return s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def _trunc(s: str, px_width: float, char_px: float = 7.0) -> str:
    max_c = max(1, int(px_width / char_px))
    return s if len(s) <= max_c else s[:max_c - 1] + "…"


# ── Primitive helpers ────────────────────────────────────────────────────────

def _rect(x: float, y: float, w: float, h: float, fill: str,
          stroke: str = "none", sw: float = 1.5, rx: float = 3,
          opacity: float = 1.0) -> str:
    return (
        f'<rect x="{x:.2f}" y="{y:.2f}" width="{max(w,1):.2f}" height="{h:.2f}" '
        f'fill="{fill}" stroke="{stroke}" stroke-width="{sw}" rx="{rx}" opacity="{opacity}"/>'
    )


def _text(x: float, y: float, s: str, size: int = 11,
          anchor: str = "middle", fill: str = "white", bold: bool = False,
          italic: bool = False) -> str:
    style = ""
    if bold:
        style += "font-weight:bold;"
    if italic:
        style += "font-style:italic;"
    style_attr = f' style="{style}"' if style else ""
    return (
        f'<text x="{x:.2f}" y="{y:.2f}" font-size="{size}" font-family="{FONT}" '
        f'text-anchor="{anchor}" fill="{fill}"{style_attr}>{_esc(s)}</text>'
    )


def _arrow(x: float, y: float, w: float, h: float,
           strand: str, fill: str, stroke: str = "white") -> str:
    """Pentagon arrow for IS elements / directional genes."""
    tip = min(h * 0.55, max(w * 0.22, 6.0))
    if w < 3:
        return _rect(x, y, w, h, fill)

    if strand == "+":
        pts = (f"{x:.2f},{y:.2f} "
               f"{x+w-tip:.2f},{y:.2f} "
               f"{x+w:.2f},{y+h/2:.2f} "
               f"{x+w-tip:.2f},{y+h:.2f} "
               f"{x:.2f},{y+h:.2f}")
    elif strand == "-":
        pts = (f"{x+tip:.2f},{y:.2f} "
               f"{x+w:.2f},{y:.2f} "
               f"{x+w:.2f},{y+h:.2f} "
               f"{x+tip:.2f},{y+h:.2f} "
               f"{x:.2f},{y+h/2:.2f}")
    else:
        return _rect(x, y, w, h, fill, stroke=stroke)

    return f'<polygon points="{pts}" fill="{fill}" stroke="{stroke}" stroke-width="0.8"/>'


def _tick(x1: float, x2: float, y: float, h: float, fill: str) -> str:
    mid = (x1 + x2) / 2
    return (
        f'<line x1="{mid:.2f}" y1="{y:.2f}" x2="{mid:.2f}" y2="{y+h:.2f}" '
        f'stroke="{fill}" stroke-width="2"/>'
    )


# ── Element renderers ────────────────────────────────────────────────────────

def _draw_element(f: MGEFeature, y: float, h: float, seq_len: int,
                  label_outside: bool = False) -> list[str]:
    """Render a single leaf element (IS arrow, AMR rect, attC tick, cassette rect)."""
    x1, x2 = _px(f.start, f.end, seq_len)
    w = max(x2 - x1, 1.0)
    fill = _c(f.element_type)
    mid_x = (x1 + x2) / 2
    label_y = y + h / 2 + 4
    parts: list[str] = [f'<g><title>{f.element_type}: {_esc(f.name)} ({f.start}–{f.end})</title>']

    if f.element_type == "attC":
        parts.append(_tick(x1, x2, y, h, fill))
    elif f.element_type in ("IS",):
        parts.append(_arrow(x1, y, w, h, f.strand, fill))
        if w >= 20:
            label = _trunc(f.name, w - 6)
            parts.append(_text(mid_x, label_y, label, size=9))
    else:
        # AMR, cassette, integrase, replicon → filled rect
        sw = 1.5 if f.element_type == "AMR" else 0.5
        parts.append(_rect(x1, y, w, h, fill, stroke="white", sw=sw))
        if w >= 24:
            label = _trunc(f.name, w - 6)
            parts.append(_text(mid_x, label_y, label,
                               size=9, bold=(f.element_type == "AMR")))

    if label_outside and w < 24:
        parts.append(_text(x1 + w / 2, y - 3, _trunc(f.name, 12),
                           size=8, fill="#444", anchor="middle"))

    parts.append("</g>")
    return parts


def _draw_compound(f: MGEFeature, y: float, seq_len: int) -> list[str]:
    """Render a compound element (transposon / integron) with children inside."""
    x1, x2 = _px(f.start, f.end, seq_len)
    w = max(x2 - x1, 2.0)
    fill = _c(f.element_type)
    parts: list[str] = [f'<g><title>{f.element_type}: {_esc(f.name)} ({f.start}–{f.end})</title>']

    # Outer rounded box — semi-transparent fill, solid stroke
    parts.append(_rect(x1, y, w, COMPOUND_OUTER_H, fill,
                       stroke=fill, sw=2.0, rx=5, opacity=0.12))
    parts.append(_rect(x1, y, w, COMPOUND_OUTER_H, "none",
                       stroke=fill, sw=2.0, rx=5))

    # Top label (name + span)
    span_kb = (f.end - f.start) / 1000
    label = f"{f.name}  ({span_kb:.1f} kb)"
    parts.append(_text(x1 + 6, y + 13, label, size=10,
                       fill=fill, anchor="start", bold=True))

    # Children — centred vertically inside the outer box
    child_y = y + COMPOUND_OUTER_H - COMPOUND_INNER_H - 8
    for child in sorted(f.children, key=lambda c: c.start):
        parts.extend(_draw_element(child, child_y, COMPOUND_INNER_H, seq_len))

    parts.append("</g>")
    return parts


# ── Scale bar ────────────────────────────────────────────────────────────────

def _scale_bar(y: float, seq_len: int) -> list[str]:
    parts = [
        f'<line x1="{ML}" y1="{y+8}" x2="{ML+DW}" y2="{y+8}" '
        f'stroke="#aaa" stroke-width="1.5"/>',
    ]
    n = 5
    for i in range(n + 1):
        pos = int(i * seq_len / n)
        sx = _x(pos, seq_len)
        parts.append(
            f'<line x1="{sx:.2f}" y1="{y+3}" x2="{sx:.2f}" y2="{y+13}" '
            f'stroke="#aaa" stroke-width="1.5"/>'
        )
        label = f"{pos // 1000}k" if pos else "0"
        parts.append(_text(sx, y + 24, label, size=9, fill="#666"))
    return parts


# ── Legend ───────────────────────────────────────────────────────────────────

def _legend(y: float) -> list[str]:
    items = [
        ("IS element",  "IS"),
        ("AMR gene",    "AMR"),
        ("Transposon",  "transposon"),
        ("Integron",    "integron"),
        ("Cassette",    "cassette"),
        ("attC site",   "attC"),
    ]
    parts: list[str] = []
    lx = ML
    for label, et in items:
        c = _c(et)
        parts.append(_rect(lx, y, 14, 14, c, rx=2))
        parts.append(_text(lx + 18, y + 11, label, size=10,
                           fill="#444", anchor="start"))
        lx += 130
    return parts


# ── Main entry point ─────────────────────────────────────────────────────────

def to_linear_svg(roots: list[MGEFeature], seq_len: int,
                  sample_name: str = "") -> str:
    """
    Render the MGE hierarchy as a scale-accurate linear SVG map.

    roots    — output of build_hierarchy()
    seq_len  — total sequence length in bp
    """
    compounds  = sorted([r for r in roots if r.children], key=lambda r: r.start)
    singletons = sorted([r for r in roots if not r.children], key=lambda r: r.start)

    # ── Vertical layout ──────────────────────────────────────────────────────
    cursor = TITLE_H + SCALEBAR_H

    # Backbone at top of element section
    backbone_y = cursor
    cursor += 10

    compound_rows: list[tuple[float, MGEFeature]] = []
    for f in compounds:
        compound_rows.append((cursor, f))
        cursor += COMPOUND_OUTER_H + TRACK_GAP

    singleton_y = cursor
    cursor += SINGLETON_H + TRACK_GAP
    legend_y = cursor + 6
    total_h = legend_y + 30

    parts: list[str] = [
        f'<svg xmlns="http://www.w3.org/2000/svg" '
        f'width="{W}" height="{total_h}" '
        f'style="background:#fff;font-family:{FONT};">',

        # Title
        _text(W / 2, TITLE_H - 6, sample_name or "MGE linear map",
              size=15, fill="#222", bold=True),
    ]

    # Scale bar
    parts.extend(_scale_bar(TITLE_H, seq_len))

    # Backbone line
    parts.append(
        f'<line x1="{ML}" y1="{backbone_y}" x2="{ML+DW}" y2="{backbone_y}" '
        f'stroke="#ccc" stroke-width="3"/>'
    )

    # Connector lines from backbone down to compound elements
    for (y, f) in compound_rows:
        x1, x2 = _px(f.start, f.end, seq_len)
        c = _c(f.element_type)
        for sx in (x1, x2):
            parts.append(
                f'<line x1="{sx:.2f}" y1="{backbone_y}" '
                f'x2="{sx:.2f}" y2="{y:.2f}" '
                f'stroke="{c}" stroke-width="1" stroke-dasharray="3,3" opacity="0.5"/>'
            )

    # Compound element rows
    for (y, f) in compound_rows:
        parts.extend(_draw_compound(f, y, seq_len))

    # Singletons on backbone track
    for f in singletons:
        parts.extend(_draw_element(f, singleton_y, SINGLETON_H, seq_len,
                                   label_outside=True))

    # Backbone element markers (small ticks above backbone for singletons)
    for f in singletons:
        x1, x2 = _px(f.start, f.end, seq_len)
        mid = (x1 + x2) / 2
        c = _c(f.element_type)
        parts.append(
            f'<line x1="{mid:.2f}" y1="{backbone_y - 4}" '
            f'x2="{mid:.2f}" y2="{backbone_y + 4}" '
            f'stroke="{c}" stroke-width="2"/>'
        )

    # Legend
    parts.extend(_legend(legend_y))

    parts.append("</svg>")
    return "\n".join(parts)
