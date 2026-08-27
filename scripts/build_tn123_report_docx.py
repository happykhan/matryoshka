#!/usr/bin/env python3
"""Build the reviewable Tn1/Tn2/Tn3 report DOCX from its Markdown source."""

from __future__ import annotations

import re
from pathlib import Path

from docx import Document
from docx.enum.table import WD_CELL_VERTICAL_ALIGNMENT, WD_TABLE_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Inches, Pt, RGBColor

REPO = Path(__file__).parents[1]
SOURCE = REPO / "docs" / "tn123-definition-and-validation-report.md"
OUTPUT = REPO / "reports" / "Matryoshka-Tn123-definition-and-validation-report.docx"

BLUE = "1F4D78"
MID_BLUE = "2E74B5"
LIGHT_BLUE = "E8EEF5"
LIGHT_GRAY = "F2F4F7"
INK = "182026"
MUTED = "5B6570"
GREEN = "D7F5E5"


def _font(run, name: str, size: float, *, bold: bool = False, color: str = INK) -> None:
    run.font.name = name
    run._element.get_or_add_rPr().rFonts.set(qn("w:ascii"), name)
    run._element.get_or_add_rPr().rFonts.set(qn("w:hAnsi"), name)
    run.font.size = Pt(size)
    run.bold = bold
    run.font.color.rgb = RGBColor.from_string(color)


def _shade(element, fill: str) -> None:
    properties = element.get_or_add_tcPr() if hasattr(element, "get_or_add_tcPr") else element
    shading = properties.find(qn("w:shd"))
    if shading is None:
        shading = OxmlElement("w:shd")
        properties.append(shading)
    shading.set(qn("w:fill"), fill)


def _set_cell_margins(cell, top: int = 80, start: int = 120, bottom: int = 80, end: int = 120) -> None:
    tc_pr = cell._tc.get_or_add_tcPr()
    tc_mar = tc_pr.first_child_found_in("w:tcMar")
    if tc_mar is None:
        tc_mar = OxmlElement("w:tcMar")
        tc_pr.append(tc_mar)
    for edge, value in (("top", top), ("start", start), ("bottom", bottom), ("end", end)):
        node = tc_mar.find(qn(f"w:{edge}"))
        if node is None:
            node = OxmlElement(f"w:{edge}")
            tc_mar.append(node)
        node.set(qn("w:w"), str(value))
        node.set(qn("w:type"), "dxa")


def _set_cell_width(cell, width: int) -> None:
    tc_pr = cell._tc.get_or_add_tcPr()
    tc_w = tc_pr.find(qn("w:tcW"))
    if tc_w is None:
        tc_w = OxmlElement("w:tcW")
        tc_pr.append(tc_w)
    tc_w.set(qn("w:w"), str(width))
    tc_w.set(qn("w:type"), "dxa")


def _set_table_geometry(table, widths: list[int]) -> None:
    table.autofit = False
    table.alignment = WD_TABLE_ALIGNMENT.CENTER
    tbl_pr = table._tbl.tblPr
    layout = tbl_pr.find(qn("w:tblLayout"))
    if layout is None:
        layout = OxmlElement("w:tblLayout")
        tbl_pr.append(layout)
    layout.set(qn("w:type"), "fixed")
    tbl_w = tbl_pr.find(qn("w:tblW"))
    if tbl_w is None:
        tbl_w = OxmlElement("w:tblW")
        tbl_pr.append(tbl_w)
    tbl_w.set(qn("w:w"), "9360")
    tbl_w.set(qn("w:type"), "dxa")
    tbl_ind = tbl_pr.find(qn("w:tblInd"))
    if tbl_ind is None:
        tbl_ind = OxmlElement("w:tblInd")
        tbl_pr.append(tbl_ind)
    tbl_ind.set(qn("w:w"), "120")
    tbl_ind.set(qn("w:type"), "dxa")
    grid = table._tbl.tblGrid
    for child in list(grid):
        grid.remove(child)
    for width in widths:
        column = OxmlElement("w:gridCol")
        column.set(qn("w:w"), str(width))
        grid.append(column)
    for row in table.rows:
        for cell, width in zip(row.cells, widths, strict=False):
            _set_cell_width(cell, width)
            _set_cell_margins(cell)
            cell.vertical_alignment = WD_CELL_VERTICAL_ALIGNMENT.CENTER


def _table_widths(count: int) -> list[int]:
    return {
        2: [2800, 6560],
        4: [2100, 1800, 3900, 1560],
        5: [1750, 1400, 1800, 2500, 1910],
    }.get(count, [9360 // count] * count)


def _add_page_number(paragraph) -> None:
    paragraph.alignment = WD_ALIGN_PARAGRAPH.RIGHT
    run = paragraph.add_run("Page ")
    _font(run, "Calibri", 9, color=MUTED)
    field = OxmlElement("w:fldSimple")
    field.set(qn("w:instr"), "PAGE")
    paragraph._p.append(field)


def _add_inline(paragraph, text: str, *, base_size: float = 11) -> None:
    pattern = re.compile(r"(\*\*[^*]+\*\*|`[^`]+`)")
    for token in pattern.split(text):
        if not token:
            continue
        if token.startswith("**") and token.endswith("**"):
            run = paragraph.add_run(token[2:-2])
            _font(run, "Calibri", base_size, bold=True)
        elif token.startswith("`") and token.endswith("`"):
            run = paragraph.add_run(token[1:-1])
            _font(run, "Consolas", base_size - 1, color=BLUE)
        else:
            run = paragraph.add_run(token)
            _font(run, "Calibri", base_size)


def _set_num(paragraph, num_id: int) -> None:
    p_pr = paragraph._p.get_or_add_pPr()
    num_pr = OxmlElement("w:numPr")
    ilvl = OxmlElement("w:ilvl")
    ilvl.set(qn("w:val"), "0")
    num = OxmlElement("w:numId")
    num.set(qn("w:val"), str(num_id))
    num_pr.extend((ilvl, num))
    p_pr.append(num_pr)


def _numbering(document: Document, fmt: str, marker: str) -> int:
    root = document.part.numbering_part.element
    abstract_ids = [int(node.get(qn("w:abstractNumId"))) for node in root.findall(qn("w:abstractNum"))]
    num_ids = [int(node.get(qn("w:numId"))) for node in root.findall(qn("w:num"))]
    abstract_id = max(abstract_ids, default=0) + 1
    num_id = max(num_ids, default=0) + 1
    abstract = OxmlElement("w:abstractNum")
    abstract.set(qn("w:abstractNumId"), str(abstract_id))
    multi = OxmlElement("w:multiLevelType")
    multi.set(qn("w:val"), "singleLevel")
    level = OxmlElement("w:lvl")
    level.set(qn("w:ilvl"), "0")
    start = OxmlElement("w:start")
    start.set(qn("w:val"), "1")
    num_fmt = OxmlElement("w:numFmt")
    num_fmt.set(qn("w:val"), fmt)
    level_text = OxmlElement("w:lvlText")
    level_text.set(qn("w:val"), marker)
    justification = OxmlElement("w:lvlJc")
    justification.set(qn("w:val"), "left")
    p_pr = OxmlElement("w:pPr")
    tabs = OxmlElement("w:tabs")
    tab = OxmlElement("w:tab")
    tab.set(qn("w:val"), "num")
    tab.set(qn("w:pos"), "720")
    tabs.append(tab)
    indent = OxmlElement("w:ind")
    indent.set(qn("w:left"), "720")
    indent.set(qn("w:hanging"), "360")
    spacing = OxmlElement("w:spacing")
    spacing.set(qn("w:after"), "160")
    spacing.set(qn("w:line"), "280")
    spacing.set(qn("w:lineRule"), "auto")
    p_pr.extend((tabs, indent, spacing))
    level.extend((start, num_fmt, level_text, justification, p_pr))
    abstract.extend((multi, level))
    root.append(abstract)
    num = OxmlElement("w:num")
    num.set(qn("w:numId"), str(num_id))
    abstract_ref = OxmlElement("w:abstractNumId")
    abstract_ref.set(qn("w:val"), str(abstract_id))
    num.append(abstract_ref)
    root.append(num)
    return num_id


def _configure(document: Document) -> tuple[int, int]:
    section = document.sections[0]
    section.page_width = Inches(8.5)
    section.page_height = Inches(11)
    section.top_margin = Inches(1)
    section.right_margin = Inches(1)
    section.bottom_margin = Inches(1)
    section.left_margin = Inches(1)
    section.header_distance = Inches(0.492)
    section.footer_distance = Inches(0.492)

    styles = document.styles
    normal = styles["Normal"]
    normal.font.name = "Calibri"
    normal._element.rPr.rFonts.set(qn("w:ascii"), "Calibri")
    normal._element.rPr.rFonts.set(qn("w:hAnsi"), "Calibri")
    normal.font.size = Pt(11)
    normal.font.color.rgb = RGBColor.from_string(INK)
    normal.paragraph_format.space_before = Pt(0)
    normal.paragraph_format.space_after = Pt(6)
    normal.paragraph_format.line_spacing = 1.1
    for name, size, color, before, after in (
        ("Heading 1", 16, MID_BLUE, 16, 8),
        ("Heading 2", 13, MID_BLUE, 12, 6),
        ("Heading 3", 12, BLUE, 8, 4),
    ):
        style = styles[name]
        style.font.name = "Calibri"
        style._element.rPr.rFonts.set(qn("w:ascii"), "Calibri")
        style._element.rPr.rFonts.set(qn("w:hAnsi"), "Calibri")
        style.font.size = Pt(size)
        style.font.bold = True
        style.font.color.rgb = RGBColor.from_string(color)
        style.paragraph_format.space_before = Pt(before)
        style.paragraph_format.space_after = Pt(after)
        style.paragraph_format.keep_with_next = True

    header = section.header.paragraphs[0]
    header.text = "MATRYOSHKA  |  Tn1/Tn2/Tn3 definition and validation"
    header.alignment = WD_ALIGN_PARAGRAPH.LEFT
    _font(header.runs[0], "Calibri", 8.5, bold=True, color=MUTED)
    p_pr = header._p.get_or_add_pPr()
    border = OxmlElement("w:pBdr")
    bottom = OxmlElement("w:bottom")
    bottom.set(qn("w:val"), "single")
    bottom.set(qn("w:sz"), "4")
    bottom.set(qn("w:space"), "4")
    bottom.set(qn("w:color"), "CCD3DA")
    border.append(bottom)
    p_pr.append(border)
    _add_page_number(section.footer.paragraphs[0])
    return _numbering(document, "bullet", "•"), _numbering(document, "decimal", "%1.")


def _add_table(document: Document, source_rows: list[list[str]]) -> None:
    if len(source_rows) < 2:
        return
    rows = [source_rows[0]] + source_rows[2:]
    table = document.add_table(rows=len(rows), cols=len(rows[0]))
    table.style = "Table Grid"
    _set_table_geometry(table, _table_widths(len(rows[0])))
    for row_index, values in enumerate(rows):
        for column_index, value in enumerate(values):
            cell = table.cell(row_index, column_index)
            paragraph = cell.paragraphs[0]
            paragraph.paragraph_format.space_after = Pt(2)
            _add_inline(paragraph, value.strip(), base_size=9.2)
            if row_index == 0:
                _shade(cell._tc, LIGHT_GRAY)
                for run in paragraph.runs:
                    run.bold = True
            elif value.strip() == "PASS":
                _shade(cell._tc, GREEN)
                for run in paragraph.runs:
                    run.bold = True
    document.add_paragraph().paragraph_format.space_after = Pt(0)


def build() -> None:
    document = Document()
    bullet_id, number_id = _configure(document)
    lines = SOURCE.read_text(encoding="utf-8").splitlines()
    index = 0
    in_code = False
    code_lines: list[str] = []
    first_title = True
    active_list: str | None = None
    while index < len(lines):
        line = lines[index]
        if line.startswith("```"):
            if in_code:
                paragraph = document.add_paragraph()
                paragraph.paragraph_format.left_indent = Inches(0.18)
                paragraph.paragraph_format.right_indent = Inches(0.18)
                paragraph.paragraph_format.space_before = Pt(4)
                paragraph.paragraph_format.space_after = Pt(8)
                run = paragraph.add_run("\n".join(code_lines))
                _font(run, "Consolas", 8.5, color=INK)
                shading = OxmlElement("w:shd")
                shading.set(qn("w:fill"), LIGHT_BLUE)
                paragraph._p.get_or_add_pPr().append(shading)
                code_lines = []
                in_code = False
            else:
                in_code = True
            index += 1
            continue
        if in_code:
            code_lines.append(line)
            index += 1
            continue
        if not line.strip():
            active_list = None
            index += 1
            continue
        if line.startswith("| "):
            active_list = None
            table_rows: list[list[str]] = []
            while index < len(lines) and lines[index].startswith("|"):
                table_rows.append([value.strip() for value in lines[index].strip().strip("|").split("|")])
                index += 1
            _add_table(document, table_rows)
            continue
        if line.startswith("# ") and first_title:
            active_list = None
            paragraph = document.add_paragraph()
            paragraph.paragraph_format.space_before = Pt(0)
            paragraph.paragraph_format.space_after = Pt(12)
            paragraph.paragraph_format.keep_with_next = True
            run = paragraph.add_run(line[2:])
            _font(run, "Calibri", 24, bold=True, color=BLUE)
            first_title = False
        elif line.startswith("## "):
            active_list = None
            if line[3:] == "Review questions for Sally":
                document.add_page_break()
            document.add_paragraph(line[3:], style="Heading 1")
        elif line.startswith("### "):
            active_list = None
            document.add_paragraph(line[4:], style="Heading 2")
        elif line.startswith("> "):
            active_list = None
            paragraph = document.add_paragraph()
            paragraph.paragraph_format.left_indent = Inches(0.3)
            paragraph.paragraph_format.right_indent = Inches(0.2)
            paragraph.paragraph_format.space_before = Pt(6)
            paragraph.paragraph_format.space_after = Pt(8)
            _add_inline(paragraph, line[2:])
            shading = OxmlElement("w:shd")
            shading.set(qn("w:fill"), LIGHT_BLUE)
            paragraph._p.get_or_add_pPr().append(shading)
        elif line.startswith("- "):
            if active_list != "bullet":
                bullet_id = _numbering(document, "bullet", "•")
            active_list = "bullet"
            paragraph = document.add_paragraph()
            _set_num(paragraph, bullet_id)
            _add_inline(paragraph, line[2:])
        elif re.match(r"^\d+\. ", line):
            if active_list != "number":
                number_id = _numbering(document, "decimal", "%1.")
            active_list = "number"
            paragraph = document.add_paragraph()
            _set_num(paragraph, number_id)
            _add_inline(paragraph, re.sub(r"^\d+\. ", "", line))
        else:
            active_list = None
            paragraph_lines = [line]
            while index + 1 < len(lines):
                next_line = lines[index + 1]
                if (
                    not next_line.strip()
                    or next_line.startswith(("#", "- ", "> ", "```", "| "))
                    or re.match(r"^\d+\. ", next_line)
                ):
                    break
                paragraph_lines.append(next_line)
                index += 1
            paragraph = document.add_paragraph()
            _add_inline(paragraph, " ".join(part.strip() for part in paragraph_lines))
        index += 1

    document.core_properties.title = "Automated Tn1, Tn2 and Tn3 definition and validation report"
    document.core_properties.subject = "Matryoshka expert-rule architecture and real-accession validation"
    document.core_properties.author = "Matryoshka project"
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    document.save(OUTPUT)
    print(OUTPUT)


if __name__ == "__main__":
    build()
