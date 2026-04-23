"""
detect.py — Run external detection tools and parse outputs into a unified
coordinate table.

Supported tools:
  - ISEScan (IS elements)
  - IntegronFinder (integrons + gene cassettes)
  - AMRFinder+ (resistance genes)
  - MOBsuite (plasmid replicons + mobility)

Each function returns a list of dicts with keys:
  element_type, family, name, start, end, strand, score, attributes
"""

from dataclasses import dataclass, field


@dataclass
class MGEFeature:
    element_type: str          # IS, transposon, integron, cassette, AMR, replicon
    family: str                # e.g. IS6, Tn3, class1_integron
    name: str                  # e.g. IS26, IntI1, blaCTX-M-15
    start: int
    end: int
    strand: str                # + or -
    tsd_length: int | None = None   # expected TSD length for this family
    tsd_seq: str | None = None      # confirmed TSD sequence if found
    ir_left: str | None = None
    ir_right: str | None = None
    score: float | None = None
    attributes: dict = field(default_factory=dict)
    children: list = field(default_factory=list)   # nested elements


# TSD lengths by IS family
TSD_LENGTHS = {
    "IS6":    8,   # IS26
    "IS1380": 5,   # ISEcp1
    "Tn3":    5,
    "Tn7":    5,
    "IS30":   2,   # too short, skip confirmation
    "IS91":   None,  # rolling circle, no TSD
    "ISCR":   None,  # rolling circle, no TSD
}
