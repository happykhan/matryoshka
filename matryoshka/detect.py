"""
detect.py — Run external detection tools and parse outputs into a unified
coordinate table.

Supported tools:
  - ISEScan (IS elements)
  - IntegronFinder (integrons + gene cassettes)
  - AMRFinder+ (resistance genes)
  - MOBsuite (plasmid replicons + mobility)

Each parser returns list[MGEFeature]. Coordinates are 1-based, inclusive (GFF3 convention).
"""

from __future__ import annotations

import csv
import io
import subprocess
from dataclasses import dataclass, field
from pathlib import Path


@dataclass
class MGEFeature:
    element_type: str          # IS, integron, cassette, attC, AMR, replicon, mob
    family: str                # e.g. IS6, Tn3, class1_integron
    name: str                  # e.g. IS26, IntI1, blaCTX-M-15
    start: int                 # 1-based
    end: int                   # 1-based, inclusive
    strand: str                # + or -
    tsd_length: int | None = None   # expected TSD length for this family
    tsd_seq: str | None = None      # confirmed TSD sequence if found
    ir_left: str | None = None
    ir_right: str | None = None
    score: float | None = None
    attributes: dict = field(default_factory=dict)
    children: list = field(default_factory=list)


# Expected TSD lengths by IS superfamily or transposon family.
# Values per Partridge et al. 2018 (PMC6148190) Tables 1-2.
TSD_LENGTHS = {
    # IS families
    "IS6":    8,      # IS26, IS257, IS1216
    "IS1380": 5,      # ISEcp1 (5-6bp, 6 occasional), ISKpn23
    "IS30":   2,      # ISApl1, ISAba125
    "IS91":   None,   # rolling circle — no TSD
    "ISCR":   None,   # rolling circle — no TSD (generic prefix)
    "ISCR1":  None,
    "ISCR2":  None,
    # Transposon families (most Tn3 superfamily = 5bp TSD)
    "Tn3":     5,
    "Tn2":     5,
    "Tn21":    5,
    "Tn1696":  5,
    "Tn1721":  5,
    "Tn1546":  5,
    "Tn1331":  5,
    "Tn5393":  5,
    "Tn6452":  5,
    "Tn7":     5,
    "Tn4401":  5,     # ISKpn7-bounded Tn3 variant
    "Tn1999":  5,
    "Tn6330":  2,     # ISApl1-mcr-1-ISApl1 (IS30 TSD)
    "Tn2006":  9,     # ISAba1-blaOXA-23 (IS4 TSD)
    "Tn125":   3,     # ISAba125-blaNDM
    # Tn402 / Tn5053 family (class-1 integron progenitor)
    "Tn402":   5,
    "Tn5053":  5,
    # Tn552 (Staph; 6-7bp TSD)
    "Tn552":   6,
}


# ---------------------------------------------------------------------------
# ISEScan parser
# ---------------------------------------------------------------------------

def parse_isescan(tsv_path: str | Path) -> list[MGEFeature]:
    """
    Parse ISEScan .tsv output into MGEFeatures.

    ISEScan columns (tab-separated, no header on first data row):
      seqID family cluster isBegin isEnd isLen ncopy4is
      start1 end1 start2 end2 score irId irLen nGaps
      orfBegin orfEnd strand orfLen E-value E-value4copy type ov tir
    """
    features: list[MGEFeature] = []
    with open(tsv_path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            family = row["family"]
            strand_raw = row.get("strand", "").strip()
            strand = strand_raw if strand_raw in ("+", "-") else "."

            # TIR sequences come as "left:right" or "-:-"
            tir_raw = row.get("tir", "-:-")
            parts = tir_raw.split(":")
            ir_left = parts[0] if parts[0] != "-" else None
            ir_right = parts[1] if len(parts) > 1 and parts[1] != "-" else None

            score_raw = row.get("score", "")
            try:
                score = float(score_raw)
            except ValueError:
                score = None

            f = MGEFeature(
                element_type="IS",
                family=family,
                name=row["cluster"],
                start=int(row["isBegin"]),
                end=int(row["isEnd"]),
                strand=strand,
                tsd_length=TSD_LENGTHS.get(family),
                ir_left=ir_left,
                ir_right=ir_right,
                score=score,
                attributes={
                    "seqid": row.get("seqID", ""),
                    "type": row.get("type", ""),
                    "ncopy": row.get("ncopy4is", ""),
                },
            )
            features.append(f)
    return features


# ---------------------------------------------------------------------------
# AMRFinder+ parser
# ---------------------------------------------------------------------------

def parse_amrfinder(tsv_path: str | Path) -> list[MGEFeature]:
    """
    Parse AMRFinder+ TSV output.

    Relevant columns: Contig id, Start, Stop, Strand, Element symbol,
                      Element name, Class, Subclass
    """
    features: list[MGEFeature] = []
    with open(tsv_path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            strand_raw = row.get("Strand", "+").strip()
            strand = strand_raw if strand_raw in ("+", "-") else "."

            f = MGEFeature(
                element_type="AMR",
                family=row.get("Class", ""),
                name=row.get("Element symbol", ""),
                start=int(row["Start"]),
                end=int(row["Stop"]),
                strand=strand,
                attributes={
                    "seqid": row.get("Contig id", ""),
                    "full_name": row.get("Element name", ""),
                    "subclass": row.get("Subclass", ""),
                    "method": row.get("Method", ""),
                },
            )
            features.append(f)
    return features


# ---------------------------------------------------------------------------
# IntegronFinder parser
# ---------------------------------------------------------------------------

def parse_integron_finder(integrons_path: str | Path) -> list[MGEFeature]:
    """
    Parse IntegronFinder .integrons output.

    Groups rows by ID_integron. For each integron:
      - One MGEFeature(integron) spanning the full range
      - Child cassette/attC features attached to it

    Strand encoding: 1 → "+", -1 → "-"
    """
    rows: list[dict] = []
    with open(integrons_path) as fh:
        # Skip comment lines, then use the first non-comment line as header
        lines = [line for line in fh if not line.startswith("#")]
    reader = csv.DictReader(io.StringIO("".join(lines)), delimiter="\t")
    for row in reader:
        rows.append(dict(row))

    # Group by integron ID
    by_integron: dict[str, list[dict]] = {}
    for row in rows:
        by_integron.setdefault(row["ID_integron"], []).append(row)

    features: list[MGEFeature] = []
    for integron_id, members in by_integron.items():
        positions = [int(r["pos_beg"]) for r in members] + [int(r["pos_end"]) for r in members]
        int_start = min(positions)
        int_end = max(positions)

        # Determine integron class from the intI row
        int_type = members[0].get("type", "complete")
        family = f"{int_type}_integron"

        # Strand from intI row (first row with annotation == "intI")
        int_strand = "."
        for r in members:
            if r.get("annotation") == "intI":
                int_strand = _if_strand(r["strand"])
                break

        integron_seqid = members[0].get("ID_replicon", "")
        # Build a canonical cassette-array string (e.g. "|dfrA17|aadA5|")
        # by ordering protein cassettes along the integron and joining names.
        protein_rows = [r for r in members if r.get("type_elt") == "protein"
                        and r.get("annotation") != "intI"]
        protein_rows.sort(key=lambda r: int(r["pos_beg"]))
        cassette_array = "|".join(r["element"] for r in protein_rows)
        integron_feat = MGEFeature(
            element_type="integron",
            family=family,
            name=integron_id,
            start=int_start,
            end=int_end,
            strand=int_strand,
            attributes={
                "seqid": integron_seqid,
                "cassette_array": f"|{cassette_array}|" if cassette_array else "",
                "cassette_count": str(len(protein_rows)),
            },
        )

        for r in members:
            type_elt = r.get("type_elt", "")
            annotation = r.get("annotation", "")
            child_strand = _if_strand(r["strand"])

            if annotation == "intI":
                child = MGEFeature(
                    element_type="integrase",
                    family=family,
                    name=f"IntI_{integron_id}",
                    start=int(r["pos_beg"]),
                    end=int(r["pos_end"]),
                    strand=child_strand,
                )
            elif type_elt == "attC":
                child = MGEFeature(
                    element_type="attC",
                    family="attC",
                    name=r["element"],
                    start=int(r["pos_beg"]),
                    end=int(r["pos_end"]),
                    strand=child_strand,
                    attributes={"evalue": r.get("evalue", "")},
                )
            elif type_elt == "protein":
                child = MGEFeature(
                    element_type="cassette",
                    family="gene_cassette",
                    name=r["element"],
                    start=int(r["pos_beg"]),
                    end=int(r["pos_end"]),
                    strand=child_strand,
                )
            else:
                continue

            integron_feat.children.append(child)

        features.append(integron_feat)

    return features


def _if_strand(raw: str) -> str:
    """Convert IntegronFinder strand (1/-1) to +/-."""
    try:
        return "+" if int(raw) >= 0 else "-"
    except ValueError:
        return raw if raw in ("+", "-") else "."


# ---------------------------------------------------------------------------
# MOBsuite parser
# ---------------------------------------------------------------------------

def parse_mobsuite(report_path: str | Path) -> list[MGEFeature]:
    """
    Parse mob_typer TSV output (mob_typer_results.txt).

    Relevant columns: sample_id, rep_type(s), relaxase_type(s), mash_nearest_neighbor
    MOBsuite operates at replicon level — returns one MGEFeature per replicon hit.

    NOTE: mob_typer does not provide per-base coordinates; start/end are set to 0.
    Use mob_recon for coordinate-level replicon placement.
    """
    features: list[MGEFeature] = []
    with open(report_path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            rep_types = row.get("rep_type(s)", "").split(",")
            relaxases = row.get("relaxase_type(s)", "-").split(",")
            for rep in rep_types:
                rep = rep.strip()
                if not rep or rep == "-":
                    continue
                f = MGEFeature(
                    element_type="replicon",
                    family=rep,
                    name=rep,
                    start=0,
                    end=0,
                    strand=".",
                    attributes={
                        "relaxase": ",".join(relaxases),
                        "mash_neighbor": row.get("mash_nearest_neighbor", ""),
                    },
                )
                features.append(f)
    return features


# ---------------------------------------------------------------------------
# Tool runners (subprocess wrappers)
# ---------------------------------------------------------------------------

_ISESCAN_HEADER = (
    "seqID\tfamily\tcluster\tisBegin\tisEnd\tisLen\tncopy4is\t"
    "start1\tend1\tstart2\tend2\tscore\tirId\tirLen\tnGaps\t"
    "orfBegin\torfEnd\tstrand\torfLen\tE-value\tE-value4copy\ttype\tov\ttir\n"
)


def run_isescan(fasta: str | Path, outdir: str | Path) -> Path:
    """Run ISEScan via its pixi environment. Returns path to .tsv output.

    ISEScan produces no output file when it finds no IS elements or when the
    input sequence is too short (< ~1 kb). In both cases we create a
    header-only TSV so downstream parsers always receive a valid file.
    """
    import warnings
    fasta = Path(fasta)
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        ["pixi", "run", "-e", "isescan", "isescan.py",
         "--seqfile", str(fasta), "--output", str(outdir), "--nthread", "4"],
        check=True,
    )
    tsv = outdir / fasta.name / (fasta.name + ".tsv")
    if not tsv.exists():
        candidates = list(outdir.rglob("*.tsv"))
        if candidates:
            tsv = candidates[0]
    if not tsv.exists():
        warnings.warn(
            f"ISEScan produced no output for {fasta.name} — "
            "sequence may be too short or contain no IS elements. "
            "Writing empty result.",
            RuntimeWarning,
            stacklevel=2,
        )
        tsv = outdir / (fasta.stem + ".isescan.tsv")
        tsv.write_text(_ISESCAN_HEADER)
    return tsv


def run_amrfinder(fasta: str | Path, outdir: str | Path) -> Path:
    """Run AMRFinder+ via its pixi environment. Returns path to TSV output."""
    fasta = Path(fasta)
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    out_tsv = outdir / (fasta.stem + ".tsv")
    subprocess.run(
        ["pixi", "run", "-e", "amrfinder", "amrfinder",
         "-n", str(fasta), "-o", str(out_tsv), "--plus"],
        check=True,
    )
    return out_tsv


def run_integron_finder(fasta: str | Path, outdir: str | Path) -> Path:
    """Run IntegronFinder via its pixi environment. Returns path to .integrons file."""
    fasta = Path(fasta)
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        ["pixi", "run", "-e", "integron", "integron_finder",
         "--outdir", str(outdir), "--circ", str(fasta)],
        check=True,
    )
    candidates = list(outdir.rglob("*.integrons"))
    if not candidates:
        raise FileNotFoundError(f"IntegronFinder produced no .integrons file in {outdir}")
    return candidates[0]


def run_mobtyper(fasta: str | Path, outdir: str | Path) -> Path:
    """Run MOBsuite mob_typer via its pixi environment. Returns path to results TSV."""
    fasta = Path(fasta)
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    out_tsv = outdir / "mob_typer_results.txt"
    subprocess.run(
        ["pixi", "run", "-e", "mobsuite", "mob_typer",
         "--infile", str(fasta), "--out_file", str(out_tsv)],
        check=True,
    )
    return out_tsv
