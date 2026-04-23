"""
reference_scan.py — Homology-based detection of known MGE structures via BLAST.

Some MGEs can't be called purely from IS+AMR patterns — they require
sequence homology to reference exemplars. This module runs blastn of the
input against bundled reference FASTAs and emits MGEFeatures from hits
above configurable identity/coverage thresholds.

Reference FASTAs live in ``matryoshka/references/*.fasta``. Each record
header encodes metadata as key=value pairs after the accession, e.g.::

    >Tn4401a EU176011 element_type=transposon family=Tn4401 variant=a
    >Tn3_res_site element_type=res_site family=Tn3

Parsed headers drive how hits map onto MGEFeature fields.
"""

from __future__ import annotations

import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path

from Bio import SeqIO

from .detect import MGEFeature


REFERENCES_DIR = Path(__file__).parent / "references"


@dataclass(frozen=True)
class BlastHit:
    qseqid: str      # query contig id
    sseqid: str      # subject (reference) id
    pident: float
    length: int
    qstart: int      # 1-based
    qend: int
    sstart: int
    send: int
    evalue: float
    qlen: int
    slen: int

    @property
    def qcovs(self) -> float:
        """Percent of subject covered by this HSP."""
        return 100.0 * self.length / max(self.slen, 1)


def blast_available() -> bool:
    return shutil.which("blastn") is not None


def _parse_header(header: str) -> dict[str, str]:
    """Parse 'name ...key=val key=val...' style FASTA headers."""
    parts = header.split()
    meta: dict[str, str] = {"_id": parts[0]}
    for p in parts[1:]:
        if "=" in p:
            k, v = p.split("=", 1)
            meta[k] = v
    return meta


def _load_reference_metadata(fasta: Path) -> dict[str, dict[str, str]]:
    """Map reference seq id -> metadata dict parsed from header."""
    meta: dict[str, dict[str, str]] = {}
    for rec in SeqIO.parse(fasta, "fasta"):
        parsed = _parse_header(rec.description)
        parsed["_length"] = str(len(rec.seq))
        meta[rec.id] = parsed
    return meta


def blastn_hits(
    query_fasta: Path | str,
    reference_fasta: Path | str,
    min_identity: float = 85.0,
    min_length: int = 500,
    min_subject_coverage: float = 0.0,
    evalue: float = 1e-10,
) -> list[BlastHit]:
    """Run blastn and return high-confidence hits.

    A hit is kept iff:
      pident >= min_identity
      AND length >= min_length
      AND (subject_coverage >= min_subject_coverage   OR  min_subject_coverage == 0)

    When references are "element-sized" (e.g. a single Tn), a high
    min_subject_coverage works well. When references are whole plasmids
    or chromosomes containing the element, use min_length only (set
    min_subject_coverage=0).
    """
    if not blast_available():
        raise RuntimeError("blastn not found on PATH")

    reference_fasta = Path(reference_fasta)
    with tempfile.TemporaryDirectory() as tmp:
        db_prefix = Path(tmp) / "refdb"
        subprocess.run(
            ["makeblastdb", "-in", str(reference_fasta),
             "-dbtype", "nucl", "-out", str(db_prefix), "-parse_seqids"],
            check=True, capture_output=True,
        )
        outfmt = (
            "6 qseqid sseqid pident length qstart qend sstart send "
            "evalue qlen slen"
        )
        result = subprocess.run(
            ["blastn", "-query", str(query_fasta), "-db", str(db_prefix),
             "-outfmt", outfmt, "-evalue", str(evalue),
             "-perc_identity", str(min_identity), "-num_threads", "2"],
            check=True, capture_output=True, text=True,
        )

    hits: list[BlastHit] = []
    for line in result.stdout.splitlines():
        cols = line.split("\t")
        if len(cols) < 11:
            continue
        h = BlastHit(
            qseqid=cols[0], sseqid=cols[1],
            pident=float(cols[2]), length=int(cols[3]),
            qstart=int(cols[4]), qend=int(cols[5]),
            sstart=int(cols[6]), send=int(cols[7]),
            evalue=float(cols[8]),
            qlen=int(cols[9]), slen=int(cols[10]),
        )
        if h.length < min_length:
            continue
        if min_subject_coverage > 0 and h.qcovs < min_subject_coverage:
            continue
        hits.append(h)
    return hits


def _hit_to_feature(hit: BlastHit, meta: dict[str, str]) -> MGEFeature:
    start, end = sorted((hit.qstart, hit.qend))
    strand = "+" if hit.qstart <= hit.qend else "-"
    element_type = meta.get("element_type", "mobile_element")
    family = meta.get("family", meta.get("_id", ""))
    name = meta.get("name", meta.get("_id", ""))
    attrs: dict[str, object] = {
        "seqid": hit.qseqid,
        "reference_id": hit.sseqid,
        "blast_identity": round(hit.pident, 2),
        "blast_coverage": round(hit.qcovs, 2),
        "source": "reference_scan",
    }
    # Tn4401: if coverage is < 98% of the Tn4401b reference, a known
    # variant-specific internal deletion is likely — flag that we can't
    # disambiguate a/b/c/d/e/f/g/h without the full variant set.
    if family == "Tn4401" and hit.qcovs < 98.0:
        deletion_bp = max(0, hit.slen - hit.length)
        attrs["variant"] = "unknown"
        attrs["deletion_bp"] = deletion_bp
        attrs["note"] = (
            "coverage <98% of Tn4401b — possible variant a/c/d/e/f/g/h "
            f"with ~{deletion_bp}bp internal deletion"
        )
        name = "Tn4401"  # drop the 'b' since we can't confirm

    for k, v in meta.items():
        if k.startswith("_") or k in ("element_type", "family", "name"):
            continue
        attrs.setdefault(k, v)

    return MGEFeature(
        element_type=element_type,
        family=family,
        name=name,
        start=start,
        end=end,
        strand=strand,
        score=hit.pident,
        attributes=attrs,
    )


def _dedupe_overlapping(hits: list[BlastHit]) -> list[BlastHit]:
    """Keep the best hit per overlapping query interval.

    When multiple references overlap on the same region (e.g. Tn4401a and
    Tn4401b share 99% of sequence), report only the highest-identity one.
    """
    if not hits:
        return []
    # Normalise orientation and sort by identity desc, length desc
    ordered = sorted(
        hits,
        key=lambda h: (-h.pident, -h.length),
    )
    kept: list[BlastHit] = []
    for h in ordered:
        qs, qe = sorted((h.qstart, h.qend))
        overlaps = False
        for k in kept:
            ks, ke = sorted((k.qstart, k.qend))
            # Overlap if ranges intersect by more than 50% of the shorter
            ovl = max(0, min(qe, ke) - max(qs, ks))
            short_len = min(qe - qs, ke - ks)
            if short_len > 0 and ovl / short_len > 0.5:
                overlaps = True
                break
        if not overlaps:
            kept.append(h)
    return kept


def scan(
    query_fasta: Path | str,
    reference_fasta: Path | str,
    min_identity: float = 85.0,
    min_length: int = 500,
    min_subject_coverage: float = 0.0,
) -> list[MGEFeature]:
    """Run a BLAST reference scan and return emitted MGEFeatures."""
    reference_fasta = Path(reference_fasta)
    if not reference_fasta.exists():
        return []
    meta = _load_reference_metadata(reference_fasta)
    hits = blastn_hits(
        query_fasta, reference_fasta,
        min_identity=min_identity,
        min_length=min_length,
        min_subject_coverage=min_subject_coverage,
    )
    hits = _dedupe_overlapping(hits)
    return [_hit_to_feature(h, meta.get(h.sseqid, {"_id": h.sseqid})) for h in hits]


# Per-reference-file search parameters. Short motifs (res site) need low
# min_length; whole-plasmid exemplars use moderate length thresholds.
REFERENCE_PARAMS: dict[str, dict] = {
    "tn3_res_site.fasta":       {"min_identity": 75.0, "min_length": 80},
    "tn4401.fasta":             {"min_identity": 95.0, "min_length": 4_000},
    "tn1546.fasta":             {"min_identity": 95.0, "min_length": 8_000},
    "tn21.fasta":               {"min_identity": 90.0, "min_length": 5_000},
    "tn1331.fasta":             {"min_identity": 95.0, "min_length": 5_000},
    "tn5393.fasta":             {"min_identity": 95.0, "min_length": 3_000},
    "acinetobacter_islands.fasta": {"min_identity": 95.0, "min_length": 10_000},
    "mcr1_exemplars.fasta":     {"min_identity": 97.0, "min_length": 5_000},
    "integron_archetypes.fasta":{"min_identity": 90.0, "min_length": 2_000},
    "gi_sul2.fasta":            {"min_identity": 95.0, "min_length": 5_000},
    "isecp1.fasta":             {"min_identity": 90.0, "min_length": 1_000},
    "tn7.fasta":                {"min_identity": 90.0, "min_length": 3_000},
    "tn552.fasta":              {"min_identity": 95.0, "min_length": 2_000},
}


def _infer_acinetobacter_islands(hits: list[MGEFeature]) -> list[MGEFeature]:
    """Detect AbGRI1 (Tn6022 + Tn6172 flanking comM insertion) and AbaR3-like
    islands (Tn6019 backbone) from the Acinetobacter Tn hit set.

    Emits one parent MGEFeature("genomic_island") per detected island.
    """
    t6022 = [h for h in hits if h.family == "Tn6022"]
    t6172 = [h for h in hits if h.family == "Tn6172"]
    t6019 = [h for h in hits if h.family == "Tn6019"]

    out: list[MGEFeature] = []

    # AbGRI1 = Tn6022 + Tn6172 within ~20kb of each other on the same contig
    for a in t6022:
        for b in t6172:
            if a.attributes.get("seqid") != b.attributes.get("seqid"):
                continue
            gap = max(a.start, b.start) - min(a.end, b.end)
            if 0 <= gap <= 20_000:
                start = min(a.start, b.start)
                end = max(a.end, b.end)
                out.append(MGEFeature(
                    element_type="genomic_island",
                    family="AbGRI1",
                    name="AbGRI1",
                    start=start, end=end, strand=".",
                    attributes={
                        "seqid": a.attributes.get("seqid", ""),
                        "source": "reference_scan",
                        "note": "Tn6022 + Tn6172 flanking pattern — comM-inserted",
                    },
                ))

    # AbaR3-like: Tn6019 backbone alone (no other marker needed)
    for a in t6019:
        out.append(MGEFeature(
            element_type="genomic_island",
            family="AbaR3",
            name="AbaR3-like",
            start=a.start, end=a.end, strand=a.strand,
            attributes={
                "seqid": a.attributes.get("seqid", ""),
                "source": "reference_scan",
                "note": "Tn6019 backbone — AbaR3-family resistance island",
            },
        ))

    return out


def scan_all(
    query_fasta: Path | str,
    references_dir: Path = REFERENCES_DIR,
) -> list[MGEFeature]:
    """Scan every .fasta under references_dir with per-file parameters.

    Adds post-processing: if Tn6022+Tn6172 or Tn6019 hits are present,
    emit a parent AbGRI1 / AbaR3 genomic_island feature.
    """
    out: list[MGEFeature] = []
    for ref in sorted(references_dir.glob("*.fasta")):
        params = REFERENCE_PARAMS.get(ref.name, {"min_identity": 85.0, "min_length": 500})
        out.extend(scan(query_fasta, ref, **params))
    # Compound-island detection (runs only over the hit set we just made)
    out.extend(_infer_acinetobacter_islands(out))
    return out
