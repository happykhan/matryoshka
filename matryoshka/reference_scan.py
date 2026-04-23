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
    """Parse 'name ...key=val key=val...' style FASTA headers.

    Also handles PlasmidFinder-style headers
    (``>RepName_variant__accession``) by synthesising element_type,
    family and name fields.
    """
    parts = header.split()
    meta: dict[str, str] = {"_id": parts[0]}

    # PlasmidFinder-style: no key=value metadata, seq id encodes it
    if len(parts) == 1 and "__" in parts[0]:
        seq_id = parts[0]
        # e.g. "IncHI1B(R27)_1_R27_AF250878" or "pKPC-CAV1321_1__CP011611"
        name_part, _, accession = seq_id.rpartition("__")
        # Strip trailing _N variant number
        canonical = name_part.rsplit("_", 1)[0] if name_part.rsplit("_", 1)[-1].isdigit() else name_part
        meta.update({
            "element_type": "replicon",
            "family": canonical,
            "name": canonical,
            "source_accession": accession,
            "source": "reference_scan",
        })
        return meta

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


def _score_variant_hits(hits: list[BlastHit]) -> dict[str, dict]:
    """Aggregate per-subject coverage across multiple HSPs.

    Returns a dict keyed by sseqid with:
        {total_align, slen, coverage_pct, hsp_count, max_pident, span_start,
         span_end, qstart_min, qend_max}

    Used for variant discrimination (Tn4401 a/b/c/d/…) where variants
    differ by short internal deletions. A single-HSP full-length match
    is the cleanest signal that query matches a given variant exactly.
    """
    by_sub: dict[str, dict] = {}
    for h in hits:
        entry = by_sub.setdefault(h.sseqid, {
            "total_align": 0,
            "slen": h.slen,
            "hsp_count": 0,
            "max_pident": 0.0,
            "qseqid": h.qseqid,
            "qstart_min": h.qstart,
            "qend_max": h.qend,
        })
        entry["total_align"] += h.length
        entry["hsp_count"] += 1
        entry["max_pident"] = max(entry["max_pident"], h.pident)
        entry["qstart_min"] = min(entry["qstart_min"], h.qstart, h.qend)
        entry["qend_max"] = max(entry["qend_max"], h.qstart, h.qend)
    for v in by_sub.values():
        v["coverage_pct"] = 100.0 * v["total_align"] / max(v["slen"], 1)
    return by_sub


def _pick_best_variant(hits: list[BlastHit]) -> BlastHit | None:
    """Given hits that all target Tn4401 variants overlapping the same query
    region, pick the best-fitting variant.

    Best = highest coverage, tie-break by fewest HSPs (cleaner match),
    then highest identity.
    """
    if not hits:
        return None
    scores = _score_variant_hits(hits)
    best_sub = max(
        scores,
        key=lambda s: (
            scores[s]["coverage_pct"],
            -scores[s]["hsp_count"],
            scores[s]["max_pident"],
        ),
    )
    # Return the top HSP for the winning subject (boundary annotation)
    return max(
        (h for h in hits if h.sseqid == best_sub),
        key=lambda h: h.length,
    )


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
    # Tn4401: when the best-matching variant still has < 98% subject
    # coverage we're looking at an unseen variant (c/d/e/f/g/h) with an
    # internal deletion not represented in the bundled reference set.
    if family == "Tn4401":
        sub_cov_s = meta.get("blast_subject_coverage")
        sub_cov = float(sub_cov_s) if sub_cov_s else hit.qcovs
        if sub_cov < 98.0:
            deletion_bp = max(0, hit.slen - hit.length)
            attrs["variant"] = "unknown"
            attrs["deletion_bp"] = deletion_bp
            attrs["note"] = (
                f"best-match coverage {sub_cov:.1f}% of known variants — "
                f"likely uncatalogued variant with ~{deletion_bp}bp deletion"
            )
            name = "Tn4401"

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


def _group_by_query_region(
    hits: list[BlastHit],
    overlap_frac: float = 0.5,
) -> list[list[BlastHit]]:
    """Cluster hits that overlap on the query, regardless of subject.

    Tn4401 variants are highly similar to each other, so one Tn4401 in
    the query produces near-identical hits against multiple variant
    references. We want to pick the *best* matching variant for each
    query region, not emit all of them.
    """
    if not hits:
        return []
    normed = [
        (min(h.qstart, h.qend), max(h.qstart, h.qend), h)
        for h in hits
    ]
    normed.sort()
    groups: list[list[BlastHit]] = []
    current: list[tuple[int, int, BlastHit]] = [normed[0]]
    cur_start, cur_end = normed[0][0], normed[0][1]
    for s, e, h in normed[1:]:
        ovl = max(0, min(cur_end, e) - max(cur_start, s))
        short = min(cur_end - cur_start, e - s)
        if short > 0 and ovl / short > overlap_frac:
            current.append((s, e, h))
            cur_end = max(cur_end, e)
        else:
            groups.append([t[2] for t in current])
            current = [(s, e, h)]
            cur_start, cur_end = s, e
    groups.append([t[2] for t in current])
    return groups


def scan(
    query_fasta: Path | str,
    reference_fasta: Path | str,
    min_identity: float = 85.0,
    min_length: int = 500,
    min_subject_coverage: float = 0.0,
    pick_best_variant: bool = False,
) -> list[MGEFeature]:
    """Run a BLAST reference scan and return emitted MGEFeatures.

    If `pick_best_variant` is True (used for Tn4401 a/b/c/d/e/f/g/h),
    hits that overlap on the query are clustered and only the best-
    fitting variant is emitted per cluster.
    """
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

    if pick_best_variant:
        features: list[MGEFeature] = []
        for group in _group_by_query_region(hits):
            best = _pick_best_variant(group)
            if best is None:
                continue
            group_meta = dict(meta.get(best.sseqid, {"_id": best.sseqid}))
            # Enrich with aggregate variant stats for downstream reporting
            scores = _score_variant_hits(group)[best.sseqid]
            group_meta["blast_subject_coverage"] = f"{scores['coverage_pct']:.1f}"
            group_meta["blast_hsp_count"] = str(scores["hsp_count"])
            features.append(_hit_to_feature(best, group_meta))
        return features

    # Default path: per-HSP dedup (good for distinct non-homologous references)
    hits = _dedupe_overlapping(hits)
    return [_hit_to_feature(h, meta.get(h.sseqid, {"_id": h.sseqid})) for h in hits]


# Per-reference-file search parameters. Short motifs (res site) need low
# min_length; whole-plasmid exemplars use moderate length thresholds.
REFERENCE_PARAMS: dict[str, dict] = {
    "tn3_res_site.fasta":       {"min_identity": 75.0, "min_length": 80},
    "tn4401.fasta":             {"min_identity": 95.0, "min_length": 1_000,
                                 "pick_best_variant": True},
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
    # Tn3-family extras — Tn3, Tn2, Tn1696, Tn1721, Tn6452
    "tn3_family_extras.fasta":  {"min_identity": 90.0, "min_length": 3_000},
    # PlasmidFinder replicon sequences are ~500-900bp per record.
    # Identity ≥95% is the standard PlasmidFinder cutoff.
    "plasmidfinder_enterobacteriales.fasta": {
        "min_identity": 95.0, "min_length": 200,
        "min_subject_coverage": 60.0,
    },
    # ter-site motifs for IS91 / ISCR rolling-circle elements.
    # Short consensus (~30bp) — needs a broad min_identity tolerance.
    # Marked as experimental pending a validated motif library.
    "rolling_circle_ter_sites.fasta": {
        "min_identity": 85.0, "min_length": 20,
    },
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
