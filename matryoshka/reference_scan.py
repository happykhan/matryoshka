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
from collections.abc import Callable
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path

from Bio import SeqIO

from .detect import TSD_LENGTHS, MGEFeature
from .tn123 import REFERENCE_METADATA

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


@dataclass(frozen=True)
class Tn123Assembly:
    """A collinear set of BLAST HSPs representing one candidate Tn locus."""

    qseqid: str
    sseqid: str
    strand: str
    qstart: int
    qend: int
    slen: int
    identity: float
    reference_coverage: float
    hsp_count: int
    query_span: int
    aligned_query_bases: int
    inserted_bases: int
    deleted_bases: int
    left_end_covered: bool
    right_end_covered: bool
    structural_status: str
    hits: tuple[BlastHit, ...]


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
        if fasta.name == "tn1_tn2_tn3.fasta" and rec.id in REFERENCE_METADATA:
            parsed.update(REFERENCE_METADATA[rec.id])
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
    num_threads: int = 1,
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
        task = "blastn-short" if min_length < 50 else "blastn"
        result = subprocess.run(
            ["blastn", "-task", task, "-query", str(query_fasta), "-db", str(db_prefix),
             "-outfmt", outfmt, "-evalue", str(evalue),
             "-perc_identity", str(min_identity), "-num_threads", str(max(1, num_threads))],
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
        v["coverage_pct"] = min(
            100.0,
            100.0 * v["total_align"] / max(v["slen"], 1),
        )
    return by_sub


def _pick_best_variant(
    hits: list[BlastHit],
    prefer_identity: bool = False,
) -> BlastHit | None:
    """Given hits that all target Tn4401 variants overlapping the same query
    region, pick the best-fitting variant.

    Best = highest coverage, tie-break by fewest HSPs (cleaner match),
    then highest identity.
    """
    if not hits:
        return None
    scores = _score_variant_hits(hits)
    if prefer_identity:
        best_sub = max(
            scores,
            key=lambda s: (
                scores[s]["max_pident"],
                -abs(100.0 - scores[s]["coverage_pct"]),
                -scores[s]["hsp_count"],
                scores[s]["total_align"],
            ),
        )
    else:
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
    strand = _relative_strand(hit)
    element_type = meta.get("element_type", "mobile_element")
    family = meta.get("family", meta.get("_id", ""))
    name = meta.get("name", meta.get("_id", ""))
    attrs: dict[str, object] = {
        "seqid": hit.qseqid,
        "reference_id": hit.sseqid,
        "blast_identity": round(hit.pident, 2),
        "blast_coverage": round(hit.qcovs, 2),
        "reference_length": hit.slen,
        "source": "reference_scan",
    }
    if element_type in {"IS", "IR", "integron_segment", "transposon_component"}:
        attrs["type"] = "c" if hit.qcovs >= 95.0 else "p"
        if hit.qcovs < 95.0:
            attrs["fragment"] = True

    if family == "ISEcp1_TPU":
        aggregate_cov_s = meta.get("blast_subject_coverage")
        aggregate_cov = float(aggregate_cov_s) if aggregate_cov_s else hit.qcovs
        complete = aggregate_cov >= 98.0 and hit.pident >= 98.0
        attrs["structural_status"] = (
            "complete_reference_match" if complete else "fragment_only"
        )
        attrs["type"] = "c" if complete else "p"
        if not complete:
            attrs["fragment"] = True
            attrs["best_match"] = name
            attrs["note"] = (
                "ISEcp1-associated homology is incomplete; IRalt and outer "
                "direct-repeat boundaries are not established"
            )
            name = "ISEcp1-associated_region#"
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

    if meta.get("tn123_canonical") == "true":
        sub_cov_s = meta.get("blast_subject_coverage")
        sub_cov = float(sub_cov_s) if sub_cov_s else hit.qcovs
        if sub_cov < 98.0 or hit.pident < 98.0:
            attrs["best_match"] = family
            attrs["fragment"] = True
            attrs["note"] = (
                "partial or divergent homology cannot reliably distinguish Tn1, Tn2 and Tn3"
            )
            family = "Tn3_family"
            name = "Tn1/2/3"

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
        tsd_length=(
            int(meta["tsd_length"])
            if meta.get("tsd_length", "").isdigit()
            else TSD_LENGTHS.get(family)
        ),
        attributes=attrs,
    )


def _dedupe_overlapping(hits: list[BlastHit]) -> list[BlastHit]:
    """Keep the best hit per overlapping query interval.

    When multiple references overlap on the same region (e.g. Tn4401a and
    Tn4401b share 99% of sequence), report only the highest-identity one.
    """
    if not hits:
        return []
    # Normalise orientation and sort by identity desc, subject coverage desc,
    # length desc — prefer hits that cover more of the reference sequence
    ordered = sorted(
        hits,
        key=lambda h: (-h.pident, -h.qcovs, -h.length),
    )
    kept: list[BlastHit] = []
    for h in ordered:
        qs, qe = sorted((h.qstart, h.qend))
        overlaps = False
        for k in kept:
            # Only compare hits on the same query contig
            if k.qseqid != h.qseqid:
                continue
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
    # Query records may all begin at coordinate 1 (a multi-FASTA reference
    # corpus is a common example). Keep records contiguous before clustering;
    # otherwise interleaved coordinates create multiple false loci per record.
    normed.sort(key=lambda t: (t[2].qseqid, t[0], t[1]))
    groups: list[list[BlastHit]] = []
    current: list[tuple[int, int, BlastHit]] = [normed[0]]
    cur_start, cur_end = normed[0][0], normed[0][1]
    cur_query = normed[0][2].qseqid
    for s, e, h in normed[1:]:
        ovl = max(0, min(cur_end, e) - max(cur_start, s))
        short = min(cur_end - cur_start, e - s)
        if h.qseqid == cur_query and short > 0 and ovl / short > overlap_frac:
            current.append((s, e, h))
            cur_end = max(cur_end, e)
        else:
            groups.append([t[2] for t in current])
            current = [(s, e, h)]
            cur_start, cur_end = s, e
            cur_query = h.qseqid
    groups.append([t[2] for t in current])
    return groups


def _relative_strand(hit: BlastHit) -> str:
    """Orientation of the query locus relative to the reference sequence."""
    query_step = 1 if hit.qend >= hit.qstart else -1
    subject_step = 1 if hit.send >= hit.sstart else -1
    return "+" if query_step == subject_step else "-"


def _query_interval(hit: BlastHit) -> tuple[int, int]:
    return min(hit.qstart, hit.qend), max(hit.qstart, hit.qend)


def _subject_axis_interval(hit: BlastHit, strand: str) -> tuple[int, int]:
    """Reference interval transformed so it increases along query coordinates."""
    start, end = sorted((hit.sstart, hit.send))
    if strand == "+":
        return start, end
    return hit.slen - end + 1, hit.slen - start + 1


def _union_length(intervals: list[tuple[int, int]]) -> int:
    if not intervals:
        return 0
    ordered = sorted((min(start, end), max(start, end)) for start, end in intervals)
    total = 0
    current_start, current_end = ordered[0]
    for start, end in ordered[1:]:
        if start <= current_end + 1:
            current_end = max(current_end, end)
        else:
            total += current_end - current_start + 1
            current_start, current_end = start, end
    return total + current_end - current_start + 1


def _tn123_chain_gaps(
    hits: list[BlastHit],
    strand: str,
) -> tuple[int, int]:
    """Return inserted and deleted bases between consecutive HSPs."""
    ordered = sorted(hits, key=_query_interval)
    inserted = 0
    deleted = 0
    for left, right in zip(ordered, ordered[1:], strict=False):
        left_qend = _query_interval(left)[1]
        right_qstart = _query_interval(right)[0]
        left_send = _subject_axis_interval(left, strand)[1]
        right_sstart = _subject_axis_interval(right, strand)[0]
        query_gap = max(0, right_qstart - left_qend - 1)
        reference_gap = max(0, right_sstart - left_send - 1)
        inserted += max(0, query_gap - reference_gap)
        deleted += max(0, reference_gap - query_gap)
    return inserted, deleted


def _chain_tn123_subject_hits(
    hits: list[BlastHit],
    max_query_insertion: int = 50_000,
    max_reference_gap: int = 500,
    max_reference_overlap: int = 100,
) -> list[list[BlastHit]]:
    """Chain collinear HSPs while keeping separate copies as separate loci."""
    grouped: dict[tuple[str, str, str], list[BlastHit]] = {}
    for hit in hits:
        key = (hit.qseqid, hit.sseqid, _relative_strand(hit))
        grouped.setdefault(key, []).append(hit)

    chains: list[list[BlastHit]] = []
    for (_qseqid, _sseqid, strand), subject_hits in grouped.items():
        ordered = sorted(subject_hits, key=lambda item: _query_interval(item))
        current: list[BlastHit] = []
        current_qend = 0
        current_send = 0
        for hit in ordered:
            qstart, qend = _query_interval(hit)
            sstart, send = _subject_axis_interval(hit, strand)
            if not current:
                current = [hit]
                current_qend = qend
                current_send = send
                continue

            query_gap = qstart - current_qend - 1
            reference_gap = sstart - current_send - 1
            collinear = (
                query_gap <= max_query_insertion
                and -max_reference_overlap <= reference_gap <= max_reference_gap
            )
            if collinear:
                current.append(hit)
                current_qend = max(current_qend, qend)
                current_send = max(current_send, send)
            else:
                chains.append(current)
                current = [hit]
                current_qend = qend
                current_send = send
        if current:
            chains.append(current)
    return chains


def _assemble_tn123_chain(hits: list[BlastHit]) -> Tn123Assembly:
    first = hits[0]
    strand = _relative_strand(first)
    q_intervals = [_query_interval(hit) for hit in hits]
    s_intervals: list[tuple[int, int]] = [
        (min(hit.sstart, hit.send), max(hit.sstart, hit.send)) for hit in hits
    ]
    qstart = min(start for start, _ in q_intervals)
    qend = max(end for _, end in q_intervals)
    subject_start = min(start for start, _ in s_intervals)
    subject_end = max(end for _, end in s_intervals)
    reference_bases = _union_length(s_intervals)
    query_bases = _union_length(q_intervals)
    aligned_weight = sum(hit.length for hit in hits)
    identity = (
        sum(hit.pident * hit.length for hit in hits) / aligned_weight
        if aligned_weight
        else 0.0
    )
    query_span = qend - qstart + 1
    left_covered = subject_start <= 50
    right_covered = subject_end >= first.slen - 49
    inserted_bases, deleted_bases = _tn123_chain_gaps(hits, strand)
    coverage = min(100.0, 100.0 * reference_bases / max(first.slen, 1))

    if left_covered and right_covered:
        if inserted_bases > 20 and deleted_bases > 20:
            status = "complete_with_indel"
        elif inserted_bases > 20:
            status = "complete_with_insertion"
        elif deleted_bases > 20:
            status = "complete_with_deletion"
        else:
            status = "intact"
    elif left_covered:
        status = "right_partial"
    elif right_covered:
        status = "left_partial"
    else:
        status = "internal_fragment"

    return Tn123Assembly(
        qseqid=first.qseqid,
        sseqid=first.sseqid,
        strand=strand,
        qstart=qstart,
        qend=qend,
        slen=first.slen,
        identity=identity,
        reference_coverage=coverage,
        hsp_count=len(hits),
        query_span=query_span,
        aligned_query_bases=query_bases,
        inserted_bases=inserted_bases,
        deleted_bases=deleted_bases,
        left_end_covered=left_covered,
        right_end_covered=right_covered,
        structural_status=status,
        hits=tuple(sorted(hits, key=lambda item: _query_interval(item))),
    )


def _group_tn123_assemblies(
    assemblies: list[Tn123Assembly],
) -> list[list[Tn123Assembly]]:
    """Group competing Tn1/Tn2/Tn3 reference assemblies at each query locus."""
    ordered = sorted(assemblies, key=lambda item: (item.qseqid, item.qstart, item.qend))
    groups: list[list[Tn123Assembly]] = []
    for assembly in ordered:
        for group in reversed(groups):
            representative = group[0]
            if representative.qseqid != assembly.qseqid:
                continue
            overlap = max(
                0,
                min(representative.qend, assembly.qend)
                - max(representative.qstart, assembly.qstart)
                + 1,
            )
            shorter = min(
                representative.qend - representative.qstart + 1,
                assembly.qend - assembly.qstart + 1,
            )
            if shorter and overlap / shorter >= 0.5:
                group.append(assembly)
                break
        else:
            groups.append([assembly])
    return groups


def _reference_segments(assembly: Tn123Assembly) -> list[dict[str, int | float]]:
    return [
        {
            "qstart": hit.qstart,
            "qend": hit.qend,
            "sstart": hit.sstart,
            "send": hit.send,
            "identity": round(hit.pident, 3),
        }
        for hit in assembly.hits
    ]


def _tn123_assembly_to_feature(
    assembly: Tn123Assembly,
    runner_up: Tn123Assembly | None,
    meta: dict[str, str],
) -> MGEFeature:
    family = meta.get("family", "Tn3_family")
    margin = (
        assembly.identity - runner_up.identity
        if runner_up is not None
        else assembly.identity
    )
    complete = assembly.left_end_covered and assembly.right_end_covered
    clear_nearest = margin >= 0.5
    exact = (
        clear_nearest
        and complete
        and assembly.structural_status == "intact"
        and assembly.identity >= 99.8
        and assembly.reference_coverage >= 99.5
    )
    named_variant = (
        clear_nearest
        and assembly.identity >= 95.0
        and assembly.reference_coverage >= 80.0
    )

    if exact:
        name = family
        variant_status = "exact_reference"
    elif named_variant:
        name = f"{family}-like"
        variant_status = (
            "minor_sequence_variant"
            if complete and assembly.identity >= 98.0
            else assembly.structural_status
        )
    else:
        family = "Tn3_family"
        name = "Tn1/2/3" if complete else "Tn1/2/3 fragment"
        variant_status = "ambiguous"

    fragment = not complete or assembly.reference_coverage < 80.0
    attrs: dict[str, object] = {
        "seqid": assembly.qseqid,
        "reference_id": assembly.sseqid,
        "source": "reference_scan",
        "source_accession": meta.get("source_accession", ""),
        "tn123_canonical": "true",
        "best_match": meta.get("family", ""),
        "blast_identity": round(assembly.identity, 2),
        "blast_coverage": round(assembly.reference_coverage, 2),
        "blast_subject_coverage": round(assembly.reference_coverage, 2),
        "blast_hsp_count": assembly.hsp_count,
        "reference_length": assembly.slen,
        "reference_segments": _reference_segments(assembly),
        "variant_margin": round(margin, 2),
        "variant_status": variant_status,
        "structural_status": assembly.structural_status,
        "query_span": assembly.query_span,
        "inserted_bases": assembly.inserted_bases,
        "deleted_bases": assembly.deleted_bases,
        "left_end_covered": assembly.left_end_covered,
        "right_end_covered": assembly.right_end_covered,
    }
    if fragment:
        attrs["fragment"] = True
    if not clear_nearest:
        attrs["note"] = "Tn1, Tn2 and Tn3 reference scores are too close to distinguish"
    elif not exact:
        attrs["note"] = f"closest reference is {meta.get('family', '')}; reported as a variant"

    return MGEFeature(
        element_type="transposon",
        family=family,
        name=name,
        start=assembly.qstart,
        end=assembly.qend,
        strand=assembly.strand,
        score=assembly.identity,
        attributes=attrs,
    )


def _scan_tn123_hits(
    hits: list[BlastHit],
    meta: dict[str, dict[str, str]],
) -> list[MGEFeature]:
    chains = _chain_tn123_subject_hits(hits)
    assemblies = [_assemble_tn123_chain(chain) for chain in chains]
    features: list[MGEFeature] = []
    for group in _group_tn123_assemblies(assemblies):
        ranked = sorted(
            group,
            key=lambda item: (
                item.identity,
                item.reference_coverage,
                item.left_end_covered and item.right_end_covered,
                -item.hsp_count,
            ),
            reverse=True,
        )
        best = ranked[0]
        runner_up = next(
            (item for item in ranked[1:] if item.sseqid != best.sseqid),
            None,
        )
        features.append(_tn123_assembly_to_feature(
            best,
            runner_up,
            meta.get(best.sseqid, {"_id": best.sseqid}),
        ))
    return features


def _generic_assembly_to_feature(
    assembly: Tn123Assembly,
    meta: dict[str, str],
) -> MGEFeature:
    """Convert a collinear assembly for a curated non-Tn1/2/3 unit."""
    family = meta.get("family", assembly.sseqid)
    base_name = meta.get("name", family)
    complete = assembly.left_end_covered and assembly.right_end_covered
    exact = (
        complete
        and assembly.structural_status == "intact"
        and assembly.identity >= 99.8
        and assembly.reference_coverage >= 99.5
    )
    named_variant = (
        complete
        and assembly.identity >= 95.0
        and assembly.reference_coverage >= 80.0
    )
    if exact:
        name = base_name
        variant_status = "exact_reference"
    elif named_variant:
        name = f"{base_name}-like"
        variant_status = (
            "minor_sequence_variant"
            if assembly.structural_status == "intact"
            else assembly.structural_status
        )
    else:
        name = f"{base_name} fragment"
        variant_status = assembly.structural_status

    attrs: dict[str, object] = {
        "seqid": assembly.qseqid,
        "reference_id": assembly.sseqid,
        "source": "reference_scan",
        "blast_identity": round(assembly.identity, 2),
        "blast_coverage": round(assembly.reference_coverage, 2),
        "blast_subject_coverage": round(assembly.reference_coverage, 2),
        "blast_hsp_count": assembly.hsp_count,
        "reference_length": assembly.slen,
        "reference_segments": _reference_segments(assembly),
        "variant_status": variant_status,
        "structural_status": assembly.structural_status,
        "query_span": assembly.query_span,
        "inserted_bases": assembly.inserted_bases,
        "deleted_bases": assembly.deleted_bases,
        "left_end_covered": assembly.left_end_covered,
        "right_end_covered": assembly.right_end_covered,
    }
    if not named_variant and not exact:
        attrs["fragment"] = True
    for key, value in meta.items():
        if key.startswith("_") or key in {"element_type", "family", "name"}:
            continue
        attrs.setdefault(key, value)

    return MGEFeature(
        element_type=meta.get("element_type", "transposon"),
        family=family,
        name=name,
        start=assembly.qstart,
        end=assembly.qend,
        strand=assembly.strand,
        score=assembly.identity,
        tsd_length=(
            int(meta["tsd_length"])
            if meta.get("tsd_length", "").isdigit()
            else TSD_LENGTHS.get(family)
        ),
        attributes=attrs,
    )


def _scan_collinear_unit_hits(
    hits: list[BlastHit],
    meta: dict[str, dict[str, str]],
) -> list[MGEFeature]:
    """Assemble interrupted/minor variants and choose one reference per locus."""
    chains = _chain_tn123_subject_hits(hits)
    assemblies = [_assemble_tn123_chain(chain) for chain in chains]
    features: list[MGEFeature] = []
    for group in _group_tn123_assemblies(assemblies):
        best = max(
            group,
            key=lambda item: (
                item.identity,
                item.reference_coverage,
                item.left_end_covered and item.right_end_covered,
                item.slen,
                -item.hsp_count,
            ),
        )
        features.append(_generic_assembly_to_feature(
            best,
            meta.get(best.sseqid, {"_id": best.sseqid}),
        ))
    # Sally's annotation order gives an exact complete unit priority over a
    # competing terminal fragment from a related reference. This matters for
    # Tn1721, whose partial duplication can otherwise produce an extra Tn1722
    # fragment at the end of the same complete locus.
    complete = [
        feature
        for feature in features
        if not feature.attributes.get("fragment")
        and feature.attributes.get("left_end_covered")
        and feature.attributes.get("right_end_covered")
    ]
    return [
        feature
        for feature in features
        if not (
            feature.attributes.get("fragment")
            and any(
                parent.attributes.get("seqid") == feature.attributes.get("seqid")
                and parent.start <= feature.start
                and feature.end <= parent.end
                for parent in complete
            )
        )
    ]


def scan(
    query_fasta: Path | str,
    reference_fasta: Path | str,
    min_identity: float = 85.0,
    min_length: int = 500,
    min_subject_coverage: float = 0.0,
    pick_best_variant: bool = False,
    prefer_identity: bool = False,
    exclude_ids: tuple[str, ...] = (),
    assemble_collinear: bool = False,
    evalue: float = 1e-10,
    blast_threads: int = 1,
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
        evalue=evalue,
        num_threads=blast_threads,
    )
    if exclude_ids:
        excluded = set(exclude_ids)
        hits = [hit for hit in hits if hit.sseqid not in excluded]

    if reference_fasta.name == "tn1_tn2_tn3.fasta":
        return _scan_tn123_hits(hits, meta)

    if assemble_collinear:
        return _scan_collinear_unit_hits(hits, meta)

    if pick_best_variant:
        features: list[MGEFeature] = []
        for group in _group_by_query_region(hits):
            best = _pick_best_variant(group, prefer_identity=prefer_identity)
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
    "mara_isecp1_tpu.fasta":   {"min_identity": 98.0, "min_length": 1_000,
                                 "min_subject_coverage": 70.0,
                                 "pick_best_variant": True,
                                 "prefer_identity": True},
    "mara_exact_is.fasta":      {"min_identity": 98.0, "min_length": 500},
    "mara_integron_segments.fasta": {"min_identity": 98.0, "min_length": 200},
    "mara_integron_ends.fasta": {"min_identity": 100.0, "min_length": 25,
                                  "evalue": 10.0},
    "mara_integron_tni.fasta": {"min_identity": 98.0, "min_length": 500,
                                 "assemble_collinear": True},
    "mara_partridge_units.fasta": {
        "min_identity": 95.0,
        "min_length": 500,
        "assemble_collinear": True,
    },
    "tn7.fasta":                {"min_identity": 90.0, "min_length": 3_000},
    "tn552.fasta":              {"min_identity": 95.0, "min_length": 2_000},
    # Canonical references selected by Sally Partridge. A lower length
    # threshold permits fragment reporting, but only near-complete hits are
    # assigned an exact Tn name.
    "tn1_tn2_tn3.fasta": {
        "min_identity": 95.0,
        "min_length": 400,
        "pick_best_variant": True,
        "prefer_identity": True,
    },
    # Tn3-family extras — the canonical Tn1/Tn2/Tn3 records above override
    # the older copies in this legacy collection.
    "tn3_family_extras.fasta": {
        "min_identity": 90.0,
        "min_length": 3_000,
        "min_subject_coverage": 60.0,
        "exclude_ids": ("Tn1", "Tn2", "Tn3"),
    },
    # PlasmidFinder replicon sequences are ~500-900bp per record.
    # Identity ≥95% is the standard PlasmidFinder cutoff.
    "plasmidfinder_enterobacteriales.fasta": {
        "min_identity": 95.0, "min_length": 200,
        "min_subject_coverage": 60.0,
    },
    # ISAba125 (IS30 family) — ISEScan fails to detect this element.
    # ~1087bp; high identity needed to avoid spurious IS30 hits.
    "isaba125.fasta":               {"min_identity": 90.0, "min_length": 800},
    # IS26 (IS6 family) — ISEScan fails to detect IS26 at sequence boundaries
    # (positions 1 or seq_length) due to insufficient flanking context.
    # ~820bp; high identity to avoid spurious IS6 family hits.
    "is26.fasta":                   {"min_identity": 95.0, "min_length": 600},
    # ter-site motifs for IS91 / ISCR rolling-circle elements.
    # Short consensus (~30bp) — needs a broad min_identity tolerance.
    # Marked as experimental pending a validated motif library.
    "rolling_circle_ter_sites.fasta": {
        "min_identity": 85.0, "min_length": 20,
    },
}


# The public alpha defaults to references validated against Sally Partridge's
# supplied examples. Broader historical/experimental references remain
# available only through the explicit ``all`` profile.
VALIDATED_REFERENCE_FILES = frozenset({
    "integron_archetypes.fasta",
    "isecp1.fasta",
    "mara_exact_is.fasta",
    "mara_integron_ends.fasta",
    "mara_integron_segments.fasta",
    "mara_integron_tni.fasta",
    "mara_isecp1_tpu.fasta",
    "mara_partridge_units.fasta",
    "tn21.fasta",
    "tn1_tn2_tn3.fasta",
})

REFERENCE_PROFILES = {
    "validated": VALIDATED_REFERENCE_FILES,
    "all": None,
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
    *,
    profile: str = "all",
    threads: int = 1,
    progress: Callable[[str, int, int], None] | None = None,
) -> list[MGEFeature]:
    """Scan every .fasta under references_dir with per-file parameters.

    Adds post-processing: if Tn6022+Tn6172 or Tn6019 hits are present,
    emit a parent AbGRI1 / AbaR3 genomic_island feature.
    """
    if profile not in REFERENCE_PROFILES:
        choices = ", ".join(sorted(REFERENCE_PROFILES))
        raise ValueError(f"unknown reference profile {profile!r}; choose {choices}")
    selected = REFERENCE_PROFILES[profile]
    refs = [
        ref
        for ref in sorted(references_dir.glob("*.fasta"))
        if selected is None or ref.name in selected
    ]
    worker_count = max(1, min(threads, len(refs) or 1))

    def _run(ref: Path) -> list[MGEFeature]:
        params = dict(
            REFERENCE_PARAMS.get(ref.name, {"min_identity": 85.0, "min_length": 500})
        )
        params["blast_threads"] = 1
        return scan(query_fasta, ref, **params)

    by_name: dict[str, list[MGEFeature]] = {}
    completed = 0
    if worker_count == 1:
        for ref in refs:
            by_name[ref.name] = _run(ref)
            completed += 1
            if progress:
                progress(ref.name, completed, len(refs))
    else:
        with ThreadPoolExecutor(max_workers=worker_count) as executor:
            futures = {executor.submit(_run, ref): ref for ref in refs}
            for future in as_completed(futures):
                ref = futures[future]
                by_name[ref.name] = future.result()
                completed += 1
                if progress:
                    progress(ref.name, completed, len(refs))

    out: list[MGEFeature] = []
    for ref in refs:
        out.extend(by_name[ref.name])
    # Compound-island detection (runs only over the hit set we just made)
    out.extend(_infer_acinetobacter_islands(out))
    return out
