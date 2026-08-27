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
from .element_definitions import tn123_definition, tn123_rules
from .tn123 import REFERENCE_METADATA, Tn123ComponentReference, component_references

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
    mismatch: int = 0
    gapopen: int = 0
    gaps: int = 0

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
    aligned_reference_bases: int
    mismatch_bases: int
    alignment_gap_bases: int
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
            "evalue qlen slen mismatch gapopen gaps"
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
        if len(cols) < 14:
            continue
        h = BlastHit(
            qseqid=cols[0], sseqid=cols[1],
            pident=float(cols[2]), length=int(cols[3]),
            qstart=int(cols[4]), qend=int(cols[5]),
            sstart=int(cols[6]), send=int(cols[7]),
            evalue=float(cols[8]),
            qlen=int(cols[9]), slen=int(cols[10]),
            mismatch=int(cols[11]), gapopen=int(cols[12]), gaps=int(cols[13]),
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
    inserted = sum(
        max(0, abs(hit.qend - hit.qstart) - abs(hit.send - hit.sstart))
        for hit in hits
    )
    deleted = sum(
        max(0, abs(hit.send - hit.sstart) - abs(hit.qend - hit.qstart))
        for hit in hits
    )
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
    mismatch_bases = sum(hit.mismatch for hit in hits)
    alignment_gap_bases = sum(hit.gaps for hit in hits)
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
        if inserted_bases and deleted_bases:
            status = "complete_with_indel"
        elif inserted_bases:
            status = "complete_with_insertion"
        elif deleted_bases:
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
        aligned_reference_bases=reference_bases,
        mismatch_bases=mismatch_bases,
        alignment_gap_bases=alignment_gap_bases,
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
            "mismatch_bases": hit.mismatch,
            "gap_bases": hit.gaps,
        }
        for hit in assembly.hits
    ]


def _tn123_assembly_to_feature(
    assembly: Tn123Assembly,
    runner_up: Tn123Assembly | None,
    meta: dict[str, str],
) -> MGEFeature:
    type_name = meta.get("type", meta.get("family", "Tn3_family"))
    definition_name = meta.get("name", type_name)
    rules = tn123_rules()["classification"]
    exact_rule = rules["exact_definition_match"]
    reference_rule = rules["reference_comparison"]
    close_rule = rules["close_variant"]
    margin = (
        assembly.identity - runner_up.identity
        if runner_up is not None
        else assembly.identity
    )
    complete = assembly.left_end_covered and assembly.right_end_covered
    clear_nearest = margin >= float(
        reference_rule["minimum_reference_margin_percent"]
    )
    exact = (
        complete
        and assembly.identity >= float(exact_rule["identity_percent"])
        and assembly.reference_coverage >= float(
            exact_rule["reference_coverage_percent"]
        )
        and assembly.mismatch_bases == 0
        and assembly.inserted_bases == 0
        and assembly.deleted_bases == 0
    )
    close_reference = (
        clear_nearest
        and complete
        and assembly.identity >= float(reference_rule["minimum_identity_percent"])
        and assembly.reference_coverage >= float(
            reference_rule["minimum_reference_coverage_percent"]
        )
    )

    if exact:
        comparison_status = "exact_reference"
    elif close_reference:
        comparison_status = (
            "minor_sequence_variant"
            if assembly.identity >= float(close_rule["minimum_identity_percent"])
            else assembly.structural_status
        )
    else:
        comparison_status = "ambiguous_reference"

    fragment = not complete or assembly.reference_coverage < float(
        reference_rule["minimum_reference_coverage_percent"]
    )
    name = (
        "Tn1/Tn2/Tn3 reference-comparison candidate"
        if complete
        else str(rules["fragment"]["visible_label"])
    )
    definition = tn123_definition(assembly.sseqid)
    attrs: dict[str, object] = {
        "seqid": assembly.qseqid,
        "reference_id": assembly.sseqid,
        "source": "reference_scan",
        "source_accession": meta.get("source_accession", ""),
        "definition_id": assembly.sseqid,
        "definition_version": meta.get("definition_version", ""),
        "defined_type": type_name,
        "defined_subtype": meta.get("subtype", ""),
        "closest_definition": definition_name,
        "closest_reference": assembly.sseqid,
        "reference_comparison_type": type_name,
        "reference_comparison_label": definition_name,
        "reference_comparison_status": comparison_status,
        "reference_comparison_only": True,
        "classification_basis": "reference_comparison_only",
        "reference_kind": meta.get("reference_kind", ""),
        "expert_rule": definition.get("expert_rule", ""),
        "known_differences_from_parent": definition.get(
            "known_differences_from_parent", []
        ),
        "review_status": definition.get("review_status", ""),
        "tn123_canonical": "true",
        "best_match": type_name,
        "blast_identity": round(assembly.identity, 2),
        "blast_coverage": round(assembly.reference_coverage, 2),
        "blast_subject_coverage": round(assembly.reference_coverage, 2),
        "blast_hsp_count": assembly.hsp_count,
        "reference_length": assembly.slen,
        "reference_segments": _reference_segments(assembly),
        "variant_margin": round(margin, 2),
        "variant_status": comparison_status,
        "structural_status": assembly.structural_status,
        "query_span": assembly.query_span,
        "reference_inserted_bases": assembly.inserted_bases,
        "reference_deleted_bases": assembly.deleted_bases,
        "reference_mismatch_bases": assembly.mismatch_bases,
        "inserted_bases": assembly.inserted_bases,
        "deleted_bases": assembly.deleted_bases,
        "mismatch_bases": assembly.mismatch_bases,
        "alignment_gap_bases": assembly.alignment_gap_bases,
        "left_end_covered": assembly.left_end_covered,
        "right_end_covered": assembly.right_end_covered,
    }
    if fragment:
        attrs["fragment"] = True
    if not clear_nearest:
        attrs["note"] = (
            "secondary Tn1, Tn2 and Tn3 reference scores are too close; "
            "component rules must classify the locus"
        )
    elif not exact:
        attrs["note"] = (
            f"secondary closest reference is {definition_name} ({type_name}); "
            "this comparison does not create the biological call"
        )

    return MGEFeature(
        element_type="transposon",
        family="Tn3_family",
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
    whole_locus_rules = tn123_rules()["component_detection"]["whole_locus"]
    chains = _chain_tn123_subject_hits(
        hits,
        max_query_insertion=int(whole_locus_rules["maximum_query_insertion_bp"]),
        max_reference_gap=int(whole_locus_rules["maximum_reference_gap_bp"]),
        max_reference_overlap=int(whole_locus_rules["maximum_reference_overlap_bp"]),
    )
    minimum_bases = int(
        whole_locus_rules["minimum_assembled_reference_bases"]
    )
    minimum_type_coverage = float(
        tn123_rules()["classification"]["reference_comparison"][
            "minimum_reference_coverage_percent"
        ]
    )
    assemblies = [
        assembly
        for chain in chains
        if (assembly := _assemble_tn123_chain(chain)).aligned_reference_bases
        >= minimum_bases
        # Long subtype references contain inserted material that can occur
        # elsewhere independently (for example ISEcp1-associated sequence in
        # Tn2.1). They are therefore eligible only as substantially complete
        # subtype matches. Partial-family evidence is always assessed against
        # the three canonical references, preventing inserted subtype context
        # from creating false Tn1/2/3 fragments.
        and (
            meta.get(assembly.sseqid, {}).get("reference_kind") == "canonical"
            or (
                assembly.left_end_covered
                and assembly.right_end_covered
                and assembly.reference_coverage >= minimum_type_coverage
            )
        )
    ]
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
        best_type = meta.get(best.sseqid, {}).get("type")
        runner_up = next(
            (
                item
                for item in ranked[1:]
                if meta.get(item.sseqid, {}).get("type") != best_type
            ),
            None,
        )
        features.append(_tn123_assembly_to_feature(
            best,
            runner_up,
            meta.get(best.sseqid, {"_id": best.sseqid}),
        ))
    ranked_features = sorted(
        features,
        key=lambda feature: (
            bool(not feature.attributes.get("fragment")),
            float(feature.attributes.get("blast_coverage", 0) or 0),
            feature.attributes.get("reference_kind") == "canonical",
            float(feature.attributes.get("blast_identity", 0) or 0),
            feature.end - feature.start,
        ),
        reverse=True,
    )
    kept: list[MGEFeature] = []
    for feature in ranked_features:
        redundant = False
        for selected in kept:
            if feature.attributes.get("seqid") != selected.attributes.get("seqid"):
                continue
            overlap = max(
                0,
                min(feature.end, selected.end) - max(feature.start, selected.start) + 1,
            )
            shorter = min(
                feature.end - feature.start + 1,
                selected.end - selected.start + 1,
            )
            if shorter and overlap / shorter >= 0.8:
                redundant = True
                break
        if not redundant:
            kept.append(feature)
    return sorted(kept, key=lambda feature: (feature.start, feature.end, feature.name))


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


def _component_strand(reference_strand: str, alignment_strand: str) -> str:
    if reference_strand == ".":
        return "."
    if alignment_strand == "+":
        return reference_strand
    return "+" if reference_strand == "-" else "-"


def _component_feature(
    assembly: Tn123Assembly,
    reference: Tn123ComponentReference,
) -> MGEFeature:
    complete = assembly.left_end_covered and assembly.right_end_covered
    attributes: dict[str, object] = {
        "seqid": assembly.qseqid,
        "source": "tn123_component_scan",
        "evidence_class": "sequence_detected",
        "component_role": reference.role,
        "component_reference": reference.reference_id,
        "reference_parent": reference.parent_reference,
        "reference_family": reference.parent_family,
        "blast_identity": round(assembly.identity, 2),
        "blast_coverage": round(assembly.reference_coverage, 2),
        "reference_length": assembly.slen,
        "reference_segments": _reference_segments(assembly),
        "component_status": "complete" if complete else "partial",
        "structural_status": assembly.structural_status,
    }
    if assembly.inserted_bases:
        attributes["interrupted"] = True
        attributes["inserted_bases"] = assembly.inserted_bases
    if assembly.deleted_bases:
        attributes["deleted_bases"] = assembly.deleted_bases
    if assembly.mismatch_bases:
        attributes["mismatch_bases"] = assembly.mismatch_bases
    if not complete:
        attributes["fragment"] = True
    name = "Tn3-family IR" if reference.role == "terminal_IR" else reference.name
    return MGEFeature(
        element_type=reference.element_type,
        family="Tn3_family",
        name=name,
        start=assembly.qstart,
        end=assembly.qend,
        strand=_component_strand(reference.strand, assembly.strand),
        score=assembly.identity,
        attributes=attributes,
    )


def _same_component_locus(left: MGEFeature, right: MGEFeature) -> bool:
    if left.attributes.get("seqid") != right.attributes.get("seqid"):
        return False
    if left.attributes.get("component_role") != right.attributes.get("component_role"):
        return False
    overlap = max(0, min(left.end, right.end) - max(left.start, right.start) + 1)
    shorter = min(left.end - left.start + 1, right.end - right.start + 1)
    return bool(shorter and overlap / shorter >= 0.7)


def _dedupe_component_features(features: list[MGEFeature]) -> list[MGEFeature]:
    groups: list[list[MGEFeature]] = []
    for feature in sorted(features, key=lambda item: (item.start, item.end, item.name)):
        for group in groups:
            if _same_component_locus(feature, group[0]):
                group.append(feature)
                break
        else:
            groups.append([feature])

    def rank(feature: MGEFeature) -> tuple[float, float, int]:
        return (
            float(feature.attributes.get("blast_identity", 0)),
            float(feature.attributes.get("blast_coverage", 0)),
            feature.end - feature.start,
        )

    selected: list[MGEFeature] = []
    for group in groups:
        best = max(group, key=rank)
        matches_by_type: dict[str, dict[str, object]] = {}
        for alternative in group:
            type_name = str(alternative.attributes.get("reference_family", ""))
            if type_name not in {"Tn1", "Tn2", "Tn3"}:
                continue
            identity = float(alternative.attributes.get("blast_identity", 0))
            coverage = float(alternative.attributes.get("blast_coverage", 0))
            profile_score = identity * coverage / 100.0
            current = matches_by_type.get(type_name)
            if current is not None and float(current["profile_score"]) >= profile_score:
                continue
            matches_by_type[type_name] = {
                "type": type_name,
                "component_reference": alternative.attributes.get(
                    "component_reference", ""
                ),
                "identity": round(identity, 3),
                "coverage": round(coverage, 3),
                "profile_score": round(profile_score, 3),
            }
        best.attributes["component_profile_matches"] = [
            matches_by_type[type_name]
            for type_name in sorted(matches_by_type)
        ]
        selected.append(best)
    return selected


def scan_tn123_components(
    query_fasta: Path | str,
    *,
    blast_threads: int = 1,
) -> list[MGEFeature]:
    """Detect canonical Tn1/Tn2/Tn3 components independently genome-wide."""
    canonical = REFERENCES_DIR / "tn1_tn2_tn3.fasta"
    references: list[Tn123ComponentReference] = []
    for record in SeqIO.parse(canonical, "fasta"):
        references.extend(component_references(record.id, str(record.seq)))
    by_id = {reference.reference_id: reference for reference in references}

    categories = {
        "long": [reference for reference in references if len(reference.sequence) >= 200],
        "res": [reference for reference in references if reference.role == "res"],
        "ir": [reference for reference in references if reference.role == "terminal_IR"],
    }
    detection_rules = tn123_rules()["component_detection"]
    parameters = {
        # Keep partial HSPs here so interrupted components can be reassembled
        # before the component-level coverage threshold is applied.
        "long": dict(
            min_identity=float(detection_rules["long"]["minimum_identity_percent"]),
            min_length=int(detection_rules["long"]["minimum_hsp_length"]),
            min_subject_coverage=0.0,
        ),
        "res": dict(
            min_identity=float(detection_rules["res"]["minimum_identity_percent"]),
            min_length=int(detection_rules["res"]["minimum_hsp_length"]),
            min_subject_coverage=float(
                detection_rules["res"]["minimum_subject_coverage_percent"]
            ),
        ),
        "ir": dict(
            min_identity=float(
                detection_rules["terminal_IR"]["minimum_identity_percent"]
            ),
            min_length=int(detection_rules["terminal_IR"]["minimum_hsp_length"]),
            min_subject_coverage=float(
                detection_rules["terminal_IR"]["minimum_subject_coverage_percent"]
            ),
            evalue=float(detection_rules["terminal_IR"]["evalue"]),
        ),
    }
    hits: list[BlastHit] = []
    with tempfile.TemporaryDirectory() as tmp:
        for category, category_references in categories.items():
            reference_path = Path(tmp) / f"tn123_{category}.fasta"
            reference_path.write_text("".join(
                f">{reference.reference_id}\n{reference.sequence}\n"
                for reference in category_references
            ))
            category_parameters = parameters[category]
            hits.extend(blastn_hits(
                query_fasta,
                reference_path,
                min_identity=float(category_parameters["min_identity"]),
                min_length=int(category_parameters["min_length"]),
                min_subject_coverage=float(
                    category_parameters["min_subject_coverage"]
                ),
                evalue=float(category_parameters.get("evalue", 1e-10)),
                num_threads=blast_threads,
            ))

    short_hits = [
        hit
        for hit in hits
        if by_id[hit.sseqid].role in {"terminal_IR", "res"}
    ]
    long_hits = [
        hit
        for hit in hits
        if by_id[hit.sseqid].role not in {"terminal_IR", "res"}
    ]
    chains = [[hit] for hit in short_hits]
    chains.extend(_chain_tn123_subject_hits(
        long_hits,
        max_reference_gap=int(detection_rules["long"]["maximum_reference_gap_bp"]),
    ))

    features: list[MGEFeature] = []
    for chain in chains:
        assembly = _assemble_tn123_chain(chain)
        reference = by_id[assembly.sseqid]
        # Short components must be near-complete; longer coding components may
        # remain useful as explicit partial or interrupted evidence.
        minimum_coverage = (
            float(detection_rules["res"]["minimum_assembled_coverage_percent"])
            if len(reference.sequence) < int(
                detection_rules["long"]["minimum_reference_length"]
            )
            else float(detection_rules["long"]["minimum_assembled_coverage_percent"])
        )
        if assembly.reference_coverage < minimum_coverage:
            continue
        features.append(_component_feature(assembly, reference))
    return _dedupe_component_features(features)


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
        "min_identity": float(
            tn123_rules()["component_detection"]["whole_locus"][
                "minimum_identity_percent"
            ]
        ),
        "min_length": int(
            tn123_rules()["component_detection"]["whole_locus"][
                "minimum_hsp_length"
            ]
        ),
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
    # Expert-rule discovery without complete Tn1/Tn2/Tn3 reference lookup.
    # The component models still come from reviewed canonical examples.
    "tn123-components": (),
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
    if (
        profile == "tn123-components"
        or selected is None
        or "tn1_tn2_tn3.fasta" in selected
    ):
        out.extend(scan_tn123_components(query_fasta, blast_threads=1))
    # Compound-island detection (runs only over the hit set we just made)
    out.extend(_infer_acinetobacter_islands(out))
    return out
