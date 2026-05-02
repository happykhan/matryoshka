"""
transposon.py — Infer composite transposons from IS element + cargo patterns.

Three inference patterns are supported, each corresponding to a distinct
biological mobilisation mechanism (Partridge et al., 2018):

1. FLANKED-CARGO composites (``FLANKED_RULES``)
       [left IS] -- cargo -- [right IS]
   Classical composite transposon / Tn3-family unit. Two complete IS
   copies flanking one or more AMR cargo genes.
   Examples: Tn4401 (ISKpn7/blaKPC/ISKpn6), Tn1999 (IS1999/blaOXA-48/IS1999),
   mcr-1 composite (ISApl1/mcr-1/ISApl1).

2. ONE-ENDED capture (``ONE_ENDED_RULES``)
       [IS] -- cargo
   A single complete IS upstream of cargo, using one IR plus a cryptic
   end-beyond-IRR to mobilise adjacent DNA. Characteristic of IS1380
   family (ISEcp1 → blaCTX-M / blaCMY). No flanking TSD.

3. ROLLING-CIRCLE capture (``ROLLING_CIRCLE_RULES``)
       [IS] -- cargo   (no TSD; ter-site mediated)
   Single IS91-family or ISCR element adjacent to cargo. Replication
   beyond the *ter* site carries adjacent DNA. No TSD expected.

4. SIGNATURE detection (``infer_tn1546_by_vana``)
       Pattern-matched cargo clusters that imply a known unit transposon
       even without flanking IS evidence (e.g. Tn1546 = vanA + vanH + vanX).

IS26 resistance islands are handled separately because they don't follow
the single-cargo flanked pattern: IS26 clusters recursively via the
"translocatable unit" mechanism, producing dense mosaics. We detect all
IS26/IS26 pairs with >=1 AMR cargo between them, then merge overlapping
intervals into a single island.

Only *complete* IS copies (ISEScan type == "c") participate in inference.
Partial copies (ISEScan "p" / name prefix "new_") are ignored — a
fragmentary IS cannot bracket a functional transposon.
"""

from __future__ import annotations

from dataclasses import dataclass

from .detect import MGEFeature

# ---------------------------------------------------------------------------
# Rule definitions
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class FlankedRule:
    """Flanked-cargo composite transposon rule."""
    family: str
    name: str
    cargo_match: str          # case-insensitive substring in AMR name
    left_family: str
    right_family: str
    upstream_max: int
    downstream_max: int


@dataclass(frozen=True)
class OneEndedRule:
    """Single-IS capture (one-ended transposition) rule."""
    family: str               # output family, e.g. "ISEcp1_capture"
    name: str
    is_family: str            # required IS family (complete copies only)
    cargo_matches: tuple[str, ...]   # any of these substrings (case-insensitive)
    max_gap: int              # max bp between IS end and cargo start


@dataclass(frozen=True)
class SignatureRule:
    """Cargo-signature detection (no IS required).

    Matches when the `required_genes` are all present within `max_span` bp
    of one another. Emits one feature spanning the full set of matched
    cargo positions.
    """
    family: str
    name: str
    required_genes: tuple[str, ...]   # all these substrings must appear
    max_span: int


# ---------------------------------------------------------------------------
# Rule tables
# ---------------------------------------------------------------------------

FLANKED_RULES: list[FlankedRule] = [
    FlankedRule(
        family="Tn4401", name="Tn4401",
        cargo_match="KPC",
        left_family="IS21", right_family="IS1182",
        upstream_max=500, downstream_max=1000,
    ),
    FlankedRule(
        family="Tn1999", name="Tn1999",
        cargo_match="OXA-48",
        left_family="IS4", right_family="IS4",
        upstream_max=2000, downstream_max=200,
    ),
    # mcr-1 colistin resistance composite — ISApl1 (IS30 family) pair.
    # Canonical Tn number: Tn6330. Partridge 2018: CP016184 (single IS),
    # KY689635 (both IS).
    FlankedRule(
        family="Tn6330", name="Tn6330",
        cargo_match="mcr-1",
        left_family="IS30", right_family="IS30",
        upstream_max=500, downstream_max=500,
    ),
    # Tn2006 — ISAba1 flanking blaOXA-23 (A. baumannii carbapenem resistance).
    # ISAba1 is IS4-family; we match on the family (ISEScan classification).
    # downstream_max=2000: element contains helicase/methyltransferase between
    # OXA-23 and the downstream ISAba1, so the gap can reach ~1600bp.
    FlankedRule(
        family="Tn2006", name="Tn2006",
        cargo_match="OXA-23",
        left_family="IS4", right_family="IS4",
        upstream_max=500, downstream_max=2000,
    ),
    # Tn125 — ISAba125 flanking blaNDM (A. baumannii NDM metallo-beta-lactamase).
    # ISAba125 is IS30-family. The canonical Tn125 is ~10kb; NDM sits
    # near the left ISAba125 (gap ~94bp), while the right ISAba125 is
    # ~7kb downstream because the element contains ble, trpF, dsbC,
    # cutA, groES, groEL, ISCR21, and aphA6 between NDM and the right IS.
    FlankedRule(
        family="Tn125", name="Tn125",
        cargo_match="NDM",
        left_family="IS30", right_family="IS30",
        upstream_max=1000, downstream_max=8000,
    ),
    # Tn10 — IS10R/IS10L (IS4-family) flanking tet(B) tetracycline resistance.
    # The canonical Tn10 is 9.3kb; IS10 elements are ~1.3kb each.
    # upstream_max=5000: tet(B) sits ~4 kb from the left IS10 because
    # the element contains jemC, orf1, orf2, tetR, and tetC between the
    # left IS10 and tet(B). downstream_max=1500: additional ORFs between
    # tet(B) and the right IS10.
    FlankedRule(
        family="Tn10", name="Tn10",
        cargo_match="tet",
        left_family="IS4", right_family="IS4",
        upstream_max=5000, downstream_max=1500,
    ),
    # Tn5 — IS50L/IS50R (IS4-family) flanking aminoglycoside resistance cluster
    # (aph(3')-IIa, aph(6)-Ic, ble). Canonical Tn5 is 5.8kb.
    # The cargo cluster (aph(3')-IIa, ble, aph(6)-Ic) spans ~2kb between
    # the two IS50 copies. infer_flanked works per-cargo-gene, so thresholds
    # must accommodate the gap from any single aph gene to both IS elements:
    # upstream_max=1500: aph(6)-Ic sits ~1.2kb from IS50L
    # downstream_max=2000: aph(3')-IIa sits ~1.9kb from IS50R
    FlankedRule(
        family="Tn5", name="Tn5",
        cargo_match="aph",
        left_family="IS4", right_family="IS4",
        upstream_max=1500, downstream_max=2000,
    ),
]


ONE_ENDED_RULES: list[OneEndedRule] = [
    # ISEcp1 captures blaCTX-M, blaCMY, blaACC etc. via IRL + cryptic end.
    # Partridge 2018 accession: FJ621588.
    # Patterns are matched case-insensitively as substrings; they are
    # chosen to be specific enough to avoid hitting blaTEM / blaSHV.
    OneEndedRule(
        family="ISEcp1_capture", name="ISEcp1_capture",
        is_family="IS1380",
        cargo_matches=("CTX-M", "CMY-", "blaACC", "blaLAT", "blaDHA", "qnrB", "qnrA", "qnrS"),
        # ISEcp1 one-ended transposition can capture adjacent DNA up to
        # several kb beyond the IRR. 3 kb is a conservative practical limit.
        max_gap=3000,
    ),
]


ROLLING_CIRCLE_RULES: list[OneEndedRule] = [
    # ISCR1 / ISCR2 rolling-circle capture units (IS91-like).
    # Partridge 2018: carry blaCTX-M-9 variants, qnr genes, sul2, etc.
    OneEndedRule(
        family="ISCR_capture", name="ISCR_capture",
        is_family="ISCR",
        cargo_matches=("CTX-M", "qnr", "sul", "mph", "arm"),
        max_gap=800,
    ),
    OneEndedRule(
        family="IS91_capture", name="IS91_capture",
        is_family="IS91",
        cargo_matches=("CTX-M", "CMY", "NDM"),
        max_gap=800,
    ),
]


SIGNATURE_RULES: list[SignatureRule] = [
    # Tn1546 — vanA vancomycin-resistance cluster.
    # Partridge 2018 accession: M97297.
    SignatureRule(
        family="Tn1546", name="Tn1546",
        required_genes=("vanA",),
        max_span=12_000,
    ),
    # Tn1331 — aac(6')-Ib + blaOXA-9 + aadA1 clustered within a Tn3-family unit.
    # Partridge 2018 accession: AF479774.
    SignatureRule(
        family="Tn1331", name="Tn1331",
        required_genes=("aac(6')-Ib", "blaOXA-9", "aadA1"),
        max_span=8_000,
    ),
]


# ---------------------------------------------------------------------------
# IS26 island parameters
# ---------------------------------------------------------------------------

IS26_MAX_PAIR_SPAN = 20_000


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _is_complete(f: MGEFeature) -> bool:
    """True if the IS element is a complete (not partial) copy."""
    return f.attributes.get("type", "c") == "c"


def _complete_is_of(features: list[MGEFeature], family: str) -> list[MGEFeature]:
    # family may be exact or a prefix (e.g. "ISCR" matches "ISCR1", "ISCR2")
    return [
        f for f in features
        if f.element_type == "IS"
        and (f.family == family or f.family.startswith(family))
        and _is_complete(f)
    ]


def _cargo_matches(f: MGEFeature, patterns: tuple[str, ...] | str) -> bool:
    if f.element_type != "AMR":
        return False
    name = f.name.upper()
    if isinstance(patterns, str):
        patterns = (patterns,)
    return any(p.upper() in name for p in patterns)


# ---------------------------------------------------------------------------
# Flanked-cargo inference
# ---------------------------------------------------------------------------

def infer_flanked(rule: FlankedRule, features: list[MGEFeature]) -> list[MGEFeature]:
    """Infer transposons matching a flanked-cargo rule."""
    left_is = _complete_is_of(features, rule.left_family)
    right_is = _complete_is_of(features, rule.right_family)
    cargo = [f for f in features if _cargo_matches(f, rule.cargo_match)]

    results: list[MGEFeature] = []
    for c in cargo:
        upstream = [f for f in left_is
                    if f.end <= c.start and c.start - f.end <= rule.upstream_max]
        if not upstream:
            continue
        up = max(upstream, key=lambda f: f.end)

        downstream = [f for f in right_is
                      if f.start >= c.end and f.start - c.end <= rule.downstream_max]
        if not downstream:
            continue
        down = min(downstream, key=lambda f: f.start)

        results.append(MGEFeature(
            element_type="transposon",
            family=rule.family,
            name=rule.name,
            start=up.start,
            end=down.end,
            strand=c.strand,
        ))
    return results


# ---------------------------------------------------------------------------
# One-ended / rolling-circle inference
# ---------------------------------------------------------------------------

def infer_one_ended(rule: OneEndedRule, features: list[MGEFeature]) -> list[MGEFeature]:
    """Infer a single-IS-plus-adjacent-cargo unit.

    IS may sit on either side of the cargo; the emitted transposon spans
    min(IS.start, cargo.start) .. max(IS.end, cargo.end).
    """
    iss = _complete_is_of(features, rule.is_family)
    cargo = [f for f in features if _cargo_matches(f, rule.cargo_matches)]

    emitted: list[MGEFeature] = []
    used_is: set[int] = set()
    for c in cargo:
        candidates: list[tuple[int, MGEFeature]] = []
        for i in iss:
            if id(i) in used_is:
                continue
            if i.end <= c.start and c.start - i.end <= rule.max_gap:
                candidates.append((c.start - i.end, i))
            elif i.start >= c.end and i.start - c.end <= rule.max_gap:
                candidates.append((i.start - c.end, i))
        if not candidates:
            continue
        candidates.sort(key=lambda t: t[0])
        _, best_is = candidates[0]
        used_is.add(id(best_is))

        span_start = min(best_is.start, c.start)
        span_end = max(best_is.end, c.end)
        emitted.append(MGEFeature(
            element_type="transposon",
            family=rule.family,
            name=rule.name,
            start=span_start,
            end=span_end,
            strand=c.strand,
            attributes={"cargo": c.name, "mediator_is": best_is.name},
        ))
    return emitted


# ---------------------------------------------------------------------------
# Signature-based inference
# ---------------------------------------------------------------------------

def infer_signature(
    rule: SignatureRule,
    features: list[MGEFeature],
) -> list[MGEFeature]:
    """Detect cargo-signature transposons (e.g. Tn1546 via vanA cluster).

    All `required_genes` must appear in a single cluster spanning no more
    than `max_span` bp. Emits one transposon per distinct cluster.
    """
    matches_by_gene: dict[str, list[MGEFeature]] = {
        g: [f for f in features
            if f.element_type == "AMR" and g.lower() in f.name.lower()]
        for g in rule.required_genes
    }
    if any(not v for v in matches_by_gene.values()):
        return []

    anchors = matches_by_gene[rule.required_genes[0]]
    results: list[MGEFeature] = []
    seen_spans: set[tuple[int, int]] = set()
    for anchor in anchors:
        # Find siblings within the span window on either side
        cluster: list[MGEFeature] = [anchor]
        for g in rule.required_genes[1:]:
            nearest = min(
                (m for m in matches_by_gene[g]
                 if abs(m.start - anchor.start) <= rule.max_span),
                key=lambda m: abs(m.start - anchor.start),
                default=None,
            )
            if nearest is None:
                cluster = []
                break
            cluster.append(nearest)
        if not cluster:
            continue
        start = min(m.start for m in cluster)
        end = max(m.end for m in cluster)
        if (start, end) in seen_spans:
            continue
        seen_spans.add((start, end))
        results.append(MGEFeature(
            element_type="transposon",
            family=rule.family,
            name=rule.name,
            start=start,
            end=end,
            strand=".",
            attributes={
                "signature_genes": ",".join(m.name for m in cluster),
            },
        ))
    return results


# ---------------------------------------------------------------------------
# Legacy per-family shims
# ---------------------------------------------------------------------------

def infer_tn4401(features: list[MGEFeature]) -> list[MGEFeature]:
    return infer_flanked(FLANKED_RULES[0], features)


def infer_tn1999(features: list[MGEFeature]) -> list[MGEFeature]:
    return infer_flanked(FLANKED_RULES[1], features)


# ---------------------------------------------------------------------------
# IS26 island inference (unchanged)
# ---------------------------------------------------------------------------

def _merge_overlapping_islands(
    islands: list[MGEFeature],
    amr: list[MGEFeature],
) -> list[MGEFeature]:
    if not islands:
        return []
    merged = sorted(islands, key=lambda f: (f.start, f.end))
    changed = True
    while changed:
        changed = False
        out: list[MGEFeature] = [merged[0]]
        for isl in merged[1:]:
            prev = out[-1]
            if isl.start <= prev.end:
                new_start = min(prev.start, isl.start)
                new_end = max(prev.end, isl.end)
                out[-1] = MGEFeature(
                    element_type="transposon",
                    family="IS26_island",
                    name=prev.name,
                    start=new_start,
                    end=new_end,
                    strand=".",
                )
                changed = True
            else:
                out.append(isl)
        merged = out

    final: list[MGEFeature] = []
    for i, isl in enumerate(merged, start=1):
        cargo = [a for a in amr if a.start >= isl.start and a.end <= isl.end]
        isl.name = f"IS26_island_{i}"
        isl.attributes["cargo_count"] = len(cargo)
        isl.attributes["cargo"] = ",".join(a.name for a in cargo)
        final.append(isl)
    return final


def infer_is26_islands(features: list[MGEFeature]) -> list[MGEFeature]:
    is26 = sorted(_complete_is_of(features, "IS6"), key=lambda f: f.start)
    amr = [f for f in features if f.element_type == "AMR"]

    candidates: list[MGEFeature] = []
    for i, left in enumerate(is26):
        for right in is26[i + 1:]:
            if right.end - left.start > IS26_MAX_PAIR_SPAN:
                break
            cargo = [a for a in amr
                     if a.start >= left.end and a.end <= right.start]
            if not cargo:
                continue
            candidates.append(MGEFeature(
                element_type="transposon",
                family="IS26_island",
                name="IS26_island",
                start=left.start,
                end=right.end,
                strand=".",
            ))
    return _merge_overlapping_islands(candidates, amr)


# Max distance for single-IS26 translocatable unit detection
IS26_TU_MAX_GAP = 1_000


def infer_is26_translocatable_units(features: list[MGEFeature]) -> list[MGEFeature]:
    """Detect single-IS26 translocatable units.

    Partridge 2018 emphasises that the *unit of mobility* for IS26 is one
    IS26 copy plus an adjacent region — the "translocatable unit" (TU) —
    which inserts preferentially next to an existing IS26 (~50× random).
    This rule emits a low-confidence TU feature when a complete IS26 is
    within IS26_TU_MAX_GAP bp of AMR cargo that is NOT already enclosed
    by a two-IS26 island.
    """
    is26 = _complete_is_of(features, "IS6")
    amr = [f for f in features if f.element_type == "AMR"]
    islands = [f for f in features
               if f.element_type == "transposon" and f.family == "IS26_island"]

    def _inside_island(s: int, e: int) -> bool:
        return any(isl.start <= s and e <= isl.end for isl in islands)

    emitted: list[MGEFeature] = []
    used_pairs: set[tuple[int, int]] = set()
    for i in is26:
        for a in amr:
            if _inside_island(a.start, a.end):
                continue
            if i.end <= a.start and a.start - i.end <= IS26_TU_MAX_GAP:
                span = (i.start, a.end)
            elif a.end <= i.start and i.start - a.end <= IS26_TU_MAX_GAP:
                span = (a.start, i.end)
            else:
                continue
            key = (span[0], span[1])
            if key in used_pairs:
                continue
            used_pairs.add(key)
            emitted.append(MGEFeature(
                element_type="transposon",
                family="IS26_TU",
                name=f"IS26_TU::{a.name}",
                start=span[0],
                end=span[1],
                strand=".",
                attributes={
                    "cargo": a.name,
                    "mediator_is": i.name,
                    "confidence": "low",
                    "note": "single-IS26 translocatable unit (Partridge 2018)",
                },
            ))
    return emitted


def infer_is26_composites(features: list[MGEFeature]) -> list[MGEFeature]:
    return infer_is26_islands(features)


# ---------------------------------------------------------------------------
# Top-level driver
# ---------------------------------------------------------------------------

# Tn3-superfamily *unit* transposons — have tnpA/tnpR/res between IRL and IRR.
# Composites (Tn1999, Tn_mcr1, IS26_island) are NOT included — they lack a
# Tn3-family res site because they don't encode their own resolvase.
TN3_FAMILY_MEMBERS = frozenset({
    "Tn3", "Tn4401", "Tn1546", "Tn21", "Tn1331", "Tn5393",
    "Tn6022", "Tn6019", "Tn6021", "Tn6172",
})


def annotate_res_sites(transposons: list[MGEFeature]) -> list[MGEFeature]:
    """Emit expected-res-site sub-features for Tn3-superfamily transposons.

    The Tn3-family res site sits between tnpA and tnpR, typically within
    the first ~1 kb of the element starting from the 3' end of tnpA. We
    can't sequence-confirm it without a gene-prediction step, but we can
    annotate the *expected* location so downstream viewers (Artemis,
    IGV, our own SVG) show it.

    Returns a list of new MGEFeature(element_type="res_site") that callers
    can append to the feature set before building the hierarchy.
    """
    out: list[MGEFeature] = []
    for t in transposons:
        if t.family not in TN3_FAMILY_MEMBERS:
            continue
        # Expected res region: ~200-400bp from IRL (start for + strand,
        # end for - strand). 120bp wide.
        if t.strand == "-":
            res_end = t.end - 200
            res_start = res_end - 120
        else:
            res_start = t.start + 200
            res_end = res_start + 120
        if res_start < 1 or res_end <= res_start:
            continue
        out.append(MGEFeature(
            element_type="res_site",
            family=t.family,
            name=f"res_{t.name}",
            start=res_start,
            end=res_end,
            strand=t.strand,
            attributes={
                "parent_transposon": t.name,
                "note": "expected_position_not_sequence_confirmed",
            },
        ))
    return out


def infer_transposons(features: list[MGEFeature]) -> list[MGEFeature]:
    """Run all transposon inference rules and return the new features.

    Rules are applied in priority order:
      1. Flanked-cargo composites (most specific structural evidence)
      2. One-ended captures (single IS + cargo)
      3. Rolling-circle captures (ISCR / IS91 + cargo)
      4. Signature clusters (cargo-only pattern matches)
      5. IS26 islands (merged pairs)

    De-duplication: emitted features are unique by (family, start, end).

    Note: res-site annotation is applied downstream (in the CLI) so it
    uses authoritative BLAST-confirmed transposon boundaries when they
    exist, rather than approximate rule-inferred ones.
    """
    out: list[MGEFeature] = []
    seen: set[tuple[str, int, int]] = set()

    def _extend(items: list[MGEFeature]) -> None:
        for f in items:
            key = (f.family, f.start, f.end)
            if key in seen:
                continue
            seen.add(key)
            out.append(f)

    for rule in FLANKED_RULES:
        _extend(infer_flanked(rule, features))
    for rule in ONE_ENDED_RULES:
        _extend(infer_one_ended(rule, features))
    for rule in ROLLING_CIRCLE_RULES:
        _extend(infer_one_ended(rule, features))
    for srule in SIGNATURE_RULES:
        _extend(infer_signature(srule, features))
    _extend(infer_is26_islands(features))
    # Single-IS26 translocatable units: low-confidence, emitted only for
    # cargo not already covered by an island
    _extend(infer_is26_translocatable_units(features + out))
    return out
