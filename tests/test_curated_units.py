"""Regression tests for expert-reviewed unit transposons outside Tn1/Tn2/Tn3."""

from __future__ import annotations

from pathlib import Path

import pytest
from Bio import SeqIO

from matryoshka.__main__ import _suppress_redundant_inference
from matryoshka.curated_units import UNIT_MAPS, annotate_curated_units
from matryoshka.detect import MGEFeature
from matryoshka.hierarchy import build_hierarchy
from matryoshka.locus_map import to_locus_map_svg
from matryoshka.locus_table import to_locus_table_svg
from matryoshka.reference_scan import REFERENCES_DIR, blast_available, scan

REVIEWED_FASTA = (
    Path(__file__).parent / "test-data" / "reviewed-examples" / "other_Tn.fasta"
)
UNIT_REFERENCE = REFERENCES_DIR / "reviewed_units.fasta"
UNIT_IDS = {
    "Tn1721_AB366441",
    "Tn1722_AB366441",
    "Tn4401_GU595196",
    "Tn5393_AF262622",
    "Tn5403_X75779",
}


def _parent(
    family: str,
    reference_id: str,
    length: int,
    strand: str = "+",
) -> MGEFeature:
    return MGEFeature(
        element_type="transposon",
        family=family,
        name=family,
        start=101,
        end=100 + length,
        strand=strand,
        attributes={
            "reference_id": reference_id,
            "blast_identity": 100.0,
            "blast_coverage": 100.0,
            "source": "reference_scan",
        },
    )


def test_every_supplied_unit_has_a_curated_internal_map():
    assert UNIT_MAPS.keys() >= UNIT_IDS


def test_in2_reference_is_the_exact_iri_to_irt_region_not_a_whole_accession():
    records = list(SeqIO.parse(REFERENCES_DIR / "integron_archetypes.fasta", "fasta"))
    assert len(records) == 1
    assert records[0].id == "In2_AF071413_4040_15039"
    assert len(records[0].seq) == 11_000
    assert "coordinates=4040..15039" in records[0].description


def test_tn4401_map_uses_exact_is_names_and_reviewed_boundaries():
    parent = _parent("Tn4401", "Tn4401_GU595196", 9_907)
    children = annotate_curated_units([parent])
    by_name = {child.name: child for child in children}
    assert (by_name["ISKpn7"].start, by_name["ISKpn7"].end) == (5_028, 6_983)
    assert by_name["ISKpn7"].family == "IS21"
    assert (by_name["ISKpn6"].start, by_name["ISKpn6"].end) == (8_166, 9_705)
    assert by_name["ISKpn6"].family == "IS1182"
    assert by_name["ISKpn6"].strand == "-"
    assert "blaKPC region" in by_name
    assert parent.attributes["curated_internal_features"] is True


def test_reverse_tn1721_projection_preserves_internal_ir_copy():
    parent = _parent("Tn1721", "Tn1721_AB366441", 11_128, strand="-")
    children = annotate_curated_units([parent])
    irrs = [child for child in children if child.name == "IRR"]
    assert len(irrs) == 2
    assert sorted((child.start, child.end) for child in irrs) == [
        (101, 138),
        (5_589, 5_626),
    ]


def test_tn21_map_contains_integron_ends_segments_cassette_and_tni_fragment():
    parent = _parent("Tn21", "Tn21", 19_672)
    children = annotate_curated_units([parent])
    names = [child.name for child in children]
    assert names == [
        "IRR", "tnp", "IRi", "5'-CS", "aadA1 cassette region",
        "3'-CS", "tni", "IRt", "mer", "IRL",
    ]
    tni = next(child for child in children if child.name == "tni")
    assert tni.attributes["fragment"] is True
    roots = build_hierarchy([parent, *children])
    svg = to_locus_map_svg(roots, parent.end, "Tn21 locus map")
    table = to_locus_table_svg(roots, "Tn21 locus map")
    assert "5&#x27;-CS" in svg
    assert "aadA1 cassette region" in svg
    assert "tni#" in svg
    assert 'fill="#9f1d2d"' in svg
    assert 'fill="#9b4b9d"' in svg
    assert ">cassette<" in table
    assert "LOCUS-REF-5CS-GGG" in table
    assert "LOCUS-REF-Tn402-tni" in table


def test_reviewed_reference_overrides_legacy_duplicate_at_the_same_locus():
    reviewed = _parent("Tn1721", "Tn1721_AB366441", 11_128)
    reviewed.attributes["provenance"] = "expert_reviewed"
    legacy = _parent("Tn3", "Tn1721", 11_128)
    legacy.name = "Tn1721"
    kept = _suppress_redundant_inference([legacy, reviewed])
    assert kept == [reviewed]


def test_curated_components_override_duplicate_reference_scans_without_disappearing():
    curated_is = MGEFeature(
        "IS", "IS1182", "ISKpn6", 8066, 9605, "-",
        attributes={"source": "curated_reference", "fid": "LOCUS-REF-ISKpn6"},
    )
    scanned_is = MGEFeature(
        "IS", "IS1182", "ISKpn6", 8066, 9605, "-",
        attributes={"source": "reference_scan", "fid": "LOCUS-REF-ISKpn6"},
    )
    curated_tni = MGEFeature(
        "transposon_component", "Tn402", "tni", 12358, 15039, ".",
        attributes={"source": "curated_reference", "fid": "LOCUS-REF-Tn402-tni"},
    )
    scanned_tni = MGEFeature(
        "transposon_component", "Tn402", "tni fragment", 12362, 15039, "+",
        attributes={"source": "reference_scan", "fid": "LOCUS-REF-Tn402-tni"},
    )
    kept = _suppress_redundant_inference(
        [scanned_is, curated_is, scanned_tni, curated_tni]
    )
    assert kept == [curated_is, curated_tni]


@pytest.mark.skipif(not blast_available(), reason="blastn not on PATH")
def test_all_five_supplied_units_are_exact_reference_calls(tmp_path: Path):
    records = [
        record
        for record in SeqIO.parse(REVIEWED_FASTA, "fasta")
        if record.id in UNIT_IDS
    ]
    query = tmp_path / "reviewed-units.fasta"
    SeqIO.write(records, query, "fasta")
    hits = scan(
        query,
        UNIT_REFERENCE,
        min_identity=95.0,
        min_length=500,
        assemble_collinear=True,
    )
    assert len(hits) == 5
    assert {feature.attributes["seqid"] for feature in hits} == UNIT_IDS
    assert {feature.family for feature in hits} == {
        "Tn1721", "Tn1722", "Tn4401", "Tn5393", "Tn5403",
    }
    assert all(feature.attributes["blast_coverage"] == 100.0 for feature in hits)
    children = annotate_curated_units(hits)
    assert {child.attributes["reference_parent"] for child in children} == UNIT_IDS


@pytest.mark.skipif(not blast_available(), reason="blastn not on PATH")
def test_internal_insertion_is_one_tn5393_like_locus_with_projected_terminal_ir(
    tmp_path: Path,
):
    source = next(
        record
        for record in SeqIO.parse(REVIEWED_FASTA, "fasta")
        if record.id == "Tn5393_AF262622"
    )
    insertion = "ACGT" * 150
    sequence = str(source.seq[:2500]) + insertion + str(source.seq[2500:])
    query = tmp_path / "inserted-tn5393.fasta"
    query.write_text(f">inserted_tn5393\n{sequence}\n")
    hits = scan(
        query,
        UNIT_REFERENCE,
        min_identity=95.0,
        min_length=500,
        assemble_collinear=True,
    )
    assert len(hits) == 1
    assert hits[0].name == "Tn5393-like"
    assert hits[0].attributes["structural_status"] == "complete_with_insertion"
    assert 590 <= int(hits[0].attributes["inserted_bases"]) <= 610
    children = annotate_curated_units(hits)
    irr = next(child for child in children if child.name == "IRR")
    assert irr.end == len(sequence)
    insertion_feature = next(
        child for child in children if child.element_type == "inserted_sequence"
    )
    assert insertion_feature.end - insertion_feature.start + 1 == 600
