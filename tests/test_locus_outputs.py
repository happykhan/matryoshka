"""Regression tests for the published locus feature vocabulary."""

from pathlib import Path

import pytest
from Bio import SeqIO

from matryoshka.__main__ import _annotate_contig
from matryoshka.detect import MGEFeature, parse_amrfinder, parse_integron_finder
from matryoshka.locus_map import to_locus_map_svg
from matryoshka.locus_table import to_locus_table_svg
from matryoshka.locus_views import extract_locus_views
from matryoshka.reference_scan import (
    REFERENCES_DIR,
    BlastHit,
    _hit_to_feature,
    blast_available,
    scan,
)


def feature(
    element_type: str,
    name: str,
    start: int,
    end: int,
    strand: str = "+",
    family: str = "",
    **attributes: object,
) -> MGEFeature:
    return MGEFeature(
        element_type=element_type,
        family=family,
        name=name,
        start=start,
        end=end,
        strand=strand,
        attributes=attributes,
    )


def test_renderer_uses_locus_symbols_for_core_feature_vocabulary():
    tn21 = feature("transposon", "Tn21", 100, 6000, family="Tn21")
    tn21.tsd_length = 5
    tn21.tsd_seq = "AATGC"
    tn21.children = [
        feature("IS", "IS1R", 250, 1000, family="IS1"),
        feature("integron_segment", "5'-CS", 1200, 1800),
        feature("cassette", "aadA1a", 1801, 2600),
        feature("attC", "attC_aadA1a", 2550, 2600, strand="."),
        feature("integron_segment", "3'-CS", 2601, 3500),
        feature("group_II_intron", "Kl.pn.I2", 3800, 4600, strand="-"),
    ]
    svg = to_locus_map_svg([tn21], 6500, "locus vocabulary")
    assert 'fill="#9f1d2d"' in svg
    assert 'fill="#ffffff"' in svg
    assert 'fill="#ed873b"' in svg
    assert 'fill="#68c5e3"' in svg
    assert 'fill="#b8b8b8"' in svg
    assert "IS1R" in svg
    assert "5&#x27;-CS" in svg
    assert "aadA1a" in svg
    assert "Kl.pn.I2" in svg
    assert "boundary-adjacent sequence-matched TSD: AATGC" in svg
    assert "gene (arrow shows strand)" in svg
    assert "terminal inverted repeat (IR)" in svg
    assert "sequence-matched DR/TSD" in svg
    assert "expected, unconfirmed DR/TSD" in svg
    assert "sequence-detected versus reference-projected evidence" in svg


def test_renderer_draws_the_extended_component_ontology():
    region = feature("multiresistance_region", "MRR", 100, 9800)
    integron = feature("integron", "class 1 integron", 500, 4000)
    array = feature("cassette_array", "cassette array", 1150, 2900)
    array.children = [
        feature("cassette_remnant", "dfrA17 remnant", 1400, 2100),
        feature("attC", "attC_dfrA17", 2050, 2100, strand="."),
    ]
    integron.children = [
        feature("attI", "attI1", 980, 1010, strand="."),
        feature("Pc_promoter", "Pc", 1050, 1080),
        array,
    ]
    region.children = [
        integron,
        feature("ISCR", "ISCR27#", 4300, 5400),
        feature("oriIS", "oriIS", 4310, 4330, strand="."),
        feature("terIS", "terIS", 5380, 5400, strand="."),
        feature("captured_segment", "captured bla region", 5401, 6500),
        feature("replicon", "IncFII", 6900, 7600),
        feature("oriV", "oriV", 7650, 7680, strand="."),
        feature("oriT", "oriT", 7900, 7930, strand="."),
        feature("ncRNA", "RNAI", 8150, 8500, strand="-"),
        feature("direct_repeat", "DR-left", 8750, 8754, strand=".", sequence="TATGA"),
        feature("unknown_fragment", "unresolved feature", 9000, 9300),
    ]
    svg = to_locus_map_svg([region], 10_000, "extended locus vocabulary")
    for label in (
        "MRR",
        "attI1",
        "Pc",
        "dfrA17 remnant#",
        "ISCR27#",
        "oriIS",
        "terIS",
        "captured bla region",
        "IncFII",
        "oriV",
        "oriT",
        "RNAI",
        "unresolved feature#",
    ):
        assert label in svg
    assert "boundary-adjacent sequence-matched TSD: TATGA" in svg
    assert 'stroke-dasharray="6,3"' in svg


def test_table_preserves_locus_feature_types():
    features = [
        feature("integron_segment", "5'-CS", 100, 600),
        feature("cassette", "aadA1a", 601, 1400),
        feature("group_II_intron", "Kl.pn.I2", 1500, 2200),
    ]
    table = to_locus_table_svg(features, "locus table")
    assert ">CS<" in table
    assert ">cassette<" in table
    assert ">intron<" in table
    assert "terminal inverted repeat (IR)" in table
    assert "Solid components are sequence-detected" in table
    assert "dashed/outlined components are reference-projected" in table


def test_partial_exact_is_is_marked_but_not_renamed_to_family():
    partial = feature("IS", "IS26", 100, 500, family="IS6", type="p", fragment=True)
    svg = to_locus_map_svg([partial], 1000)
    assert "IS26#" in svg
    assert "IS6" not in svg


def test_reference_hit_orientation_accounts_for_reverse_subject_alignment():
    hit = BlastHit("contig", "IRt", 100.0, 25, 100, 124, 25, 1, 0.0, 500, 25)
    irt = _hit_to_feature(hit, {
        "element_type": "IR", "family": "Tn402", "name": "IRt",
    })
    assert irt.strand == "-"


def test_complete_isecp1_tpu_gets_a_readable_locus_and_table_type():
    tpu = feature(
        "transposition_unit",
        "ISEcp1-blaCMY-2_TPU",
        10_001,
        14_479,
        family="ISEcp1_TPU",
        structural_status="complete_reference_match",
    )
    tpu.tsd_length = 5
    loci = extract_locus_views([tpu], 50_000, flank=2_000)
    assert len(loci) == 1
    assert (loci[0].view_start, loci[0].view_end) == (8_001, 16_479)
    svg = to_locus_map_svg(loci[0].roots, loci[0].view_length, "ISEcp1 TPU")
    table = to_locus_table_svg(loci[0].roots, "ISEcp1 TPU")
    assert "ISEcp1-blaCMY-2_TPU" in svg
    assert "expected 5 bp TSD" in svg
    assert ">TPU<" in table


def test_broad_region_and_mobile_child_get_views_without_duplicate_integron_view():
    region = feature("multiresistance_region", "MRR", 10_000, 40_000)
    tn1 = feature("transposon", "Tn1", 15_000, 20_000, family="Tn1")
    tn1.children = [feature("integron", "nested integron", 16_000, 19_000)]
    region.children = [tn1]
    loci = extract_locus_views([region], 100_000, flank=1_000)
    assert [locus.target.name for locus in loci] == ["MRR", "Tn1"]


def test_bundled_isecp1_reference_is_the_exact_is_not_a_whole_accession():
    records = list(SeqIO.parse(REFERENCES_DIR / "isecp1.fasta", "fasta"))
    assert len(records) == 1
    assert records[0].id == "ISEcp1_FJ621588"
    assert len(records[0].seq) == 1_656
    assert "name=ISEcp1" in records[0].description


def test_integron_cassette_includes_attc_and_survives_cli_hierarchy(tmp_path: Path):
    integrons = tmp_path / "example.integrons"
    integrons.write_text(
        "ID_integron\tID_replicon\tpos_beg\tpos_end\tstrand\ttype\t"
        "type_elt\tannotation\telement\tevalue\n"
        "integron_1\tcontig\t100\t900\t1\tcomplete\tprotein\tintI\tintI1\t0\n"
        "integron_1\tcontig\t1000\t1800\t1\tcomplete\tprotein\tprotein\taadA1\t0\n"
        "integron_1\tcontig\t1760\t1900\t1\tcomplete\tattC\tattC\tattC_aadA1\t1e-20\n"
    )
    parsed = parse_integron_finder(integrons)
    cassette = next(child for child in parsed[0].children if child.element_type == "cassette")
    assert (cassette.start, cassette.end) == (1000, 1900)
    assert cassette.attributes["boundary_status"] == "complete_with_attC"
    roots, all_features = _annotate_contig("N" * 2_000, "contig", parsed, True)
    assert {feature.element_type for feature in all_features} >= {
        "integron", "integrase", "cassette", "attC",
    }
    root = next(feature for feature in roots if feature.element_type == "integron")
    cassette = next(child for child in root.children if child.element_type == "cassette")
    assert any(child.element_type == "attC" for child in cassette.children)


def test_generic_blaoxa_is_explicitly_unresolved_not_presented_as_an_allele(
    tmp_path: Path,
):
    report = tmp_path / "amrfinder.tsv"
    report.write_text(
        "Contig id\tStart\tStop\tStrand\tElement symbol\tElement name\tClass\t"
        "Subclass\tMethod\n"
        "contig\t100\t900\t+\tblaOXA\tOXA beta-lactamase\tBETA-LACTAM\t"
        "CARBAPENEM\tALLELE\n"
    )
    generic = parse_amrfinder(report)[0]
    assert generic.name == "blaOXA-family (unresolved)"
    assert generic.attributes["requires_exact_allele"] is True
    assert generic.attributes["confidence"] == "low"


@pytest.mark.skipif(not blast_available(), reason="blastn not on PATH")
def test_expert_supplied_exact_is_names_are_recovered():
    supplied = (
        Path(__file__).parent
        / "test-data"
        / "reviewed-examples"
        / "other_Tn.fasta"
    )
    hits = scan(
        supplied,
        REFERENCES_DIR / "locus_exact_is.fasta",
        min_identity=98.0,
        min_length=500,
    )
    assert {feature.name for feature in hits} == {"ISKpn6", "ISKpn7"}
    assert all(feature.attributes["type"] == "c" for feature in hits)


@pytest.mark.skipif(not blast_available(), reason="blastn not on PATH")
def test_expert_supplied_complete_isecp1_tpus_are_recovered():
    supplied = (
        Path(__file__).parent
        / "test-data"
        / "reviewed-examples"
        / "TPU_CMY-2-like.fasta"
    )
    hits = scan(
        supplied,
        REFERENCES_DIR / "locus_isecp1_tpu.fasta",
        min_identity=98.0,
        min_length=1_000,
        min_subject_coverage=70.0,
        pick_best_variant=True,
        prefer_identity=True,
    )
    assert len(hits) == 7
    assert {feature.attributes["seqid"] for feature in hits} == {
        "TPU_CMY-2_3078+503_FM246884",
        "TPU_CMY-2_3910_CP001121",
        "TPU_CMY-2_4025_constructed",
        "TPU_CMY-2_4087_AY509004",
        "TPU_CMY-2_4479_FJ621588",
        "TPU_CMY-31_4025_EU331425",
        "TPU_CMY-36_4025_EU331426",
    }
    assert all(feature.family == "ISEcp1_TPU" for feature in hits)
    assert all(
        feature.attributes["structural_status"] == "complete_reference_match"
        for feature in hits
    )
    assert all(feature.start == 1 for feature in hits)
    assert all(
        feature.end == int(feature.attributes["reference_length"])
        for feature in hits
    )
    assert all(feature.attributes["reference_id"] == feature.attributes["seqid"] for feature in hits)
