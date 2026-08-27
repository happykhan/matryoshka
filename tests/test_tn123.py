"""Tests for Sally Partridge's canonical Tn1/Tn2/Tn3 workflow."""

from __future__ import annotations

import json
import random
from pathlib import Path

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from click.testing import CliRunner

from matryoshka.__main__ import cli
from matryoshka.detect import MGEFeature
from matryoshka.hierarchy import build_hierarchy
from matryoshka.mara_table import to_mara_table_svg
from matryoshka.mara_viz import to_mara_svg
from matryoshka.reference_scan import (
    REFERENCES_DIR,
    BlastHit,
    _group_by_query_region,
    _hit_to_feature,
    blast_available,
    scan,
    scan_tn123_components,
)
from matryoshka.tn123 import (
    REFERENCE_METADATA,
    annotate_tn123,
    assemble_tn123_components,
    curated_internal_features,
)

PARTRIDGE_FASTA = (
    Path(__file__).parent / "test-data" / "partridge-examples" / "Tn1-Tn2-Tn3.fasta"
)
TN123_REFERENCE = REFERENCES_DIR / "tn1_tn2_tn3.fasta"
REFERENCE_SEQUENCES = {
    record.id: str(record.seq)
    for record in SeqIO.parse(TN123_REFERENCE, "fasta")
}


def _random_dna(length: int, seed: int) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(length))


def _write_fasta(path: Path, name: str, sequence: str) -> Path:
    path.write_text(f">{name}\n{sequence}\n")
    return path


def _mutate(sequence: str, count: int) -> str:
    bases = list(sequence)
    positions = random.Random(71).sample(range(50, len(bases) - 50), count)
    replacement = {"A": "C", "C": "G", "G": "T", "T": "A"}
    for position in positions:
        bases[position] = replacement[bases[position]]
    return "".join(bases)


def _parent(family: str = "Tn1", strand: str = "+") -> MGEFeature:
    length = {"Tn1": 4949, "Tn2": 4950, "Tn3": 4948}[family]
    return MGEFeature(
        element_type="transposon",
        family=family,
        name=family,
        start=101,
        end=100 + length,
        strand=strand,
        attributes={
            "seqid": "contig",
            "blast_identity": 100.0,
            "blast_subject_coverage": 100.0,
            "reference_length": length,
            "source_accession": REFERENCE_METADATA[
                {"Tn1": "Tn1_NC_008357", "Tn2": "Tn2_AY123253", "Tn3": "Tn3_HM749966"}[family]
            ]["source_accession"],
        },
    )


def test_canonical_reference_accessions():
    assert REFERENCE_METADATA["Tn1_NC_008357"]["source_accession"] == "NC_008357"
    assert REFERENCE_METADATA["Tn2_AY123253"]["source_accession"] == "AY123253"
    assert REFERENCE_METADATA["Tn3_HM749966"]["source_accession"] == "HM749966"
    assert (REFERENCES_DIR / "tn1_tn2_tn3.fasta").exists()


def test_curated_internal_layout():
    parent = _parent("Tn2")
    features = curated_internal_features(parent)
    assert [feature.name for feature in features] == [
        "IRL", "blaTEM-1b", "tnpR", "res", "tnpA", "IRR",
    ]
    assert features[0].start == parent.start
    assert features[-1].end == parent.end
    assert features[0].end - features[0].start + 1 == 38
    assert features[-1].end - features[-1].start + 1 == 38
    assert parent.attributes["curated_internal_features"] is True
    assert parent.tsd_length == 5


def test_reverse_strand_projection():
    parent = _parent("Tn3", "-")
    features = curated_internal_features(parent)
    bla = next(feature for feature in features if feature.name == "blaTEM-1")
    tnpa = next(feature for feature in features if feature.name == "tnpA")
    assert bla.strand == "+"
    assert tnpa.strand == "-"
    assert features[0].end == parent.end


def test_short_match_is_reported_as_ambiguous_fragment():
    hit = BlastHit("contig", "Tn1_NC_008357", 99.2, 800, 50, 849, 1, 800, 0.0, 2000, 4949)
    meta = dict(REFERENCE_METADATA["Tn1_NC_008357"])
    meta["blast_subject_coverage"] = "16.2"
    feature = _hit_to_feature(hit, meta)
    assert feature.name == "Tn1/2/3"
    assert feature.family == "Tn3_family"
    assert feature.attributes["fragment"] is True
    assert feature.attributes["best_match"] == "Tn1"


def test_variant_groups_do_not_mix_contigs():
    one = BlastHit("contig_a", "Tn1", 100.0, 1000, 1, 1000, 1, 1000, 0.0, 1000, 1000)
    two = BlastHit("contig_b", "Tn2", 100.0, 1000, 1, 1000, 1, 1000, 0.0, 1000, 1000)
    assert len(_group_by_query_region([one, two])) == 2


def test_mara_renderer_is_separate_single_line_output():
    parent = _parent("Tn1")
    parent.children = curated_internal_features(parent)
    svg = to_mara_svg([parent], parent.end, "Tn1 example")
    assert svg.startswith("<svg")
    assert "Tn1 example" in svg
    assert "blaTEM-2" in svg
    assert "tnpR" in svg
    assert "tnpA" in svg
    assert 'fill="#009b55"' in svg
    assert 'fill="#000"' in svg
    assert "expected 5 bp TSD; flanking sequence unavailable" in svg


def test_mara_table_contains_hierarchy_and_evidence_columns():
    parent = _parent("Tn2")
    parent.children = curated_internal_features(parent)
    svg = to_mara_table_svg([parent], "Tn2 example")
    assert svg.startswith("<svg")
    assert "Position" in svg
    assert "Name*" in svg
    assert "FID" in svg
    assert "Type" in svg
    assert "Notes" in svg
    assert "AY123253" in svg
    assert "blaTEM-1b" in svg
    assert "expected DR=5 bp" in svg


@pytest.mark.skipif(not blast_available(), reason="blastn not on PATH")
def test_components_are_independently_detected_across_canonical_records():
    components = scan_tn123_components(PARTRIDGE_FASTA)
    by_contig: dict[str, list[MGEFeature]] = {}
    for component in components:
        by_contig.setdefault(str(component.attributes["seqid"]), []).append(component)
    assert set(by_contig) == {
        "Tn1_NC_008357", "Tn2_AY123253", "Tn3_HM749966",
    }
    for detected in by_contig.values():
        roles = [str(feature.attributes["component_role"]) for feature in detected]
        assert roles.count("terminal_IR") == 2
        assert set(roles) == {"terminal_IR", "blaTEM", "tnpR", "res", "tnpA"}
        assert all(
            feature.attributes["evidence_class"] == "sequence_detected"
            for feature in detected
        )


@pytest.mark.skipif(not blast_available(), reason="blastn not on PATH")
def test_component_grammar_can_emit_parent_without_whole_locus_call():
    components = [
        feature
        for feature in scan_tn123_components(PARTRIDGE_FASTA)
        if feature.attributes["seqid"] == "Tn1_NC_008357"
    ]
    parents = assemble_tn123_components(components)
    assert len(parents) == 1
    assert parents[0].name == "Tn3-family unit"
    assert parents[0].attributes["source"] == "component_assembly"
    assert parents[0].attributes["component_assembly_status"] == "complete"
    assert parents[0].attributes["component_order_valid"] is True
    assert parents[0].attributes["detected_component_count"] == 6


@pytest.mark.skipif(not blast_available(), reason="blastn not on PATH")
def test_component_grammar_rejects_missing_required_res_site():
    components = [
        feature
        for feature in scan_tn123_components(PARTRIDGE_FASTA)
        if feature.attributes["seqid"] == "Tn1_NC_008357"
        and feature.attributes["component_role"] != "res"
    ]
    assert assemble_tn123_components(components) == []


@pytest.mark.skipif(not blast_available(), reason="blastn not on PATH")
def test_random_sequence_has_no_tn123_component_calls(tmp_path: Path):
    query = _write_fasta(
        tmp_path / "negative.fasta",
        "random_negative",
        _random_dna(50_000, 90210),
    )
    assert scan_tn123_components(query) == []


@pytest.mark.skipif(not blast_available(), reason="blastn not on PATH")
def test_tn1_is_detected_inside_arbitrary_long_sequence(tmp_path: Path):
    prefix = _random_dna(20_000, 1)
    suffix = _random_dna(20_000, 2)
    tn1 = REFERENCE_SEQUENCES["Tn1_NC_008357"]
    query = _write_fasta(tmp_path / "embedded.fasta", "arbitrary_contig", prefix + tn1 + suffix)
    hits = scan(query, TN123_REFERENCE, min_identity=95.0, min_length=400)
    assert len(hits) == 1
    hit = hits[0]
    assert (hit.name, hit.start, hit.end, hit.strand) == (
        "Tn1", len(prefix) + 1, len(prefix) + len(tn1), "+",
    )
    assert hit.attributes["variant_status"] == "exact_reference"


@pytest.mark.skipif(not blast_available(), reason="blastn not on PATH")
def test_minor_sequence_variation_is_named_as_nearest_variant(tmp_path: Path):
    tn1 = REFERENCE_SEQUENCES["Tn1_NC_008357"]
    variant = _mutate(tn1, 25)
    query = _write_fasta(tmp_path / "variant.fasta", "minor_variant", variant)
    hits = scan(query, TN123_REFERENCE, min_identity=95.0, min_length=400)
    assert len(hits) == 1
    hit = hits[0]
    assert hit.name == "Tn1-like"
    assert hit.attributes["best_match"] == "Tn1"
    assert hit.attributes["variant_status"] == "minor_sequence_variant"
    assert 99.0 < float(hit.attributes["blast_identity"]) < 100.0


@pytest.mark.skipif(not blast_available(), reason="blastn not on PATH")
def test_inserted_sequence_assembles_into_one_complete_locus(tmp_path: Path):
    tn2 = REFERENCE_SEQUENCES["Tn2_AY123253"]
    insertion = _random_dna(800, 3)
    interrupted = tn2[:2500] + insertion + tn2[2500:]
    query = _write_fasta(tmp_path / "inserted.fasta", "tn2_with_insertion", interrupted)
    hits = scan(query, TN123_REFERENCE, min_identity=95.0, min_length=400)
    assert len(hits) == 1
    hit = hits[0]
    assert hit.name == "Tn2-like"
    assert (hit.start, hit.end) == (1, len(interrupted))
    assert hit.attributes["structural_status"] == "complete_with_insertion"
    assert 790 <= int(hit.attributes["inserted_bases"]) <= 810
    assert int(hit.attributes["blast_hsp_count"]) >= 2

    roots = build_hierarchy(hits + annotate_tn123(hits))
    diagram = to_mara_svg(roots, len(interrupted), "interrupted Tn2")
    table = to_mara_table_svg(roots, "interrupted Tn2")
    assert "Tn2-like" in diagram
    assert "insertion ≈" in diagram
    assert "inserted sequence" in table


@pytest.mark.skipif(not blast_available(), reason="blastn not on PATH")
def test_internal_deletion_is_reported_as_one_variant_locus(tmp_path: Path):
    tn2 = REFERENCE_SEQUENCES["Tn2_AY123253"]
    deleted = tn2[:2200] + tn2[2500:]
    query = _write_fasta(tmp_path / "deleted.fasta", "tn2_with_deletion", deleted)
    hits = scan(query, TN123_REFERENCE, min_identity=95.0, min_length=400)
    assert len(hits) == 1
    hit = hits[0]
    assert hit.name == "Tn2-like"
    assert (hit.start, hit.end) == (1, len(deleted))
    assert hit.attributes["structural_status"] == "complete_with_deletion"
    assert 290 <= int(hit.attributes["deleted_bases"]) <= 310


@pytest.mark.skipif(not blast_available(), reason="blastn not on PATH")
def test_reverse_locus_with_insertion_is_assembled_and_oriented(tmp_path: Path):
    tn3 = REFERENCE_SEQUENCES["Tn3_HM749966"]
    insertion = _random_dna(600, 31)
    interrupted = tn3[:2600] + insertion + tn3[2600:]
    query = _write_fasta(
        tmp_path / "reverse_inserted.fasta",
        "reverse_inserted_tn3",
        str(Seq(interrupted).reverse_complement()),
    )
    hits = scan(query, TN123_REFERENCE, min_identity=95.0, min_length=400)
    assert len(hits) == 1
    hit = hits[0]
    assert hit.name == "Tn3-like"
    assert hit.strand == "-"
    assert hit.attributes["structural_status"] == "complete_with_insertion"
    assert 590 <= int(hit.attributes["inserted_bases"]) <= 610


@pytest.mark.skipif(not blast_available(), reason="blastn not on PATH")
def test_reverse_complement_is_detected_on_reverse_strand(tmp_path: Path):
    tn3 = REFERENCE_SEQUENCES["Tn3_HM749966"]
    query = _write_fasta(
        tmp_path / "reverse.fasta",
        "reverse_tn3",
        _random_dna(1000, 4) + str(Seq(tn3).reverse_complement()) + _random_dna(1000, 5),
    )
    hits = scan(query, TN123_REFERENCE, min_identity=95.0, min_length=400)
    assert len(hits) == 1
    assert hits[0].name == "Tn3"
    assert hits[0].strand == "-"


@pytest.mark.skipif(not blast_available(), reason="blastn not on PATH")
def test_partial_locus_keeps_covered_end_and_fragment_status(tmp_path: Path):
    tn1 = REFERENCE_SEQUENCES["Tn1_NC_008357"]
    partial = tn1[700:]
    query = _write_fasta(tmp_path / "partial.fasta", "partial_tn1", partial)
    hits = scan(query, TN123_REFERENCE, min_identity=95.0, min_length=400)
    assert len(hits) == 1
    hit = hits[0]
    assert hit.name == "Tn1-like"
    assert hit.attributes["fragment"] is True
    assert hit.attributes["structural_status"] == "left_partial"
    assert hit.attributes["left_end_covered"] is False
    assert hit.attributes["right_end_covered"] is True


@pytest.mark.skipif(not blast_available(), reason="blastn not on PATH")
def test_two_loci_on_one_contig_get_separate_mara_outputs(tmp_path: Path):
    tn1 = REFERENCE_SEQUENCES["Tn1_NC_008357"]
    tn3 = REFERENCE_SEQUENCES["Tn3_HM749966"]
    sequence = (
        _random_dna(3000, 6) + tn1 + _random_dna(4000, 7)
        + tn3 + _random_dna(3000, 8)
    )
    query = _write_fasta(tmp_path / "two.fasta", "two_loci", sequence)
    out_dir = tmp_path / "mara-loci"
    result = CliRunner().invoke(
        cli,
        ["annotate", str(query), "--format", "mara", "--mara-flank", "1000",
         "--out", str(out_dir)],
    )
    assert result.exit_code == 0, result.output
    outputs = sorted(out_dir.glob("*.svg"))
    assert len(outputs) == 2
    assert any("Tn1" in output.name for output in outputs)
    assert any("Tn3" in output.name for output in outputs)


@pytest.mark.skipif(not blast_available(), reason="blastn not on PATH")
def test_reference_only_cli_detects_all_three_and_writes_mara(tmp_path: Path):
    json_out = tmp_path / "tn123.json"
    runner = CliRunner()
    result = runner.invoke(
        cli,
        ["annotate", str(PARTRIDGE_FASTA), "--format", "json", "-o", str(json_out)],
    )
    assert result.exit_code == 0, result.output
    detected = json.loads(json_out.read_text())
    assert {
        contig: [(feature["name"], feature["attributes"]["blast_identity"])
                 for feature in features]
        for contig, features in detected.items()
    } == {
        "Tn1_NC_008357": [("Tn1", 100.0)],
        "Tn2_AY123253": [("Tn2", 100.0)],
        "Tn3_HM749966": [("Tn3", 100.0)],
    }
    for features in detected.values():
        parent = features[0]
        assert parent["attributes"]["component_assembly_status"] == "complete"
        assert parent["attributes"]["component_order_valid"] is True
        assert parent["attributes"]["detected_component_count"] == 6
        assert parent["attributes"]["reference_projected_component_count"] == 0
        assert sum(
            child["element_type"] == "res_site"
            for child in parent["children"]
        ) == 1
        assert {
            child["attributes"]["evidence_class"] for child in parent["children"]
            if child["element_type"] in {"IR", "AMR", "gene", "res_site"}
        } == {"sequence_detected"}

    mara_dir = tmp_path / "mara"
    result = runner.invoke(
        cli,
        ["annotate", str(PARTRIDGE_FASTA), "--format", "mara", "-o", str(mara_dir)],
    )
    assert result.exit_code == 0, result.output
    svgs = sorted(mara_dir.glob("*.svg"))
    assert len(svgs) == 3
    assert all('fill="#009b55"' in svg.read_text() for svg in svgs)

    table_dir = tmp_path / "mara-table"
    result = runner.invoke(
        cli,
        ["annotate", str(PARTRIDGE_FASTA), "--format", "mara-table", "-o", str(table_dir)],
    )
    assert result.exit_code == 0, result.output
    tables = sorted(table_dir.glob("*.svg"))
    assert len(tables) == 3
    assert all("Position" in table.read_text() for table in tables)
