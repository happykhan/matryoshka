"""Tests for component-aware locus class-1-integron reconstruction."""

from matryoshka.detect import MGEFeature
from matryoshka.integron_structures import infer_integron_structures


def part(element_type, family, name, start, end, strand="+", **attributes):
    return MGEFeature(
        element_type=element_type,
        family=family,
        name=name,
        start=start,
        end=end,
        strand=strand,
        attributes={"seqid": "contig", **attributes},
    )


def test_complete_tn402_requires_both_ends_segments_and_complete_tni():
    features = [
        part("IR", "Tn402", "IRi", 100, 124),
        part("integron_segment", "class1_integron", "5'-CS", 100, 1454),
        part("integron_segment", "class1_integron", "3'-CS", 2500, 4738),
        part("transposon_component", "Tn402", "tni", 5000, 9732),
        part("IR", "Tn402", "IRt", 9708, 9732, strand="-"),
    ]
    inferred = infer_integron_structures(features)
    integron = next(feature for feature in inferred if feature.element_type == "integron")
    tn402 = next(feature for feature in inferred if feature.family == "Tn402")
    assert (integron.start, integron.end) == (100, 4738)
    assert integron.attributes["iri_status"] == "confirmed"
    assert integron.attributes["irt_status"] == "confirmed"
    assert (tn402.start, tn402.end) == (100, 9732)
    assert tn402.attributes["structural_status"] == "complete_component_set"


def test_partial_tni_does_not_create_complete_tn402():
    features = [
        part("IR", "Tn402", "IRi", 100, 124),
        part("integron_segment", "class1_integron", "5'-CS", 100, 1454),
        part("integron_segment", "class1_integron", "3'-CS", 2500, 4738),
        part("transposon_component", "Tn402", "tni", 5000, 7000, fragment=True),
        part("IR", "Tn402", "IRt", 6976, 7000, strand="-"),
    ]
    inferred = infer_integron_structures(features)
    assert any(feature.element_type == "integron" for feature in inferred)
    assert not any(feature.family == "Tn402" for feature in inferred)


def test_existing_reference_integron_prevents_duplicate_candidate():
    features = [
        part("integron", "class1_integron", "In2", 50, 5000),
        part("integron_segment", "class1_integron", "5'-CS", 100, 1454),
        part("integron_segment", "class1_integron", "3'-CS", 2500, 4738),
    ]
    assert infer_integron_structures(features) == []
