import json

from matryoshka.detect import MGEFeature
from matryoshka.report import annotation_gff3, annotation_json


def _feature(name: str = "Tn test") -> MGEFeature:
    child = MGEFeature("AMR", "beta lactam", "blaTEST", 20, 80, "+")
    return MGEFeature(
        "transposon",
        "Tn3",
        name,
        1,
        100,
        "+",
        attributes={"note": "minor variant; review=true"},
        children=[child],
    )


def test_json_envelope_and_identifiers_are_stable_and_contig_scoped(tmp_path):
    fasta = tmp_path / "input.fasta"
    fasta.write_text(">contig A\n" + "A" * 100 + "\n>contig-B\n" + "C" * 100 + "\n")
    records = [("contig", 100, [_feature()]), ("contig-B", 100, [_feature()])]

    first = json.loads(annotation_json(records, input_path=fasta, profile="validated"))
    second = json.loads(annotation_json(records, input_path=fasta, profile="validated"))

    assert first["schema_version"] == "1.0"
    assert first["reference_database"]["version"]
    first_ids = [record["features"][0]["id"] for record in first["records"]]
    second_ids = [record["features"][0]["id"] for record in second["records"]]
    assert first_ids == second_ids
    assert len(set(first_ids)) == 2
    assert first["records"][0]["features"][0]["children"][0]["parent_id"] == first_ids[0]


def test_gff3_is_one_document_with_encoded_attributes_and_unique_ids():
    gff = annotation_gff3([
        ("contig-A", 100, [_feature()]),
        ("contig-B", 100, [_feature()]),
    ])
    assert gff.count("##gff-version 3") == 1
    assert gff.count("##sequence-region") == 2
    assert "minor%20variant%3B%20review%3Dtrue" in gff
    rows = [line.split("\t") for line in gff.splitlines() if not line.startswith("#")]
    ids = [field.split("=", 1)[1] for row in rows for field in row[8].split(";") if field.startswith("ID=")]
    assert len(ids) == len(set(ids))
    assert all(len(row) == 9 for row in rows)
