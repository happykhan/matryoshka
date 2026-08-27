"""Contract tests for the source-backed locus component ontology."""

import json

from click.testing import CliRunner

from matryoshka.__main__ import cli
from matryoshka.component_catalog import REQUIRED_COMPONENT_FIELDS, load_component_catalog
from matryoshka.locus_map import SUPPORTED_RENDER_SYMBOLS


def test_catalog_has_complete_component_and_assembly_contracts():
    catalog = load_component_catalog()
    assert len(catalog["component_types"]) == 30
    assert len(catalog["assembly_grammars"]) == 18
    for component in catalog["component_types"].values():
        assert component.keys() >= REQUIRED_COMPONENT_FIELDS
        assert component["render"]["symbol"] in SUPPORTED_RENDER_SYMBOLS
        assert component["current_status"] in {
            "detected",
            "partially_detected",
            "reference_needed",
        }


def test_catalog_preserves_reviewed_exact_names_and_corrections():
    catalog = load_component_catalog()
    examples = catalog["named_examples"]
    assert {"ISKpn6", "ISKpn7", "ISKpn31", "IS1999", "IS26", "ISEcp1"} <= set(
        examples["insertion_sequences"]
    )
    assert {"Tn1", "Tn2", "Tn3", "Tn7", "Tn21", "Tn1721", "Tn1722", "Tn4401"} <= set(
        examples["transposons"]
    )
    assert {"blaOXA-48", "blaOXA-181", "blaCTX-M-15", "blaCMY-2"} <= set(
        examples["cassette_and_resistance_examples"]
    )
    is26_exclusions = " ".join(
        catalog["assembly_grammars"]["is26_pseudocomposite"]["excludes"]
    )
    assert "oppositely oriented IS26" in is26_exclusions
    assert "generic IS6" in is26_exclusions
    tpu_exclusions = " ".join(catalog["assembly_grammars"]["isecp1_tpu"]["excludes"])
    assert "adjacent Tn1721" in tpu_exclusions


def test_catalog_cli_exports_tsv_and_complete_json(tmp_path):
    runner = CliRunner()
    tsv = runner.invoke(cli, ["catalog", "--format", "tsv"])
    assert tsv.exit_code == 0
    assert tsv.output.startswith("component\tlevel\tdefinition")
    assert "resistance_gene\tatomic\t" in tsv.output

    target = tmp_path / "catalog.json"
    result = runner.invoke(
        cli,
        ["catalog", "--format", "json", "--out", str(target)],
    )
    assert result.exit_code == 0
    exported = json.loads(target.read_text())
    assert exported["assembly_grammars"]["gene_cassette"]["requires"][-1] == "cognate attC"
