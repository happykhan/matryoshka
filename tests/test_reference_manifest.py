from pathlib import Path

from matryoshka.reference_manifest import load_reference_manifest, sha256_file
from matryoshka.reference_scan import REFERENCES_DIR, VALIDATED_REFERENCE_FILES


def test_manifest_covers_and_pins_every_bundled_fasta():
    manifest = load_reference_manifest()
    entries = {entry["file"]: entry for entry in manifest["references"]}
    fasta_names = {path.name for path in REFERENCES_DIR.glob("*.fasta")}
    assert manifest["schema_version"] == "1.0"
    assert manifest["reference_count"] == len(entries) == len(fasta_names)
    assert set(entries) == fasta_names
    for filename, entry in entries.items():
        assert entry["sha256"] == sha256_file(Path(REFERENCES_DIR) / filename)
        assert entry["records"] == len(entry["record_ids"])


def test_validated_profile_matches_manifest():
    manifest = load_reference_manifest()
    validated = {
        entry["file"] for entry in manifest["references"] if entry["profile"] == "validated"
    }
    assert validated == set(VALIDATED_REFERENCE_FILES)


def test_public_acceptance_fixtures_match_checksums():
    directory = Path(__file__).parent / "test-data" / "acceptance"
    for line in (directory / "SHA256SUMS").read_text().splitlines():
        expected, filename = line.split()
        assert sha256_file(directory / filename) == expected
