"""Evidence bundle proving the Tn1/Tn2/Tn3 annotation chain.

The proof is intentionally redundant with the full annotation document: it
extracts the small set of facts a reviewer needs to audit the path from raw
component sequence matches, through assembly grammar, to a named element and
its generated figures.
"""

from __future__ import annotations

import csv
import html
import io
import json
from pathlib import Path
from typing import Any

from .detect import MGEFeature

PROOF_SCHEMA_VERSION = "1.0"
TN123_FAMILIES = frozenset({"Tn1", "Tn2", "Tn3", "Tn3_family"})


def _walk(features: list[MGEFeature]) -> list[MGEFeature]:
    found: list[MGEFeature] = []
    for feature in features:
        found.append(feature)
        found.extend(_walk(feature.children))
    return found


def _tn123_parents(roots: list[MGEFeature]) -> list[MGEFeature]:
    return [
        feature
        for feature in _walk(roots)
        if feature.element_type == "transposon"
        and (
            feature.family in TN123_FAMILIES
            or feature.attributes.get("tn123_canonical") == "true"
        )
        and feature.attributes.get("component_assembly_status")
    ]


def _detected_components(parent: MGEFeature) -> list[MGEFeature]:
    return sorted(
        (
            feature
            for feature in _walk(parent.children)
            if feature.attributes.get("source") == "tn123_component_scan"
        ),
        key=lambda feature: (feature.start, feature.end, feature.name),
    )


def _component_document(feature: MGEFeature) -> dict[str, Any]:
    attributes = feature.attributes
    return {
        "role": attributes.get("component_role"),
        "name": feature.name,
        "start": feature.start,
        "end": feature.end,
        "strand": feature.strand,
        "identity": attributes.get("blast_identity"),
        "coverage": attributes.get("blast_coverage"),
        "status": attributes.get("structural_status"),
        "component_status": attributes.get("component_status"),
        "evidence_class": attributes.get("evidence_class"),
        "reference": attributes.get("component_reference"),
        "reference_segments": attributes.get("reference_segments", []),
        "inserted_bases": int(attributes.get("inserted_bases", 0) or 0),
    }


def _component_documents(parent: MGEFeature) -> list[dict[str, Any]]:
    """Return detected components even when hierarchy precedence made them siblings."""
    detected = _detected_components(parent)
    if detected:
        return [_component_document(component) for component in detected]
    evidence = parent.attributes.get("component_evidence", [])
    if not isinstance(evidence, list):
        return []
    documents: list[dict[str, Any]] = []
    for component in evidence:
        if not isinstance(component, dict):
            continue
        documents.append({
            "role": component.get("role"),
            "name": component.get("name"),
            "start": component.get("start"),
            "end": component.get("end"),
            "strand": component.get("strand"),
            "identity": component.get("identity"),
            "coverage": component.get("coverage"),
            "status": component.get("status"),
            "component_status": component.get("component_status"),
            "evidence_class": component.get("evidence_class"),
            "reference": component.get("reference"),
            "reference_segments": component.get("reference_segments", []),
            "inserted_bases": int(component.get("inserted_bases", 0) or 0),
        })
    return documents


def build_tn123_proof(
    records: list[tuple[str, int, list[MGEFeature]]],
    output_paths: dict[tuple[str, int, int, str], dict[str, str]],
) -> dict[str, Any]:
    """Build a machine-checkable proof document for every Tn1/2/3 locus."""
    record_documents: list[dict[str, Any]] = []
    pass_count = 0
    partial_count = 0
    fail_count = 0

    for seqid, length, roots in records:
        loci: list[dict[str, Any]] = []
        for parent in _tn123_parents(roots):
            attributes = parent.attributes
            requirements = attributes.get("component_requirements", {})
            if not isinstance(requirements, dict):
                requirements = {}
            components = _component_documents(parent)
            best_match = attributes.get("best_match")
            naming_evidence = attributes.get("naming_evidence", [])
            if not isinstance(naming_evidence, list):
                naming_evidence = []
            checks = {
                "components_detected_from_sequence": (
                    len(components) >= 6
                    and all(
                        component.get("evidence_class") == "sequence_detected"
                        for component in components
                    )
                ),
                "all_required_components_present": (
                    bool(requirements)
                    and all(requirements.get(role) is True for role in requirements)
                ),
                "component_order_and_orientation_valid": (
                    attributes.get("component_order_valid") is True
                ),
                "component_grammar_complete": (
                    attributes.get("component_assembly_status") == "complete"
                ),
                "known_element_reference_match": (
                    best_match in {"Tn1", "Tn2", "Tn3"}
                    and "whole_locus_alignment" in naming_evidence
                ),
                "no_reference_projected_components": (
                    int(attributes.get("reference_projected_component_count", 0) or 0)
                    == 0
                ),
            }
            expected_partial = (
                attributes.get("fragment") is True
                and bool(components)
                and checks["known_element_reference_match"]
                and checks["no_reference_projected_components"]
            )
            verdict = (
                "PASS" if all(checks.values())
                else "PARTIAL" if expected_partial
                else "FAIL"
            )
            if verdict == "PASS":
                pass_count += 1
            elif verdict == "PARTIAL":
                partial_count += 1
            else:
                fail_count += 1
            key = (seqid, parent.start, parent.end, parent.name)
            loci.append({
                "verdict": verdict,
                "record": seqid,
                "call": parent.name,
                "family": parent.family,
                "start": parent.start,
                "end": parent.end,
                "strand": parent.strand,
                "known_element_match": {
                    "best_match": best_match,
                    "closest_definition": attributes.get("closest_definition"),
                    "definition_id": attributes.get("definition_id"),
                    "definition_version": attributes.get("definition_version"),
                    "defined_type": attributes.get("defined_type"),
                    "defined_subtype": attributes.get("defined_subtype"),
                    "reference_id": attributes.get("reference_id"),
                    "source_accession": attributes.get("source_accession"),
                    "whole_locus_identity": attributes.get("blast_identity"),
                    "whole_locus_coverage": attributes.get("blast_coverage"),
                    "mismatch_bases": attributes.get("mismatch_bases", 0),
                    "inserted_bases": attributes.get("inserted_bases", 0),
                    "deleted_bases": attributes.get("deleted_bases", 0),
                    "variant_status": attributes.get("variant_status"),
                    "structural_status": attributes.get("structural_status"),
                    "expert_rule": attributes.get("expert_rule"),
                    "known_differences_from_parent": attributes.get(
                        "known_differences_from_parent", []
                    ),
                    "naming_evidence": naming_evidence,
                },
                "assembly": {
                    "status": attributes.get("component_assembly_status"),
                    "order_valid": attributes.get("component_order_valid"),
                    "requirements": requirements,
                    "detected_component_count": attributes.get(
                        "detected_component_count", 0
                    ),
                    "reference_projected_component_count": attributes.get(
                        "reference_projected_component_count", 0
                    ),
                    "reference_projected_context_feature_count": attributes.get(
                        "reference_projected_context_feature_count", 0
                    ),
                },
                "checks": checks,
                "components": [
                    component for component in components
                ],
                "outputs": output_paths.get(key, {}),
            })
        record_documents.append({"id": seqid, "length": length, "loci": loci})

    total = pass_count + partial_count + fail_count
    if total == 0:
        status = "NO_TN123_LOCI"
    elif fail_count:
        status = "FAIL"
    elif pass_count and partial_count:
        status = "PASS_WITH_PARTIAL_TN123_EVIDENCE"
    elif partial_count:
        status = "PARTIAL_TN123_EVIDENCE"
    else:
        status = "PASS"
    return {
        "schema_version": PROOF_SCHEMA_VERSION,
        "purpose": (
            "Audit arbitrary FASTA sequence through independent component detection, "
            "grammar assembly, known-element matching and generated outputs."
        ),
        "summary": {
            "status": status,
            "tn123_loci": total,
            "passed": pass_count,
            "partial": partial_count,
            "failed": fail_count,
        },
        "records": record_documents,
    }


def proof_components_tsv(proof: dict[str, Any]) -> str:
    """Render one row per independently sequence-detected component."""
    output = io.StringIO()
    fields = [
        "record", "locus_call", "locus_start", "locus_end", "locus_strand",
        "role", "component", "start", "end", "strand", "identity", "coverage",
        "status", "component_status", "evidence_class", "inserted_bases", "reference",
    ]
    writer = csv.DictWriter(output, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    for record in proof["records"]:
        for locus in record["loci"]:
            for component in locus["components"]:
                writer.writerow({
                    "record": record["id"],
                    "locus_call": locus["call"],
                    "locus_start": locus["start"],
                    "locus_end": locus["end"],
                    "locus_strand": locus["strand"],
                    "role": component["role"],
                    "component": component["name"],
                    "start": component["start"],
                    "end": component["end"],
                    "strand": component["strand"],
                    "identity": component["identity"],
                    "coverage": component["coverage"],
                    "status": component["status"],
                    "component_status": component["component_status"],
                    "evidence_class": component["evidence_class"],
                    "inserted_bases": component["inserted_bases"],
                    "reference": component["reference"],
                })
    return output.getvalue()


def proof_matches_tsv(proof: dict[str, Any]) -> str:
    """Render one row per assembled and matched Tn1/2/3 locus."""
    output = io.StringIO()
    fields = [
        "verdict", "record", "call", "start", "end", "strand", "best_match",
        "defined_type", "defined_subtype", "definition_id", "definition_version",
        "reference_id", "source_accession", "whole_locus_identity",
        "whole_locus_coverage", "mismatch_bases", "inserted_bases",
        "deleted_bases", "variant_status", "structural_status",
        "grammar_status", "order_valid", "detected_components",
        "projected_components", "projected_context_features",
    ]
    writer = csv.DictWriter(output, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    for record in proof["records"]:
        for locus in record["loci"]:
            match = locus["known_element_match"]
            assembly = locus["assembly"]
            writer.writerow({
                "verdict": locus["verdict"],
                "record": record["id"],
                "call": locus["call"],
                "start": locus["start"],
                "end": locus["end"],
                "strand": locus["strand"],
                "best_match": match["best_match"],
                "defined_type": match["defined_type"],
                "defined_subtype": match["defined_subtype"],
                "definition_id": match["definition_id"],
                "definition_version": match["definition_version"],
                "reference_id": match["reference_id"],
                "source_accession": match["source_accession"],
                "whole_locus_identity": match["whole_locus_identity"],
                "whole_locus_coverage": match["whole_locus_coverage"],
                "mismatch_bases": match["mismatch_bases"],
                "inserted_bases": match["inserted_bases"],
                "deleted_bases": match["deleted_bases"],
                "variant_status": match["variant_status"],
                "structural_status": match["structural_status"],
                "grammar_status": assembly["status"],
                "order_valid": assembly["order_valid"],
                "detected_components": assembly["detected_component_count"],
                "projected_components": assembly[
                    "reference_projected_component_count"
                ],
                "projected_context_features": assembly[
                    "reference_projected_context_feature_count"
                ],
            })
    return output.getvalue()


def proof_html(proof: dict[str, Any], title: str = "Matryoshka proof report") -> str:
    """Render a portable report whose images link to the result directory."""
    summary = proof["summary"]
    cards: list[str] = []
    for record in proof["records"]:
        for locus in record["loci"]:
            match = locus["known_element_match"]
            assembly = locus["assembly"]
            outputs = locus["outputs"]
            rows = "".join(
                "<tr>"
                f"<td>{html.escape(str(component['role']))}</td>"
                f"<td>{html.escape(str(component['name']))}</td>"
                f"<td>{component['start']}..{component['end']}</td>"
                f"<td>{html.escape(str(component['strand']))}</td>"
                f"<td>{component['identity']}%</td>"
                f"<td>{component['coverage']}%</td>"
                f"<td>{html.escape(str(component['status']))}</td>"
                "</tr>"
                for component in locus["components"]
            )
            links = " ".join(
                f'<a href="../{html.escape(path)}">{html.escape(name.replace("_", " "))}</a>'
                for name, path in outputs.items()
            )
            mara = outputs.get("mara")
            hierarchy = outputs.get("hierarchy")
            figures = ""
            if mara:
                figures += (
                    '<h3>MARA locus map</h3>'
                    f'<img src="../{html.escape(mara)}" alt="MARA locus map">'
                )
            if hierarchy:
                figures += (
                    '<h3>Original hierarchical view</h3>'
                    f'<img src="../{html.escape(hierarchy)}" alt="Hierarchical sequence view">'
                )
            differences = match.get("known_differences_from_parent", [])
            difference_text = (
                "<p><b>Reviewed subtype differences:</b> "
                + html.escape("; ".join(str(item) for item in differences))
                + "</p>"
                if isinstance(differences, list) and differences
                else ""
            )
            cards.append(
                '<section class="card">'
                f'<div class="verdict {locus["verdict"].lower()}">{locus["verdict"]}</div>'
                f'<h2>{html.escape(str(locus["call"]))}</h2>'
                f'<p>{html.escape(str(record["id"]))}: {locus["start"]}..{locus["end"]} '
                f'({html.escape(str(locus["strand"]))})</p>'
                '<div class="facts">'
                f'<span>known match <b>{html.escape(str(match["best_match"]))}</b></span>'
                f'<span>definition <b>{html.escape(str(match["definition_id"]))}</b></span>'
                f'<span>subtype <b>{html.escape(str(match["defined_subtype"] or "canonical"))}</b></span>'
                f'<span>whole-locus identity <b>{match["whole_locus_identity"]}%</b></span>'
                f'<span>differences <b>{match["mismatch_bases"]} substitutions; '
                f'{match["inserted_bases"]} inserted bp; {match["deleted_bases"]} deleted bp</b></span>'
                f'<span>grammar <b>{html.escape(str(assembly["status"]))}</b></span>'
                f'<span>detected components <b>{assembly["detected_component_count"]}</b></span>'
                '</div>'
                f'{difference_text}'
                '<table><thead><tr><th>Role</th><th>Call</th><th>Position</th>'
                '<th>Strand</th><th>Identity</th><th>Coverage</th><th>Status</th>'
                f'</tr></thead><tbody>{rows}</tbody></table>'
                f'<p class="links">{links}</p>{figures}</section>'
            )
    embedded_proof = json.dumps(proof, sort_keys=True).replace("</", "<\\/")
    return f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width">
<title>{html.escape(title)}</title>
<style>
body{{font:15px/1.45 system-ui,sans-serif;margin:0;background:#f4f5f6;color:#182026}}
main{{max-width:1280px;margin:auto;padding:32px}}h1,h2,h3{{margin:.3em 0}}
.summary,.card{{background:white;border:1px solid #d9dde1;border-radius:10px;padding:20px;margin:0 0 22px}}
.verdict{{float:right;font-weight:800;padding:5px 10px;border-radius:999px}}
.pass{{background:#d7f5e5;color:#096b3b}}.partial{{background:#fff0c2;color:#765600}}.fail{{background:#fde0e0;color:#9b1c1c}}
.facts{{display:flex;gap:12px;flex-wrap:wrap;margin:16px 0}}.facts span{{background:#eef1f3;padding:7px 10px;border-radius:5px}}
table{{border-collapse:collapse;width:100%;font-size:13px}}th,td{{border:1px solid #d9dde1;padding:7px;text-align:left}}th{{background:#eef1f3}}
img{{display:block;width:100%;height:auto;border:1px solid #d9dde1;background:white}}
a{{color:#075fa8}}.links{{display:flex;gap:14px;flex-wrap:wrap}}
</style></head><body><main>
<div class="summary"><h1>{html.escape(title)}</h1>
<p><b>{summary['status']}</b>: {summary['passed']} complete loci passed every component, grammar, match and projection check; {summary['partial']} expected partial loci retained incomplete evidence; {summary['failed']} loci failed validation.</p>
<p>This report is generated from the annotation evidence. It is not a hand-authored interpretation.</p></div>
{''.join(cards)}
<script type="application/json" id="matryoshka-proof">{embedded_proof}</script>
</main></body></html>
"""


def write_proof_bundle(
    directory: Path,
    proof: dict[str, Any],
    *,
    title: str,
) -> None:
    """Write the JSON, TSV and human-readable proof artefacts."""
    directory.mkdir(parents=True, exist_ok=True)
    (directory / "proof.json").write_text(
        json.dumps(proof, indent=2) + "\n", encoding="utf-8"
    )
    (directory / "components.tsv").write_text(
        proof_components_tsv(proof), encoding="utf-8"
    )
    (directory / "matches.tsv").write_text(
        proof_matches_tsv(proof), encoding="utf-8"
    )
    (directory / "report.html").write_text(
        proof_html(proof, title), encoding="utf-8"
    )
