"""Plain-language report generated from the executable expert-rule YAML."""

from __future__ import annotations

from typing import Any

from markdown_it import MarkdownIt

from .element_definitions import load_tn123_definitions
from .unit_definitions import load_unit_definitions

ROLE_EXPLANATIONS = {
    "terminal_IR": "Short terminal inverted repeat marking an element end",
    "blaTEM": "Beta-lactam resistance gene",
    "tnpA": "Transposase gene used to move the element",
    "tnpR": "Resolvase gene used during transposition",
    "res": "Resolution site acted on by the resolvase",
    "tnpM": "Tn21-family transposition-region gene",
    "mcp": "Methyl-accepting chemotaxis-like protein gene",
    "tetR": "Regulator of tetracycline resistance",
    "tetA": "Tetracycline resistance gene",
    "pecM": "Membrane-protein gene in the Tn1721 variable region",
    "tet_tnp_junction": "Sequence junction linking the tet region to duplicated transposition DNA",
}


def _role_description(role: str, components: list[dict[str, Any]]) -> str:
    if role in ROLE_EXPLANATIONS:
        return ROLE_EXPLANATIONS[role]
    names = sorted({str(item["name"]) for item in components if item["role"] == role})
    return ", ".join(names) if names else role.replace("_", " ")


def _component_table(
    required_counts: dict[str, int],
    components: list[dict[str, Any]],
) -> list[str]:
    lines = [
        "| Component | Minimum number | What it represents |",
        "| --- | ---: | --- |",
    ]
    for role, count in required_counts.items():
        lines.append(f"| `{role}` | {count} | {_role_description(role, components)} |")
    return lines


def _tn123_section(document: dict[str, Any]) -> list[str]:
    family = document["family"]
    grammar = document["grammar"]
    classification = document["classification"]
    types = document["types"]
    definitions = document["definitions"]
    components = [
        component
        for type_definition in types.values()
        for component in type_definition["components"]
    ]
    required_counts = {
        item["role"]: int(item["minimum_count"])
        for item in grammar["required_components"]
    }
    lines = [
        "## Tn1/Tn2/Tn3 group",
        "",
        "### Plain-language decision",
        "",
        "The software first looks for the parts of the transposon independently. "
        "It calls a complete Tn1/Tn2/Tn3-group locus only when all required parts "
        "occur in the correct order between two ends. It then compares the local "
        "versions of the genes and resolution site to decide whether the locus is "
        "most consistent with Tn1, Tn2 or Tn3. A whole-element match is checked "
        "afterwards and is not allowed to create the initial family call.",
        "",
        f"Expected element boundaries: two {family['terminal_ir_length']} bp terminal "
        f"inverted repeats. A {family['tsd_length']} bp target-site duplication is "
        "supporting boundary evidence, not a substitute for the required internal parts.",
        "",
        "### Evidence that must be found",
        "",
        *_component_table(required_counts, components),
        "",
        "Required order on the forward strand:",
        "",
        f"`{' → '.join(grammar['forward_order'])}`",
        "",
        "The reverse-complement order is also accepted.",
        "",
        "### How Tn1, Tn2 and Tn3 are separated",
        "",
        "| Reported type | Canonical example | Resistance allele | Local evidence compared |",
        "| --- | --- | --- | --- |",
    ]
    for type_name, definition in types.items():
        lines.append(
            f"| {type_name} | `{definition['source_accession']}` | "
            f"{definition['bla_allele']} | blaTEM, tnpR, res and tnpA profiles |"
        )
    assignment = classification["type_assignment"]
    lines.extend([
        "",
        f"The best local component profile must score at least "
        f"{assignment['minimum_component_profile_score_percent']:g}% and beat the "
        f"next-best type by at least {assignment['minimum_type_margin_percent']:g} percentage "
        "points. The res and tnpA profiles receive the greatest weight. If this is "
        "ambiguous, the locus remains an unresolved member of the group.",
        "",
        "### Reviewed exact names and subtypes",
        "",
        "| Name | Parent type | Accession | What makes this exact definition distinct |",
        "| --- | --- | --- | --- |",
    ])
    for definition in definitions.values():
        differences = definition.get("known_differences_from_parent", [])
        distinction = "; ".join(str(value) for value in differences) or str(
            definition["expert_rule"]
        )
        lines.append(
            f"| {definition['display_name']} | {definition['type']} | "
            f"`{definition['source_accession']}` | {distinction} |"
        )
    lines.extend([
        "",
        "### What the program will not do",
        "",
        f"{document['related_element_policy']['expert_rule']}",
        "",
        "Examples deliberately not forced into a Tn1/Tn2/Tn3 name: "
        + ", ".join(document["related_element_policy"]["examples_expected_not_to_be_named_tn123"])
        + ".",
    ])
    return lines


def _unit_section(family: str, definition: dict[str, Any]) -> list[str]:
    reference = definition["canonical_reference"]
    required_counts = {
        str(role): int(count) for role, count in definition["required_counts"].items()
    }
    optional = [str(value) for value in definition.get("optional_components", [])]
    span = definition["candidate_span_bp"]
    lines = [
        f"## {family}",
        "",
        "### Plain-language decision",
        "",
        str(definition["expert_rule"]),
        "",
        f"The candidate must span approximately {int(span['minimum']):,}–"
        f"{int(span['maximum']):,} bp and contain {definition['terminal_ir_length']} bp "
        f"terminal inverted repeats. A {definition['tsd_length']} bp target-site "
        "duplication may support the inferred boundaries.",
        "",
        "### Evidence that must be found",
        "",
        *_component_table(required_counts, definition["components"]),
        "",
        "Required order on the forward strand:",
        "",
        f"`{' → '.join(definition['required_order'])}`",
        "",
    ]
    if optional:
        lines.extend([
            "Optional features are drawn when detected but do not create the family call: "
            + ", ".join(f"`{value}`" for value in optional)
            + ".",
            "",
        ])
    else:
        lines.extend(["No optional components are declared for this rule.", ""])
    lines.extend([
        "### Naming decision",
        "",
        f"Complete component grammar produces **{family}-like**. The unqualified name "
        f"**{family}** requires an exact, complete secondary match to the reviewed "
        f"example `{reference['accession']}`. Missing required components remain visible "
        "but are not promoted to a complete family call.",
    ])
    return lines


def expert_rules_as_markdown() -> str:
    """Render all executable expert rules for a non-specialist reviewer."""
    tn123 = load_tn123_definitions()
    units = load_unit_definitions()
    lines = [
        "# Expert rule review",
        "",
        f"Rule versions: Tn1/Tn2/Tn3 `{tn123['definition_version']}`; reviewed units "
        f"`{units['definition_version']}`.",
        "",
        "This report is generated directly from the same YAML definitions used by "
        "the program. It translates those executable rules into review language; it "
        "does not maintain a separate manual description.",
        "",
        "## Scope of this review",
        "",
        "This edition covers the six families with executable component-driven rules: "
        "Tn1, Tn2, Tn3, Tn21, Tn1721 and Tn1722. Other roadmap targets are not included "
        "here until their component grammars are implemented; reference-supported names "
        "must not be mistaken for reviewed expert rules.",
        "",
        "## How a sequence becomes a call",
        "",
        "1. **Detect parts independently.** Search for terminal repeats, genes, "
        "resolution sites, resistance genes and diagnostic junctions.",
        "2. **Assemble a candidate.** Check copy number, order, orientation and maximum span.",
        "3. **Apply the family rule.** A complete component grammar can create a "
        "family or family-like call without a whole-element match.",
        "4. **Assign a type where rules permit it.** Compare the locally detected "
        "components with the profiles declared for each type.",
        "5. **Compare with canonical sequences second.** This confirms an exact known "
        "definition and records differences; it does not drive discovery.",
        "6. **Keep incomplete evidence.** Components and fragments are drawn even when "
        "there is not enough evidence to name a complete element.",
        "",
        "## Meaning of result labels",
        "",
        "| Label | Meaning |",
        "| --- | --- |",
        "| Exact name, for example **Tn1** | Complete rule plus a 100% complete, "
        "mismatch-free match to a declared exact definition |",
        "| **-like** | Complete component rule, but not identical to a declared exact definition |",
        "| Unresolved group member | Complete family structure but insufficient evidence "
        "to choose a type |",
        "| Fragment | Some supported parts are present, but required grammar is incomplete |",
        "",
        *_tn123_section(tn123),
    ]
    for family, definition in units["units"].items():
        lines.extend(["", *_unit_section(family, definition)])
    lines.extend([
        "",
        "## Questions for expert review",
        "",
        "For each family, please check:",
        "",
        "- Are the required components and their minimum copy numbers correct?",
        "- Is the required order correct on the element-relative forward strand?",
        "- Should any currently optional component be required, or vice versa?",
        "- Are the acceptable candidate-size limits biologically reasonable?",
        "- Is the canonical accession appropriate only for exact naming?",
        "- Which known variants need their own exact subtype definition?",
        "- Which related elements should be explicit negative controls?",
        "",
        "## Source YAML",
        "",
        "The executable sources are `matryoshka/tn123_definitions.yaml` and "
        "`matryoshka/unit_definitions.yaml`. Edit and validate those files rather than "
        "editing this generated report.",
        "",
    ])
    return "\n".join(lines)


def expert_rules_as_html() -> str:
    """Render the plain-language review as a standalone HTML document."""
    body = MarkdownIt("commonmark", {"html": False}).enable("table").render(
        expert_rules_as_markdown()
    )
    return f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Expert rule review</title>
<style>
:root {{ --ink:#24313a; --muted:#65727b; --accent:#3b7b82; --paper:#fffdf8;
  --line:#dce3df; --soft:#edf5f2; }}
* {{ box-sizing:border-box; }}
body {{ margin:0; background:#eef1ee; color:var(--ink); font:17px/1.58 system-ui,sans-serif; }}
main {{ max-width:1080px; margin:32px auto; padding:48px 60px; background:var(--paper);
  box-shadow:0 8px 30px #26352f18; }}
h1 {{ margin-top:0; color:var(--accent); font-size:2.25rem; }}
h2 {{ margin-top:2.4rem; padding-top:1rem; border-top:2px solid var(--line); color:#315f65; }}
h3 {{ color:#446a6e; }}
table {{ width:100%; border-collapse:collapse; margin:1rem 0 1.5rem; font-size:.94rem; }}
th,td {{ padding:.65rem .75rem; border:1px solid var(--line); vertical-align:top; text-align:left; }}
th {{ background:var(--soft); }}
code {{ background:#edf0ed; padding:.12rem .3rem; border-radius:3px; }}
li {{ margin:.35rem 0; }}
@media (max-width:700px) {{ main {{ margin:0; padding:24px; }} table {{ font-size:.82rem; }} }}
@media print {{ body {{ background:white; }} main {{ margin:0; box-shadow:none; max-width:none; }} }}
</style>
</head>
<body><main>{body}</main></body>
</html>
"""
