"""Curated annotations for the canonical Tn1, Tn2 and Tn3 references.

The reference choices and feature layouts follow the material supplied by
expert-reviewed source material. Coordinates are 1-based and relative to each complete
transposon reference sequence.
"""

from __future__ import annotations

from dataclasses import dataclass

from .detect import MGEFeature
from .element_definitions import (
    tn123_components,
    tn123_definition,
    tn123_reference_metadata,
    tn123_rules,
    tn123_type,
)


@dataclass(frozen=True)
class InternalFeature:
    element_type: str
    name: str
    start: int
    end: int
    strand: str


@dataclass(frozen=True)
class Tn123ComponentReference:
    """Independently searchable component extracted from a canonical locus."""

    reference_id: str
    parent_reference: str
    parent_family: str
    role: str
    element_type: str
    name: str
    start: int
    end: int
    strand: str
    sequence: str


REFERENCE_METADATA: dict[str, dict[str, str]] = tn123_reference_metadata()


def _internal_features(reference_id: str) -> tuple[InternalFeature, ...]:
    return tuple(
        InternalFeature(
            str(component["element_type"]),
            str(component["name"]),
            int(component["start"]),
            int(component["end"]),
            str(component["strand"]),
        )
        for component in tn123_components(reference_id)
    )


def component_references(
    reference_id: str,
    sequence: str,
) -> list[Tn123ComponentReference]:
    """Return canonical components as independent nucleotide references."""
    metadata = REFERENCE_METADATA[reference_id]
    definition = tn123_definition(reference_id)
    if not definition["component_reference"]:
        return []
    family = metadata["family"]
    references: list[Tn123ComponentReference] = []
    definition_components = tn123_components(reference_id)
    for item, component in zip(
        _internal_features(reference_id), definition_components, strict=True
    ):
        start = _absolute_coordinate(item.start, len(sequence))
        end = _absolute_coordinate(item.end, len(sequence))
        name = item.name
        role = str(component["role"])
        safe_name = "".join(character if character.isalnum() else "_" for character in name)
        references.append(Tn123ComponentReference(
            reference_id=f"{reference_id}__{safe_name}__{start}_{end}",
            parent_reference=reference_id,
            parent_family=family,
            role=role,
            element_type=item.element_type,
            name=name,
            start=start,
            end=end,
            strand=item.strand,
            sequence=sequence[start - 1:end],
        ))
    return references


def _absolute_coordinate(value: int, reference_length: int) -> int:
    return value if value > 0 else reference_length + value + 1


def _project(parent: MGEFeature, start: int, end: int) -> tuple[int, int]:
    if parent.strand == "-":
        return parent.end - end + 1, parent.end - start + 1
    return parent.start + start - 1, parent.start + end - 1


def _map_reference_position(parent: MGEFeature, position: int) -> int | None:
    segments = parent.attributes.get("reference_segments")
    if not isinstance(segments, list):
        return None
    for segment in segments:
        if not isinstance(segment, dict):
            continue
        try:
            qstart = int(segment["qstart"])
            qend = int(segment["qend"])
            sstart = int(segment["sstart"])
            send = int(segment["send"])
        except (KeyError, TypeError, ValueError):
            continue
        if not min(sstart, send) <= position <= max(sstart, send):
            continue
        if send == sstart:
            return qstart
        slope = (qend - qstart) / (send - sstart)
        return round(qstart + (position - sstart) * slope)
    return None


def _project_reference_interval(
    parent: MGEFeature,
    start: int,
    end: int,
) -> tuple[int, int] | None:
    mapped_start = _map_reference_position(parent, start)
    mapped_end = _map_reference_position(parent, end)
    if mapped_start is not None and mapped_end is not None:
        interval_start, interval_end = sorted((mapped_start, mapped_end))
        return interval_start, interval_end
    if parent.attributes.get("reference_segments"):
        return None
    return _project(parent, start, end)


def _interruption_size(parent: MGEFeature, start: int, end: int) -> int:
    segments = parent.attributes.get("reference_segments")
    if not isinstance(segments, list) or len(segments) < 2:
        return 0
    parsed: list[tuple[int, int, int, int]] = []
    reference_length = int(parent.attributes.get("reference_length", 0))
    for segment in segments:
        if not isinstance(segment, dict):
            continue
        try:
            qstart, qend = sorted((int(segment["qstart"]), int(segment["qend"])))
            sstart, send = sorted((int(segment["sstart"]), int(segment["send"])))
        except (KeyError, TypeError, ValueError):
            continue
        if parent.strand == "-":
            sstart, send = reference_length - send + 1, reference_length - sstart + 1
        parsed.append((qstart, qend, sstart, send))
    parsed.sort()
    largest = 0
    for left, right in zip(parsed, parsed[1:], strict=False):
        query_gap = max(0, right[0] - left[1] - 1)
        reference_gap = max(0, right[2] - left[3] - 1)
        inserted = max(0, query_gap - reference_gap)
        junction = left[3] if parent.strand != "-" else reference_length - left[3] + 1
        if start <= junction <= end:
            largest = max(largest, inserted)
    return largest


def _project_strand(parent: MGEFeature, strand: str) -> str:
    if strand == "." or parent.strand != "-":
        return strand
    return "+" if strand == "-" else "-"


def curated_internal_features(parent: MGEFeature) -> list[MGEFeature]:
    """Project canonical internal features onto an exact full-length call."""
    reference_id = str(parent.attributes.get("reference_id", ""))
    if reference_id in REFERENCE_METADATA:
        definition = tn123_definition(reference_id)
        layout = _internal_features(reference_id)
    elif parent.family in {"Tn1", "Tn2", "Tn3"}:
        type_definition = tn123_type(parent.family)
        reference_id = str(type_definition["canonical_reference"])
        definition = tn123_definition(reference_id)
        layout = _internal_features(reference_id)
    else:
        return []
    try:
        coverage = float(parent.attributes.get("blast_subject_coverage", 0.0))
        identity = float(parent.attributes.get("blast_identity", 0.0))
        reference_length = int(parent.attributes["reference_length"])
    except (KeyError, TypeError, ValueError):
        return []
    complete = bool(parent.attributes.get("left_end_covered", True)) and bool(
        parent.attributes.get("right_end_covered", True)
    )
    if coverage < 95.0 or identity < 95.0 or not complete:
        return []

    seqid = parent.attributes.get("seqid", "")
    source_accession = parent.attributes.get("source_accession", "")
    children: list[MGEFeature] = []
    for item in layout:
        rel_start = _absolute_coordinate(item.start, reference_length)
        rel_end = _absolute_coordinate(item.end, reference_length)
        projected = _project_reference_interval(parent, rel_start, rel_end)
        if projected is None:
            continue
        start, end = projected
        name = item.name
        interruption = _interruption_size(parent, rel_start, rel_end)
        attributes: dict[str, object] = {
            "seqid": seqid,
            "source": "curated_reference",
            "source_accession": source_accession,
            "parent_transposon": parent.name,
        }
        if interruption:
            attributes["interrupted"] = True
            attributes["inserted_bases"] = interruption
        children.append(MGEFeature(
            element_type=item.element_type,
            family=parent.family,
            name=name,
            start=start,
            end=end,
            strand=_project_strand(parent, item.strand),
            attributes=attributes,
        ))
    for item in definition.get("additional_features", []):
        rel_start = _absolute_coordinate(int(item["start"]), reference_length)
        rel_end = _absolute_coordinate(int(item["end"]), reference_length)
        projected = _project_reference_interval(parent, rel_start, rel_end)
        if projected is None:
            continue
        start, end = projected
        context_attributes: dict[str, object] = {
            "seqid": seqid,
            "source": "expert_definition",
            "source_accession": source_accession,
            "parent_transposon": parent.name,
            "evidence_class": item.get("evidence", "definition_projection"),
            "definition_id": reference_id,
        }
        if item.get("status"):
            context_attributes["structural_status"] = item["status"]
            if item["status"] == "partial":
                context_attributes["fragment"] = True
        children.append(MGEFeature(
            element_type=str(item["element_type"]),
            family=str(item["family"]),
            name=str(item["name"]),
            start=start,
            end=end,
            strand=_project_strand(parent, str(item["strand"])),
            attributes=context_attributes,
        ))
    parent.attributes["curated_internal_features"] = True
    rules = tn123_rules()["family"]
    parent.attributes["ir_length"] = int(rules["terminal_ir_length"])
    parent.tsd_length = int(rules["tsd_length"])
    return children


def inserted_sequence_features(parent: MGEFeature) -> list[MGEFeature]:
    """Emit query gaps from whole-locus or independently detected component HSPs."""
    evidence_sources: list[tuple[object, int, str, str]] = [(
        parent.attributes.get("reference_segments"),
        int(parent.attributes.get("reference_length", 0) or 0),
        parent.strand,
        "reference_hsp_gap",
    )]
    component_evidence = parent.attributes.get("component_evidence", [])
    if isinstance(component_evidence, list):
        for component in component_evidence:
            if not isinstance(component, dict):
                continue
            evidence_sources.append((
                component.get("reference_segments"),
                int(component.get("reference_length", 0) or 0),
                str(component.get("strand", parent.strand)),
                "component_hsp_gap",
            ))
    out: list[MGEFeature] = []
    seen: set[tuple[int, int]] = set()
    for segments, reference_length, strand, source in evidence_sources:
        if not isinstance(segments, list) or len(segments) < 2:
            continue
        parsed: list[tuple[int, int, int, int]] = []
        for segment in segments:
            if not isinstance(segment, dict):
                continue
            try:
                qstart, qend = sorted((int(segment["qstart"]), int(segment["qend"])))
                sstart, send = sorted((int(segment["sstart"]), int(segment["send"])))
            except (KeyError, TypeError, ValueError):
                continue
            if strand == "-" and reference_length:
                sstart, send = (
                    reference_length - send + 1,
                    reference_length - sstart + 1,
                )
            parsed.append((qstart, qend, sstart, send))
        parsed.sort()
        for left, right in zip(parsed, parsed[1:], strict=False):
            query_gap = right[0] - left[1] - 1
            reference_gap = max(0, right[2] - left[3] - 1)
            inserted = query_gap - reference_gap
            if inserted <= 20:
                continue
            start = left[1] + 1
            end = right[0] - 1
            if (start, end) in seen:
                continue
            seen.add((start, end))
            out.append(MGEFeature(
                element_type="inserted_sequence",
                family="insertion",
                name=f"inserted sequence {len(out) + 1}",
                start=start,
                end=end,
                strand=".",
                attributes={
                    "seqid": parent.attributes.get("seqid", ""),
                    "source": source,
                    "parent_transposon": parent.name,
                    "inserted_bases": inserted,
                    "note": "sequence present between collinear component matches",
                },
            ))
    return out


def _required_component_counts() -> dict[str, int]:
    grammar = tn123_rules()["grammar"]
    return {
        str(component["role"]): int(component["minimum_count"])
        for component in grammar["required_components"]
    }


def _component_role(feature: MGEFeature) -> str:
    role = feature.attributes.get("component_role")
    if isinstance(role, str):
        return role
    if feature.element_type == "IR":
        return "terminal_IR"
    if feature.element_type == "AMR" and feature.name.startswith("blaTEM"):
        return "blaTEM"
    if feature.element_type == "res_site":
        return "res"
    return feature.name if feature.name in {"tnpA", "tnpR"} else ""


def _overlap_fraction(left: MGEFeature, right: MGEFeature) -> float:
    overlap = max(0, min(left.end, right.end) - max(left.start, right.start) + 1)
    shorter = min(left.end - left.start + 1, right.end - right.start + 1)
    return overlap / shorter if shorter else 0.0


def _components_for_parent(
    parent: MGEFeature,
    features: list[MGEFeature],
) -> list[MGEFeature]:
    return [
        feature
        for feature in features
        if feature.attributes.get("source") == "tn123_component_scan"
        and feature.attributes.get("seqid") == parent.attributes.get("seqid")
        and parent.start <= feature.start
        and feature.end <= parent.end
    ]


def _select_component_path(
    parent: MGEFeature,
    components: list[MGEFeature],
) -> tuple[list[MGEFeature], dict[str, bool], bool]:
    by_role: dict[str, list[MGEFeature]] = {}
    for component in components:
        role = _component_role(component)
        if role:
            by_role.setdefault(role, []).append(component)

    selected: list[MGEFeature] = []
    terminal_irs = sorted(by_role.get("terminal_IR", []), key=lambda item: item.start)
    if len(terminal_irs) >= 2:
        selected.extend([terminal_irs[0], terminal_irs[-1]])
    elif terminal_irs:
        selected.extend(terminal_irs)

    for role in ("blaTEM", "tnpR", "res", "tnpA"):
        candidates = by_role.get(role, [])
        if candidates:
            selected.append(max(
                candidates,
                key=lambda item: (
                    float(item.attributes.get("blast_coverage", 0)),
                    float(item.attributes.get("blast_identity", 0)),
                ),
            ))

    selected.sort(key=lambda item: (item.start, item.end))
    roles = [_component_role(feature) for feature in selected]
    grammar = tn123_rules()["grammar"]
    forward = [str(role) for role in grammar["forward_order"]]
    expected = forward if parent.strand != "-" else list(reversed(forward))
    required_counts = _required_component_counts()
    requirements = {
        role: len(by_role.get(role, [])) >= count
        for role, count in required_counts.items()
    }
    order_valid = roles == expected
    return selected, requirements, order_valid


def _backbone_role_group_calls(
    selected: list[MGEFeature],
) -> dict[str, dict[str, object]]:
    """Assign each discriminator gene to an expert-declared profile group."""
    type_rule = tn123_rules()["classification"]["type_assignment"]
    components = {_component_role(component): component for component in selected}
    minimum_score = float(type_rule["minimum_role_profile_score_percent"])
    minimum_margin = float(type_rule["minimum_role_group_margin_percent"])
    calls: dict[str, dict[str, object]] = {}
    for role in type_rule["required_discriminator_roles"]:
        component = components.get(str(role))
        if component is None:
            calls[str(role)] = {"call": "missing", "group_scores": {}, "margin": 0.0}
            continue
        matches = component.attributes.get("component_profile_matches", [])
        by_type = {
            str(match.get("type")): float(match.get("profile_score", 0.0))
            for match in matches
            if isinstance(match, dict)
        }
        group_scores = {
            str(group): round(max(by_type.get(str(member), 0.0) for member in members), 3)
            for group, members in type_rule["role_profile_groups"][role].items()
        }
        ranked = sorted(group_scores.items(), key=lambda item: item[1], reverse=True)
        best_group, best_score = ranked[0]
        runner_up = ranked[1][1] if len(ranked) > 1 else 0.0
        margin = round(best_score - runner_up, 3)
        call = (
            best_group
            if best_score >= minimum_score and margin >= minimum_margin
            else "ambiguous"
        )
        calls[str(role)] = {
            "call": call,
            "group_scores": group_scores,
            "margin": margin,
            "best_score": best_score,
        }
    return calls


def _apply_component_rule_classification(
    parent: MGEFeature,
    selected: list[MGEFeature],
    component_complete: bool,
) -> None:
    """Make the biological call from component rules, then attach ref context."""
    rules = tn123_rules()["classification"]
    type_rule = rules["type_assignment"]
    role_calls = _backbone_role_group_calls(selected)
    haplotype = {
        role: str(role_calls.get(role, {}).get("call", "missing"))
        for role in type_rule["required_discriminator_roles"]
    }
    best_type = next((
        str(type_name)
        for type_name, declared in type_rule["type_haplotypes"].items()
        if haplotype == {str(role): str(group) for role, group in declared.items()}
    ), "unresolved")
    all_roles_called = all(
        call not in {"missing", "ambiguous"} for call in haplotype.values()
    )
    qualified = component_complete and best_type != "unresolved"
    mosaic = component_complete and all_roles_called and not qualified

    parent.attributes.update({
        "classification_basis": "expert_component_rules",
        "rule_based_family_call": "Tn1/Tn2/Tn3 group",
        "rule_based_type_call": best_type if qualified else "unresolved",
        "backbone_role_group_calls": role_calls,
        "backbone_haplotype": haplotype,
        "backbone_haplotype_match": best_type,
        "naming_evidence": [
            "component_grammar",
            "backbone_haplotype",
        ],
    })

    if not qualified:
        parent.family = "Tn3_family"
        parent.name = (
            str(rules["fragment"]["visible_label"])
            if parent.attributes.get("fragment")
            else str(type_rule[
                "mosaic_haplotype_label" if mosaic else "unresolved_haplotype_label"
            ])
        )
        parent.attributes["variant_status"] = (
            "rule_incomplete" if not component_complete
            else "rule_haplotype_mosaic" if mosaic
            else "rule_type_unresolved"
        )
        return

    parent.family = best_type
    parent.name = str(rules["close_variant"]["label_template"]).format(
        type=best_type
    )
    parent.attributes["variant_status"] = "rule_based_type_candidate"

    reference_type = parent.attributes.get("reference_comparison_type")
    reference_status = parent.attributes.get("reference_comparison_status")
    if reference_type and reference_type != best_type:
        parent.attributes["reference_classification_conflict"] = True
        parent.attributes["note"] = (
            f"component rules support {best_type}; secondary whole-locus "
            f"comparison is closest to {reference_type}"
        )
        return
    if reference_type == best_type and reference_status == "exact_reference":
        parent.name = str(
            parent.attributes.get("reference_comparison_label", best_type)
        )
        parent.attributes["variant_status"] = "exact_reference_confirmed"
        parent.attributes["reference_confirmation"] = "exact_known_definition"
    elif reference_type == best_type:
        parent.attributes["reference_confirmation"] = "closest_reference_context"


def _record_component_assembly(
    parent: MGEFeature,
    components: list[MGEFeature],
) -> bool:
    selected, requirements, order_valid = _select_component_path(parent, components)
    component_complete = all(requirements.values()) and order_valid
    inserted_bases = sum(
        int(feature.attributes.get("inserted_bases", 0) or 0)
        for feature in selected
    )
    deleted_bases = sum(
        int(feature.attributes.get("deleted_bases", 0) or 0)
        for feature in selected
    )
    if component_complete and inserted_bases:
        structural_status = "complete_with_insertion"
    elif component_complete and deleted_bases:
        structural_status = "complete_with_deletion"
    elif component_complete:
        structural_status = "complete"
    else:
        structural_status = "partial_or_conflicting"
    if selected:
        irs = [feature for feature in selected if _component_role(feature) == "terminal_IR"]
        if len(irs) == 2:
            if parent.strand == "-":
                irs[0].name, irs[1].name = "IRR", "IRL"
            else:
                irs[0].name, irs[1].name = "IRL", "IRR"
    parent.attributes.update({
        "component_assembly_status": (
            "complete" if component_complete else "partial_or_conflicting"
        ),
        "component_order_valid": order_valid,
        "component_requirements": requirements,
        "detected_component_count": len(selected),
        "structural_status": structural_status,
        "inserted_bases": inserted_bases,
        "deleted_bases": deleted_bases,
        "component_evidence": [
            {
                "role": _component_role(feature),
                "name": feature.name,
                "start": feature.start,
                "end": feature.end,
                "strand": feature.strand,
                "evidence_class": feature.attributes.get("evidence_class"),
                "identity": feature.attributes.get("blast_identity"),
                "coverage": feature.attributes.get("blast_coverage"),
                "status": feature.attributes.get("structural_status"),
                "component_status": feature.attributes.get("component_status"),
                "reference": feature.attributes.get("component_reference"),
                "reference_length": feature.attributes.get("reference_length", 0),
                "reference_segments": feature.attributes.get("reference_segments", []),
                "inserted_bases": feature.attributes.get("inserted_bases", 0),
                "deleted_bases": feature.attributes.get("deleted_bases", 0),
            }
            for feature in selected
        ],
        "naming_evidence": [
            "component_grammar",
            *(
                ["whole_locus_alignment"]
                if parent.attributes.get("source") == "reference_scan"
                else []
            ),
        ],
    })
    if component_complete:
        parent.attributes["evidence_class"] = "sequence_detected_and_assembled"
    _apply_component_rule_classification(parent, selected, component_complete)
    return component_complete


def assemble_tn123_components(features: list[MGEFeature]) -> list[MGEFeature]:
    """Validate reference parents and emit reference-independent candidates."""
    parents = [
        feature
        for feature in features
        if feature.element_type == "transposon"
        and (
            feature.family in {"Tn1", "Tn2", "Tn3", "Tn3_family"}
            or feature.attributes.get("tn123_canonical") == "true"
        )
    ]
    components = [
        feature
        for feature in features
        if feature.attributes.get("source") == "tn123_component_scan"
    ]
    for parent in parents:
        _record_component_assembly(
            parent,
            _components_for_parent(parent, components),
        )

    emitted: list[MGEFeature] = []
    emitted_loci: set[tuple[object, int, int]] = set()
    grammar = tn123_rules()["grammar"]
    window = int(grammar["candidate_component_window_bp"])
    maximum_span = int(grammar["candidate_maximum_span_bp"])
    for tnpa in [feature for feature in components if _component_role(feature) == "tnpA"]:
        if any(_overlap_fraction(parent, tnpa) >= 0.8 for parent in parents):
            continue
        seqid = tnpa.attributes.get("seqid")
        nearby = [
            component
            for component in components
            if component.attributes.get("seqid") == seqid
            and component.start >= tnpa.start - int(
                window
            )
            and component.end <= tnpa.end + int(
                window
            )
        ]
        irs = sorted(
            (feature for feature in nearby if _component_role(feature) == "terminal_IR"),
            key=lambda item: item.start,
        )
        valid_candidates: list[MGEFeature] = []
        for left_index, left in enumerate(irs):
            for right in irs[left_index + 1:]:
                span = right.end - left.start + 1
                if span > maximum_span or not (
                    left.start <= tnpa.start and tnpa.end <= right.end
                ):
                    continue
                contained = [
                    component
                    for component in nearby
                    if left.start <= component.start and component.end <= right.end
                ]
                candidate = MGEFeature(
                    element_type="transposon",
                    family="Tn3_family",
                    name="Tn3-family unit",
                    start=left.start,
                    end=right.end,
                    strand=tnpa.strand,
                    tsd_length=5,
                    attributes={
                        "seqid": seqid,
                        "source": "component_assembly",
                        "best_match": "unresolved",
                        "variant_status": "component_assembled",
                    },
                )
                if _record_component_assembly(candidate, contained):
                    valid_candidates.append(candidate)
        if valid_candidates:
            candidate = min(valid_candidates, key=lambda item: item.end - item.start)
            locus_key = (seqid, candidate.start, candidate.end)
            if locus_key not in emitted_loci:
                emitted.append(candidate)
                emitted_loci.add(locus_key)
    # Repeated homologous IRs can create a valid larger pair around a complete
    # shorter locus. Prefer the tightest component-complete boundary while
    # retaining genuinely separate loci on the same record.
    selected: list[MGEFeature] = []
    for candidate in sorted(emitted, key=lambda item: item.end - item.start):
        if any(
            inner.attributes.get("seqid") == candidate.attributes.get("seqid")
            and candidate.start <= inner.start
            and inner.end <= candidate.end
            for inner in selected
        ):
            continue
        selected.append(candidate)
    return selected


def annotate_tn123(features: list[MGEFeature]) -> list[MGEFeature]:
    """Project only components that were not independently sequence-detected."""
    out: list[MGEFeature] = []
    parents = [feature for feature in features if feature.element_type == "transposon"]
    for parent in parents:
        projected = curated_internal_features(parent)
        detected = _components_for_parent(parent, features)
        detected_roles = {_component_role(feature) for feature in detected}
        retained = [
            feature
            for feature in projected
            if _component_role(feature) not in detected_roles
        ]
        for feature in retained:
            feature.attributes["evidence_class"] = "reference_projected"
        required_roles = set(_required_component_counts())
        projected_components = [
            feature for feature in retained
            if _component_role(feature) in required_roles
        ]
        projected_context = [
            feature for feature in retained
            if _component_role(feature) not in required_roles
        ]
        parent.attributes["curated_internal_features"] = bool(retained)
        parent.attributes["independently_detected_internal_features"] = bool(detected)
        # Only required grammar parts count as projected components. Reviewed
        # subtype context (for example the ISEcp1-associated insertion in
        # Tn2.1) is deliberately separate: it enriches the drawing but is not
        # used to manufacture the core Tn1/2/3 component proof.
        parent.attributes["reference_projected_component_count"] = len(
            projected_components
        )
        parent.attributes["reference_projected_context_feature_count"] = len(
            projected_context
        )
        out.extend(retained)
        out.extend(inserted_sequence_features(parent))
    return out
