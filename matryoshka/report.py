"""Versioned, reproducible output documents for the public alpha workflow."""

from __future__ import annotations

import hashlib
import json
import platform
import re
import shlex
import sys
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Any
from urllib.parse import quote

from .detect import MGEFeature
from .reference_manifest import REFERENCE_DATABASE_VERSION

OUTPUT_SCHEMA_VERSION = "1.0"


def _package_version() -> str:
    try:
        return version("matryoshka")
    except PackageNotFoundError:
        return "0.0.0+local"


def file_sha256(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _slug(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "-", value).strip("-") or "feature"


def _feature_identifier(
    seqid: str,
    feature: MGEFeature,
    parent_id: str | None,
    ordinal: int,
) -> str:
    identity = "|".join([
        seqid,
        parent_id or "root",
        feature.element_type,
        feature.family,
        feature.name,
        str(feature.start),
        str(feature.end),
        feature.strand,
        str(ordinal),
    ])
    suffix = hashlib.sha1(identity.encode("utf-8")).hexdigest()[:10]
    return f"mge-{_slug(seqid)}-{_slug(feature.element_type)}-{suffix}"


def _feature_document(
    feature: MGEFeature,
    seqid: str,
    parent_id: str | None,
    ordinal: int,
) -> dict[str, Any]:
    feature_id = _feature_identifier(seqid, feature, parent_id, ordinal)
    return {
        "id": feature_id,
        "parent_id": parent_id,
        "element_type": feature.element_type,
        "family": feature.family,
        "name": feature.name,
        "start": feature.start,
        "end": feature.end,
        "strand": feature.strand,
        "tsd_length": feature.tsd_length,
        "tsd_seq": feature.tsd_seq,
        "ir_left": feature.ir_left,
        "ir_right": feature.ir_right,
        "score": feature.score,
        "attributes": dict(feature.attributes) if feature.attributes else {},
        "children": [
            _feature_document(child, seqid, feature_id, child_ordinal)
            for child_ordinal, child in enumerate(feature.children, start=1)
        ],
    }


def annotation_document(
    records: list[tuple[str, int, list[MGEFeature]]],
    *,
    input_path: str | Path,
    profile: str,
    command: list[str] | None = None,
    detector_outputs: dict[str, str] | None = None,
) -> dict[str, Any]:
    """Build the stable v1 JSON annotation envelope."""
    path = Path(input_path)
    return {
        "schema_version": OUTPUT_SCHEMA_VERSION,
        "generated_by": {
            "name": "matryoshka",
            "version": _package_version(),
        },
        "reference_database": {
            "version": REFERENCE_DATABASE_VERSION,
            "profile": profile,
        },
        "run": {
            "input": str(path),
            "input_sha256": file_sha256(path),
            "command": shlex.join(command or sys.argv),
            "python": platform.python_version(),
            "platform": platform.platform(),
            "detector_outputs": detector_outputs or {},
        },
        "records": [
            {
                "id": seqid,
                "length": length,
                "features": [
                    _feature_document(feature, seqid, None, ordinal)
                    for ordinal, feature in enumerate(roots, start=1)
                ],
            }
            for seqid, length, roots in records
        ],
    }


def annotation_json(*args: Any, **kwargs: Any) -> str:
    return json.dumps(annotation_document(*args, **kwargs), indent=2) + "\n"


def _gff_escape(value: object) -> str:
    """Percent-encode a GFF3 field while retaining common identifier punctuation."""
    return quote(str(value), safe="._:-+|")


def _gff_rows(
    feature: MGEFeature,
    seqid: str,
    parent_id: str | None,
    ordinal: int,
) -> list[str]:
    feature_id = _feature_identifier(seqid, feature, parent_id, ordinal)
    attrs = [
        f"ID={_gff_escape(feature_id)}",
        f"Name={_gff_escape(feature.name)}",
        f"family={_gff_escape(feature.family)}",
    ]
    if parent_id:
        attrs.append(f"Parent={_gff_escape(parent_id)}")
    if feature.tsd_seq:
        attrs.append(f"tsd={_gff_escape(feature.tsd_seq)}")
    if feature.ir_left:
        attrs.append(f"ir_left={_gff_escape(feature.ir_left)}")
    if feature.ir_right:
        attrs.append(f"ir_right={_gff_escape(feature.ir_right)}")
    for key, value in feature.attributes.items():
        if key != "seqid" and value not in (None, ""):
            attrs.append(f"{_gff_escape(key)}={_gff_escape(value)}")
    score = str(feature.score) if feature.score is not None else "."
    rows = [
        "\t".join([
            seqid,
            "Matryoshka",
            feature.element_type,
            str(feature.start),
            str(feature.end),
            score,
            feature.strand,
            ".",
            ";".join(attrs),
        ])
    ]
    for child_ordinal, child in enumerate(feature.children, start=1):
        rows.extend(_gff_rows(child, seqid, feature_id, child_ordinal))
    return rows


def annotation_gff3(records: list[tuple[str, int, list[MGEFeature]]]) -> str:
    """Render one valid multi-record GFF3 document with globally stable IDs."""
    rows = ["##gff-version 3"]
    for seqid, length, roots in records:
        rows.append(f"##sequence-region {_gff_escape(seqid)} 1 {length}")
        for ordinal, feature in enumerate(roots, start=1):
            rows.extend(_gff_rows(feature, seqid, None, ordinal))
    return "\n".join(rows) + "\n"


def count_features(features: list[MGEFeature]) -> int:
    return sum(1 + count_features(feature.children) for feature in features)
