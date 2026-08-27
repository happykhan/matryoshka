"""
Matryoshka CLI — detect and annotate nested MGEs in bacterial genomes.

Usage:
    matryoshka annotate <fasta> [--isescan <tsv>] [--amrfinder <tsv>]
                                [--integrons <file>]
                                [--format gff3|json|wolvercote|linear|mara|mara-table]
                                [--out <path>]
"""

from __future__ import annotations

import json
import os
import sys
from importlib.metadata import PackageNotFoundError
from importlib.metadata import version as _pkg_version
from pathlib import Path

import click
from Bio import SeqIO

from .boundaries import confirm_boundaries
from .component_catalog import catalog_as_json, catalog_as_tsv
from .confidence import assign_confidence
from .detect import (
    MGEFeature,
    detector_available,
    parse_amrfinder,
    parse_integron_finder,
    parse_isescan,
    run_amrfinder,
    run_integron_finder,
    run_isescan,
)
from .element_definitions import (
    definitions_as_json,
    definitions_as_markdown,
    definitions_as_yaml,
)
from .hierarchy import build_hierarchy
from .integron_structures import infer_integron_structures
from .mara_loci import extract_mara_loci
from .mara_table import to_mara_table_svg
from .mara_viz import to_mara_svg
from .output import to_genbank, to_gff3, to_json, to_wolvercote
from .partridge_units import annotate_partridge_units
from .proof import build_tn123_proof, write_proof_bundle
from .reference_scan import REFERENCE_PROFILES, blast_available, scan_all
from .report import annotation_document, annotation_gff3, count_features
from .tn123 import annotate_tn123, assemble_tn123_components
from .transposon import annotate_res_sites, infer_transposons
from .viz import to_linear_svg

try:
    _VERSION = _pkg_version("matryoshka")
except PackageNotFoundError:
    _VERSION = "0.0.0-dev"


@click.group()
@click.version_option(_VERSION, prog_name="matryoshka")
def cli() -> None:
    pass


@cli.command()
@click.option(
    "--format",
    "fmt",
    type=click.Choice(["tsv", "json"]),
    default="tsv",
    show_default=True,
    help="Catalogue output format.",
)
@click.option("--out", "-o", default="-", help="Output file (default: stdout).")
def catalog(fmt: str, out: str) -> None:
    """List the MARA raw components and compound assembly grammar."""
    data = catalog_as_json() if fmt == "json" else catalog_as_tsv()
    _emit(data, out, False)


@cli.command("definitions")
@click.option(
    "--format",
    "fmt",
    type=click.Choice(["yaml", "json", "markdown"]),
    default="markdown",
    show_default=True,
    help="Expert-definition output format.",
)
@click.option("--out", "-o", default="-", help="Output file (default: stdout).")
def definitions_command(fmt: str, out: str) -> None:
    """Show the biological rules used to call Tn1, Tn2 and Tn3."""
    exporters = {
        "yaml": definitions_as_yaml,
        "json": definitions_as_json,
        "markdown": definitions_as_markdown,
    }
    _emit(exporters[fmt](), out, False)


def _features_on_contig(
    all_features: list[MGEFeature], contig_id: str
) -> list[MGEFeature]:
    """Return features whose attributes.seqid matches contig_id.

    Parsers that don't record a seqid (older parsers, synthetic features)
    are treated as belonging to every contig — caller is responsible for
    single-contig assemblies in that case.
    """
    out: list[MGEFeature] = []
    for f in all_features:
        sid = f.attributes.get("seqid")
        if sid is None or sid == contig_id:
            out.append(f)
    return out


def _flatten_feature_inputs(features: list[MGEFeature]) -> list[MGEFeature]:
    """Preserve parser-supplied children before rebuilding the hierarchy."""
    flattened: list[MGEFeature] = []

    def _visit(feature: MGEFeature) -> None:
        flattened.append(feature)
        for child in feature.children:
            _visit(child)

    for feature in features:
        _visit(feature)
    return flattened


def _suppress_redundant_inference(features: list[MGEFeature]) -> list[MGEFeature]:
    """Drop rule-inferred transposons that overlap a BLAST-confirmed call of
    the same family. The BLAST boundaries are authoritative.
    """
    blast_by_family: dict[str, list[MGEFeature]] = {}
    exact_reference_is = [
        feature
        for feature in features
        if feature.element_type == "IS"
        and feature.attributes.get("source") == "reference_scan"
    ]
    for f in features:
        if f.attributes.get("source") == "reference_scan" and f.element_type == "transposon":
            blast_by_family.setdefault(f.family, []).append(f)

    def _overlap(a: MGEFeature, b: MGEFeature) -> bool:
        if a.end < b.start or b.end < a.start:
            return False
        ovl = min(a.end, b.end) - max(a.start, b.start)
        short = min(a.end - a.start, b.end - b.start)
        return short > 0 and ovl / short > 0.5

    def _same_reference_unit(a: MGEFeature, b: MGEFeature) -> bool:
        if a.attributes.get("seqid") != b.attributes.get("seqid"):
            return False
        if not _overlap(a, b):
            return False
        return (
            a.family == b.family
            or a.name == b.name
            or (a.family == "Tn3" and a.name in {b.family, b.name})
            or (b.family == "Tn3" and b.name in {a.family, a.name})
        )

    def _reference_rank(feature: MGEFeature) -> tuple[float, float, float, int]:
        provenance = 1.0 if feature.attributes.get("provenance") == "Sally_Partridge" else 0.0
        coverage = float(
            feature.attributes.get("blast_subject_coverage")
            or feature.attributes.get("blast_coverage")
            or 0
        )
        identity = float(feature.attributes.get("blast_identity") or 0)
        return provenance, coverage, identity, feature.end - feature.start

    reference_units = [
        feature
        for feature in features
        if feature.element_type == "transposon"
        and feature.attributes.get("source") == "reference_scan"
    ]
    selected_reference_ids: set[int] = set()
    for candidate in sorted(reference_units, key=_reference_rank, reverse=True):
        if any(
            _same_reference_unit(candidate, selected)
            for selected in reference_units
            if id(selected) in selected_reference_ids
        ):
            continue
        selected_reference_ids.add(id(candidate))

    kept: list[MGEFeature] = []
    curated_components = [
        feature
        for feature in features
        if feature.attributes.get("source") == "curated_reference"
    ]
    for f in features:
        if (
            f.element_type == "transposon"
            and f.attributes.get("source") == "reference_scan"
            and id(f) not in selected_reference_ids
        ):
            continue
        if (
            f.attributes.get("source") == "reference_scan"
            and any(
                f.element_type == curated.element_type
                and (
                    f.name == curated.name
                    or (
                        f.attributes.get("fid")
                        and f.attributes.get("fid") == curated.attributes.get("fid")
                    )
                )
                and _overlap(f, curated)
                for curated in curated_components
            )
        ):
            continue
        if (
            f.element_type == "IS"
            and f.attributes.get("source") not in {"reference_scan", "curated_reference"}
            and any(_overlap(f, exact) for exact in exact_reference_is)
        ):
            continue
        if (
            f.element_type == "transposon"
            and f.attributes.get("source") != "reference_scan"
            and f.family in blast_by_family
            and any(_overlap(f, b) for b in blast_by_family[f.family])
        ):
            continue
        kept.append(f)
    return kept


def _annotate_contig(
    seq: str,
    contig_id: str,
    features: list[MGEFeature],
    skip_boundaries: bool,
) -> tuple[list[MGEFeature], list[MGEFeature]]:
    """Run inference + boundary confirmation; return (roots, all_features)."""
    base_features = _flatten_feature_inputs(features)
    inferred = infer_transposons(base_features) + infer_integron_structures(base_features)
    all_feats = base_features + inferred
    all_feats = _suppress_redundant_inference(all_feats)

    # Tn1/Tn2/Tn3 parents are accepted from independently detected IR, gene,
    # res and cargo components. Whole-locus homology remains corroborating
    # naming evidence rather than the source of their internal annotation.
    all_feats.extend(assemble_tn123_components(all_feats))

    # Projection is now a labelled fallback only for components that the
    # independent scan could not recover.
    all_feats.extend(annotate_tn123(all_feats))
    all_feats.extend(annotate_partridge_units(all_feats))
    all_feats = _suppress_redundant_inference(all_feats)

    # Res sites are positioned from the surviving (post-dedup) transposons
    transposons = [f for f in all_feats if f.element_type == "transposon"]
    all_feats.extend(annotate_res_sites(transposons))

    if not skip_boundaries:
        checkable = [f for f in all_feats
                     if f.element_type in ("IS", "transposon")]
        confirm_boundaries(seq, checkable)

    assign_confidence(all_feats)
    roots = build_hierarchy(all_feats)
    return roots, all_feats


def _load_detector_outputs(
    isescan: str | None,
    amrfinder: str | None,
    integrons: str | None,
) -> list[MGEFeature]:
    features: list[MGEFeature] = []
    if isescan:
        features.extend(parse_isescan(isescan))
    if amrfinder:
        features.extend(parse_amrfinder(amrfinder))
    if integrons:
        features.extend(parse_integron_finder(integrons))
    return features


@cli.command(name="run")
@click.argument("fasta", type=click.Path(exists=True, path_type=Path))
@click.option(
    "--out",
    "-o",
    type=click.Path(path_type=Path),
    required=True,
    help="Results directory.",
)
@click.option(
    "--profile",
    type=click.Choice(sorted(REFERENCE_PROFILES)),
    default="validated",
    show_default=True,
    help="Validated Sally-backed references, or every experimental/legacy reference.",
)
@click.option(
    "--threads",
    type=click.IntRange(min=1),
    default=min(4, os.cpu_count() or 1),
    show_default=True,
    help="Concurrent reference scans.",
)
@click.option(
    "--detectors",
    type=click.Choice(["none", "available", "all"]),
    default="none",
    show_default=True,
    help="Optionally run ISEScan, AMRFinder+ and IntegronFinder from PATH/Pixi.",
)
@click.option("--isescan", type=click.Path(exists=True), help="Precomputed ISEScan TSV.")
@click.option("--amrfinder", type=click.Path(exists=True), help="Precomputed AMRFinder+ TSV.")
@click.option("--integrons", type=click.Path(exists=True), help="Precomputed IntegronFinder file.")
@click.option(
    "--mara-flank",
    type=click.IntRange(min=0),
    default=5_000,
    show_default=True,
    help="Flanking bases included around each MARA locus.",
)
def run_workflow(
    fasta: Path,
    out: Path,
    profile: str,
    threads: int,
    detectors: str,
    isescan: str | None,
    amrfinder: str | None,
    integrons: str | None,
    mara_flank: int,
) -> None:
    """Run a reproducible annotation and write a complete result directory."""
    if not blast_available():
        raise click.UsageError(
            "blastn and makeblastdb are required. Install BLAST+ or use the supplied Pixi environment."
        )
    out.mkdir(parents=True, exist_ok=True)
    detector_dir = out / "detectors"
    detector_paths: dict[str, str | None] = {
        "isescan": isescan,
        "amrfinder": amrfinder,
        "integron": integrons,
    }
    runners = {
        "isescan": run_isescan,
        "amrfinder": run_amrfinder,
        "integron": run_integron_finder,
    }
    if detectors != "none":
        detector_dir.mkdir(parents=True, exist_ok=True)
        for name, runner in runners.items():
            if detector_paths[name]:
                continue
            if not detector_available(name):
                if detectors == "all":
                    raise click.UsageError(
                        f"--detectors all requested, but {name} is not installed."
                    )
                click.echo(f"Detector unavailable; skipping {name}.", err=True)
                continue
            click.echo(f"Running {name}…", err=True)
            detector_paths[name] = str(runner(fasta, detector_dir / name))

    all_features = _load_detector_outputs(
        detector_paths["isescan"],
        detector_paths["amrfinder"],
        detector_paths["integron"],
    )
    click.echo(
        f"Reference profile: {profile}; scanning with {threads} worker(s)…",
        err=True,
    )

    def _progress(reference: str, completed: int, total: int) -> None:
        click.echo(f"  [{completed}/{total}] {reference}", err=True)

    all_features.extend(
        scan_all(fasta, profile=profile, threads=threads, progress=_progress)
    )
    fasta_records = list(SeqIO.parse(fasta, "fasta"))
    if not fasta_records:
        raise click.UsageError(f"No FASTA records found in {fasta}")

    annotated: list[tuple[str, int, list[MGEFeature]]] = []
    mara_dir = out / "mara"
    table_dir = out / "mara-table"
    hierarchy_dir = out / "hierarchy"
    mara_dir.mkdir(exist_ok=True)
    table_dir.mkdir(exist_ok=True)
    hierarchy_dir.mkdir(exist_ok=True)
    proof_output_paths: dict[tuple[str, int, int, str], dict[str, str]] = {}
    mara_locus_outputs: list[dict[str, object]] = []
    locus_count = 0
    for record in fasta_records:
        sequence = str(record.seq)
        contig_features = _features_on_contig(all_features, record.id)
        roots, _ = _annotate_contig(sequence, record.id, contig_features, False)
        annotated.append((record.id, len(sequence), roots))
        safe_record_id = record.id.replace("/", "_").replace(" ", "_")
        hierarchy_path = hierarchy_dir / f"{safe_record_id}.svg"
        hierarchy_path.write_text(
            to_linear_svg(roots, len(sequence), record.id),
            encoding="utf-8",
        )
        for locus in extract_mara_loci(roots, len(sequence), mara_flank):
            locus_name = f"{record.id}__{locus.suffix}"
            label = (
                f"{record.id} • {locus.target.name} "
                f"{locus.target.start}..{locus.target.end} • "
                f"view {locus.view_start}..{locus.view_end}"
            )
            (mara_dir / f"{locus_name}.svg").write_text(
                to_mara_svg(locus.roots, locus.view_length, label),
                encoding="utf-8",
            )
            (table_dir / f"{locus_name}.svg").write_text(
                to_mara_table_svg(locus.roots, label),
                encoding="utf-8",
            )
            proof_output_paths[(
                record.id,
                locus.target.start,
                locus.target.end,
                locus.target.name,
            )] = {
                "mara": f"mara/{locus_name}.svg",
                "mara_table": f"mara-table/{locus_name}.svg",
                "hierarchy": f"hierarchy/{safe_record_id}.svg",
                "cell_format": "annotation.cell",
                "annotation_json": "annotation.json",
                "annotation_gff3": "annotation.gff3",
            }
            mara_locus_outputs.append({
                "record": record.id,
                "call": locus.target.name,
                "family": locus.target.family,
                "element_type": locus.target.element_type,
                "start": locus.target.start,
                "end": locus.target.end,
                "view_start": locus.view_start,
                "view_end": locus.view_end,
                "mara": f"mara/{locus_name}.svg",
                "mara_table": f"mara-table/{locus_name}.svg",
                "hierarchy": f"hierarchy/{safe_record_id}.svg",
            })
            locus_count += 1

    detector_provenance = {
        name: str(Path(path).resolve())
        for name, path in detector_paths.items()
        if path
    }
    document = annotation_document(
        annotated,
        input_path=fasta,
        profile=profile,
        command=[
            "matryoshka",
            "run",
            str(fasta),
            "--profile",
            profile,
            "--threads",
            str(threads),
            "--out",
            str(out),
        ],
        detector_outputs=detector_provenance,
    )
    (out / "annotation.json").write_text(
        json.dumps(document, indent=2) + "\n",
        encoding="utf-8",
    )
    (out / "annotation.gff3").write_text(
        annotation_gff3(annotated),
        encoding="utf-8",
    )
    (out / "annotation.cell").write_text(
        "\n".join(
            to_wolvercote(roots, [], seqid)
            for seqid, _, roots in annotated
        ) + "\n",
        encoding="utf-8",
    )
    proof = build_tn123_proof(annotated, proof_output_paths)
    write_proof_bundle(
        out / "proof",
        proof,
        title=f"{fasta.name} Tn1/Tn2/Tn3 proof",
    )
    summary = {
        "schema_version": "1.0",
        "profile": profile,
        "records": len(annotated),
        "features": sum(count_features(roots) for _, _, roots in annotated),
        "mara_loci": locus_count,
        "mara_locus_outputs": mara_locus_outputs,
        "proof_status": proof["summary"]["status"],
        "outputs": {
            "annotation_json": "annotation.json",
            "annotation_gff3": "annotation.gff3",
            "annotation_cell": "annotation.cell",
            "hierarchy_directory": "hierarchy",
            "mara_directory": "mara",
            "mara_table_directory": "mara-table",
            "proof_json": "proof/proof.json",
            "proof_report": "proof/report.html",
            "component_ledger": "proof/components.tsv",
            "match_ledger": "proof/matches.tsv",
        },
    }
    (out / "run.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    click.echo(
        f"Complete: {summary['features']} features and {locus_count} MARA loci in {out}",
        err=True,
    )


def _render(
    fmt: str,
    roots: list[MGEFeature],
    seq: str,
    sample_name: str,
) -> str | bytes:
    if fmt == "json":
        return to_json(roots)
    if fmt == "gff3":
        return to_gff3(roots, seqid=sample_name)
    if fmt == "genbank":
        return to_genbank(roots, seq, sample_name)
    if fmt == "wolvercote":
        return to_wolvercote(roots, [], sample_name)
    if fmt == "linear":
        return to_linear_svg(roots, len(seq), sample_name)
    if fmt == "mara":
        return to_mara_svg(roots, len(seq), sample_name)
    if fmt == "mara-table":
        return to_mara_table_svg(roots, sample_name)
    raise click.BadParameter(f"unknown format: {fmt}")


@cli.command()
@click.argument("fasta", type=click.Path(exists=True))
@click.option("--isescan",   type=click.Path(exists=True), help="ISEScan .tsv output")
@click.option("--amrfinder", type=click.Path(exists=True), help="AMRFinder+ .tsv output")
@click.option("--integrons", type=click.Path(exists=True), help="IntegronFinder .integrons file")
@click.option(
    "--format", "fmt", default="json",
    type=click.Choice([
        "json", "gff3", "genbank", "wolvercote", "linear", "mara",
        "mara-table",
    ]),
    show_default=True, help="Output format",
)
@click.option("--no-boundaries", is_flag=True, help="Skip TSD/IR confirmation step")
@click.option("--reference-scan/--no-reference-scan", default=True,
              help="Run BLAST scan against bundled MGE references (Tn1546, Tn4401 variants, "
                   "Acinetobacter islands, etc.). Requires blastn on PATH.")
@click.option(
    "--profile",
    type=click.Choice(sorted(REFERENCE_PROFILES)),
    default="validated",
    show_default=True,
    help="Reference set used by the BLAST scan.",
)
@click.option(
    "--threads",
    type=click.IntRange(min=1),
    default=min(4, os.cpu_count() or 1),
    show_default=True,
    help="Concurrent reference scans.",
)
@click.option(
    "--mara-flank",
    type=click.IntRange(min=0),
    default=5_000,
    show_default=True,
    help="Flanking bases included around each validated MARA target locus.",
)
@click.option("--out", "-o", default="-", help="Output file or directory (default: stdout)")
def annotate(
    fasta: str,
    isescan: str | None,
    amrfinder: str | None,
    integrons: str | None,
    fmt: str,
    no_boundaries: bool,
    reference_scan: bool,
    profile: str,
    threads: int,
    mara_flank: int,
    out: str,
) -> None:
    """Combine detection tool outputs into a nested MGE annotation.

    For multi-contig FASTAs each contig is processed independently. JSON,
    GFF3 and Wolvercote outputs are concatenated; SVG outputs are
    written per-contig to the directory given by --out (required for
    those formats with multi-contig input).
    """
    has_detector_output = any([isescan, amrfinder, integrons])
    if not has_detector_output and not reference_scan:
        raise click.UsageError(
            "Provide at least one detection tool output or enable the reference scan."
        )

    all_features = _load_detector_outputs(isescan, amrfinder, integrons)

    records = list(SeqIO.parse(fasta, "fasta"))
    if not records:
        raise click.UsageError(f"No FASTA records found in {fasta}")

    click.echo(
        f"Loaded {len(all_features)} features across {len(records)} contig(s)",
        err=True,
    )

    if reference_scan and blast_available():
        click.echo(f"Reference profile: {profile}; {threads} worker(s)", err=True)
        ref_hits = scan_all(fasta, profile=profile, threads=threads)
        if ref_hits:
            click.echo(
                f"BLAST reference-scan: {len(ref_hits)} hit(s)", err=True,
            )
            all_features.extend(ref_hits)
    elif reference_scan:
        if not has_detector_output:
            raise click.UsageError(
                "Reference-only annotation requires blastn on PATH."
            )
        click.echo("blastn not on PATH — skipping reference scan.", err=True)

    multi = len(records) > 1
    is_binary = False
    outputs: list[tuple[str, str | bytes]] = []

    for rec in records:
        contig_id = rec.id
        seq = str(rec.seq)
        contig_feats = _features_on_contig(all_features, contig_id) if multi else all_features
        roots, _ = _annotate_contig(seq, contig_id, contig_feats, no_boundaries)
        click.echo(
            f"  {contig_id} ({len(seq):,} bp): "
            f"{len(roots)} root-level elements from {len(contig_feats)} features",
            err=True,
        )
        if fmt in {"mara", "mara-table"}:
            loci = extract_mara_loci(roots, len(seq), mara_flank)
            if loci:
                for locus in loci:
                    locus_name = f"{contig_id}__{locus.suffix}"
                    label = (
                        f"{contig_id} • {locus.target.name} "
                        f"{locus.target.start}..{locus.target.end} • "
                        f"view {locus.view_start}..{locus.view_end}"
                    )
                    rendered: str | bytes = (
                        to_mara_svg(locus.roots, locus.view_length, label)
                        if fmt == "mara"
                        else to_mara_table_svg(locus.roots, label)
                    )
                    outputs.append((locus_name, rendered))
                continue
            if multi:
                # A multi-FASTA may contain records without a validated target
                # (for example the standalone ISEcp1 control beside TPU records).
                # Do not create a misleading blank per-locus diagram for them.
                continue
        rendered = _render(fmt, roots, seq, contig_id)
        outputs.append((contig_id, rendered))

    directory_requested = (
        fmt in {"mara", "mara-table"}
        and out != "-"
        and Path(out).suffix.lower() != ".svg"
    )
    _write_outputs(
        outputs,
        out,
        fmt,
        is_binary,
        multi or len(outputs) > 1 or directory_requested,
    )


def _write_outputs(
    outputs: list[tuple[str, str | bytes]],
    out: str,
    fmt: str,
    is_binary: bool,
    multi: bool,
) -> None:
    # Single-contig: write the one output to `out`
    if not multi:
        _, data = outputs[0]
        _emit(data, out, is_binary)
        return

    # Multi-contig with aggregatable text formats: concatenate to one stream
    aggregatable = fmt in {"json", "gff3", "wolvercote", "genbank"}
    if aggregatable:
        if fmt == "json":
            import json as _json
            combined = "{\n" + ",\n".join(
                f'  {_json.dumps(cid)}: '
                f'{data.decode("utf-8") if isinstance(data, bytes) else data}'
                for cid, data in outputs
            ) + "\n}"
        else:
            combined = "\n".join(str(data) for _, data in outputs)
        _emit(combined, out, False)
        return

    # Per-contig visual formats require an output directory.
    if out == "-":
        raise click.UsageError(
            f"--format {fmt} on multi-contig input requires --out <dir>"
        )
    out_dir = Path(out)
    out_dir.mkdir(parents=True, exist_ok=True)
    ext = "png" if fmt == "png" else "svg"
    for cid, data in outputs:
        safe = cid.replace("/", "_").replace(" ", "_")
        path = out_dir / f"{safe}.{ext}"
        if is_binary:
            path.write_bytes(data)  # type: ignore[arg-type]
        else:
            path.write_text(str(data))
    click.echo(f"Wrote {len(outputs)} files to {out_dir}", err=True)


def _emit(data: str | bytes, out: str, is_binary: bool) -> None:
    if out == "-":
        if is_binary:
            sys.stdout.buffer.write(data)  # type: ignore[arg-type]
        else:
            click.echo(data)
    else:
        path = Path(out)
        if is_binary:
            path.write_bytes(data)  # type: ignore[arg-type]
        else:
            path.write_text(str(data))
        click.echo(f"Written to {out}", err=True)


if __name__ == "__main__":
    cli()
