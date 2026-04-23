"""
Matryoshka CLI — detect and annotate nested MGEs in bacterial genomes.

Usage:
    matryoshka annotate <fasta> [--isescan <tsv>] [--amrfinder <tsv>]
                                [--integrons <file>]
                                [--format gff3|json|wolvercote|svg|png|linear]
                                [--out <path>]
"""

from __future__ import annotations

import sys
from importlib.metadata import version as _pkg_version, PackageNotFoundError
from pathlib import Path

import click
from Bio import SeqIO

from .boundaries import confirm_boundaries
from .confidence import assign_confidence
from .detect import (
    MGEFeature,
    parse_amrfinder,
    parse_integron_finder,
    parse_isescan,
)
from .hierarchy import build_hierarchy
from .output import to_genbank, to_gff3, to_json, to_png, to_svg, to_wolvercote
from .reference_scan import blast_available, scan_all
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


def _suppress_redundant_inference(features: list[MGEFeature]) -> list[MGEFeature]:
    """Drop rule-inferred transposons that overlap a BLAST-confirmed call of
    the same family. The BLAST boundaries are authoritative.
    """
    blast_by_family: dict[str, list[MGEFeature]] = {}
    for f in features:
        if f.attributes.get("source") == "reference_scan" and f.element_type == "transposon":
            blast_by_family.setdefault(f.family, []).append(f)

    def _overlap(a: MGEFeature, b: MGEFeature) -> bool:
        if a.end < b.start or b.end < a.start:
            return False
        ovl = min(a.end, b.end) - max(a.start, b.start)
        short = min(a.end - a.start, b.end - b.start)
        return short > 0 and ovl / short > 0.5

    kept: list[MGEFeature] = []
    for f in features:
        if (
            f.element_type == "transposon"
            and f.attributes.get("source") != "reference_scan"
            and f.family in blast_by_family
        ):
            if any(_overlap(f, b) for b in blast_by_family[f.family]):
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
    inferred = infer_transposons(features)
    all_feats = features + inferred
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
    if fmt == "svg":
        return to_svg(roots, sample_name)
    if fmt == "png":
        return to_png(roots, sample_name)
    if fmt == "linear":
        return to_linear_svg(roots, len(seq), sample_name)
    raise click.BadParameter(f"unknown format: {fmt}")


@cli.command()
@click.argument("fasta", type=click.Path(exists=True))
@click.option("--isescan",   type=click.Path(exists=True), help="ISEScan .tsv output")
@click.option("--amrfinder", type=click.Path(exists=True), help="AMRFinder+ .tsv output")
@click.option("--integrons", type=click.Path(exists=True), help="IntegronFinder .integrons file")
@click.option(
    "--format", "fmt", default="json",
    type=click.Choice(["json", "gff3", "genbank", "wolvercote", "svg", "png", "linear"]),
    show_default=True, help="Output format",
)
@click.option("--no-boundaries", is_flag=True, help="Skip TSD/IR confirmation step")
@click.option("--reference-scan/--no-reference-scan", default=True,
              help="Run BLAST scan against bundled MGE references (Tn1546, Tn4401 variants, "
                   "Acinetobacter islands, etc.). Requires blastn on PATH.")
@click.option("--out", "-o", default="-", help="Output file or directory (default: stdout)")
def annotate(
    fasta: str,
    isescan: str | None,
    amrfinder: str | None,
    integrons: str | None,
    fmt: str,
    no_boundaries: bool,
    reference_scan: bool,
    out: str,
) -> None:
    """Combine detection tool outputs into a nested MGE annotation.

    For multi-contig FASTAs each contig is processed independently. JSON,
    GFF3 and Wolvercote outputs are concatenated; SVG/PNG outputs are
    written per-contig to the directory given by --out (required for
    those formats with multi-contig input).
    """
    if not any([isescan, amrfinder, integrons]):
        raise click.UsageError(
            "Provide at least one detection tool output "
            "(--isescan, --amrfinder, --integrons)."
        )

    all_features: list[MGEFeature] = []
    if isescan:
        all_features.extend(parse_isescan(isescan))
    if amrfinder:
        all_features.extend(parse_amrfinder(amrfinder))
    if integrons:
        all_features.extend(parse_integron_finder(integrons))

    records = list(SeqIO.parse(fasta, "fasta"))
    if not records:
        raise click.UsageError(f"No FASTA records found in {fasta}")

    click.echo(
        f"Loaded {len(all_features)} features across {len(records)} contig(s)",
        err=True,
    )

    if reference_scan and blast_available():
        ref_hits = scan_all(fasta)
        if ref_hits:
            click.echo(
                f"BLAST reference-scan: {len(ref_hits)} hit(s)", err=True,
            )
            all_features.extend(ref_hits)
    elif reference_scan:
        click.echo("blastn not on PATH — skipping reference scan.", err=True)

    multi = len(records) > 1
    is_binary = fmt == "png"
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
        rendered = _render(fmt, roots, seq, contig_id)
        outputs.append((contig_id, rendered))

    _write_outputs(outputs, out, fmt, is_binary, multi)


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
                f'  {_json.dumps(cid)}: {data}' for cid, data in outputs
            ) + "\n}"
        else:
            combined = "\n".join(str(data) for _, data in outputs)
        _emit(combined, out, False)
        return

    # Per-contig formats (svg/png/linear): require --out directory
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
