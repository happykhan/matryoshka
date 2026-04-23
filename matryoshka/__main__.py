"""
Matryoshka CLI — detect and annotate nested MGEs in bacterial genomes.

Usage:
    matryoshka annotate <fasta> --isescan <tsv> --amrfinder <tsv>
                                --integrons <file> [--out-dir <dir>]
                                [--format gff3|json|wolvercote]
"""

from __future__ import annotations

import sys
from pathlib import Path

import click
from Bio import SeqIO

from .boundaries import confirm_boundaries
from .detect import (
    parse_amrfinder,
    parse_integron_finder,
    parse_isescan,
    parse_mobsuite,
)
from .hierarchy import build_hierarchy
from .output import to_gff3, to_json, to_png, to_svg, to_wolvercote
from .transposon import infer_transposons


@click.group()
def cli() -> None:
    pass


@cli.command()
@click.argument("fasta", type=click.Path(exists=True))
@click.option("--isescan",    type=click.Path(exists=True), help="ISEScan .tsv output")
@click.option("--amrfinder",  type=click.Path(exists=True), help="AMRFinder+ .tsv output")
@click.option("--integrons",  type=click.Path(exists=True), help="IntegronFinder .integrons file")
@click.option("--mobsuite",   type=click.Path(exists=True), help="MOBsuite mob_typer results (optional)")
@click.option("--format", "fmt", default="json",
              type=click.Choice(["json", "gff3", "wolvercote", "svg", "png"]),
              show_default=True, help="Output format")
@click.option("--no-boundaries", is_flag=True, help="Skip TSD/IR confirmation step")
@click.option("--out", "-o", default="-", help="Output file (default: stdout)")
def annotate(
    fasta: str,
    isescan: str | None,
    amrfinder: str | None,
    integrons: str | None,
    mobsuite: str | None,
    fmt: str,
    no_boundaries: bool,
    out: str,
) -> None:
    """Combine detection tool outputs into a nested MGE annotation."""
    if not any([isescan, amrfinder, integrons, mobsuite]):
        raise click.UsageError("Provide at least one detection tool output (--isescan, --amrfinder, --integrons).")

    seq_record = SeqIO.read(fasta, "fasta")
    seq = str(seq_record.seq)
    sample_name = seq_record.id

    all_features = []
    if isescan:
        all_features.extend(parse_isescan(isescan))
    if amrfinder:
        all_features.extend(parse_amrfinder(amrfinder))
    if integrons:
        all_features.extend(parse_integron_finder(integrons))
    if mobsuite:
        all_features.extend(parse_mobsuite(mobsuite))

    click.echo(f"Loaded {len(all_features)} features", err=True)

    inferred = infer_transposons(all_features)
    if inferred:
        click.echo(f"Inferred {len(inferred)} composite transposons", err=True)
    all_features = all_features + inferred

    if not no_boundaries:
        # IS elements + inferred transposons (Tn3 superfamily have 5bp TSDs)
        checkable = [f for f in all_features
                     if f.element_type in ("IS", "transposon")]
        confirm_boundaries(seq, checkable)
        confirmed = sum(1 for f in checkable if f.tsd_seq or f.ir_left)
        click.echo(
            f"Boundary evidence on {confirmed}/{len(checkable)} elements", err=True
        )

    roots = build_hierarchy(all_features)
    click.echo(f"Root-level elements: {len(roots)}", err=True)

    if fmt == "png":
        data = to_png(roots, sample_name)
        if out == "-":
            sys.stdout.buffer.write(data)
        else:
            Path(out).write_bytes(data)
            click.echo(f"Written to {out}", err=True)
        return

    if fmt == "json":
        output = to_json(roots)
    elif fmt == "gff3":
        output = to_gff3(roots, seqid=sample_name)
    elif fmt == "svg":
        output = to_svg(roots, sample_name)
    else:
        output = to_wolvercote(roots, [], sample_name)

    if out == "-":
        click.echo(output)
    else:
        Path(out).write_text(output)
        click.echo(f"Written to {out}", err=True)


if __name__ == "__main__":
    cli()
