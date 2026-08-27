#!/usr/bin/env python3
"""Fetch and run the fixed real-accession Tn1 evaluation panel."""

from __future__ import annotations

import argparse
import subprocess
import urllib.parse
import urllib.request
from pathlib import Path

ACCESSIONS = (
    "PP591960.1", "PP591959.1", "AM261760.1", "NZ_LK391770.1",
    "AP018939.2", "MN148427.1", "MG764534.2", "KX156773.1",
    "KX156772.1", "HQ451074.1", "AF550679.1", "AY150843.2",
    "GQ160960.2", "HM804085.1", "AY925200.1", "GQ293500.1",
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=Path("results/tn1-accessions"))
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    query = urllib.parse.urlencode({
        "db": "nuccore",
        "id": ",".join(ACCESSIONS),
        "rettype": "fasta",
        "retmode": "text",
    })
    request = urllib.request.Request(
        f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?{query}",
        headers={"User-Agent": "matryoshka-accession-evaluation/1.0"},
    )
    fasta = args.out / "accessions.fasta"
    with urllib.request.urlopen(request, timeout=120) as response:
        fasta.write_bytes(response.read())
    subprocess.run([
        "matryoshka", "run", str(fasta),
        "--profile", "component-rules",
        "--detectors", "none",
        "--out", str(args.out / "component-rule-results"),
    ], check=True)


if __name__ == "__main__":
    main()
