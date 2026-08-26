#!/usr/bin/env python3
"""Create a deterministic arbitrary-contig Tn1/Tn2/Tn3 demonstration FASTA."""

from __future__ import annotations

import argparse
import random
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq


def random_dna(length: int, seed: int) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(length))


def mutate(sequence: str, count: int, seed: int) -> str:
    bases = list(sequence)
    rng = random.Random(seed)
    replacement = {"A": "C", "C": "G", "G": "T", "T": "A"}
    for position in rng.sample(range(50, len(bases) - 50), count):
        bases[position] = replacement[bases[position]]
    return "".join(bases)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path)
    args = parser.parse_args()

    references = (
        Path(__file__).parents[1]
        / "matryoshka"
        / "references"
        / "tn1_tn2_tn3.fasta"
    )
    sequences = {
        record.id: str(record.seq)
        for record in SeqIO.parse(references, "fasta")
    }
    tn1_like = mutate(sequences["Tn1_NC_008357"], 25, 71)
    tn2 = sequences["Tn2_AY123253"]
    tn2_inserted = tn2[:2500] + random_dna(800, 72) + tn2[2500:]
    reverse_tn3 = str(Seq(sequences["Tn3_HM749966"]).reverse_complement())
    contig = "".join([
        random_dna(8_000, 73),
        tn1_like,
        random_dna(7_000, 74),
        tn2_inserted,
        random_dna(6_000, 75),
        reverse_tn3,
        random_dna(8_000, 76),
    ])
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(f">arbitrary_demo_contig\n{contig}\n")


if __name__ == "__main__":
    main()
