"""Regenerate orf/data/demo.fasta.

The demo sequences are synthetic: random DNA with a handful of "genes" planted
on both strands, so the report has something worth looking at. Run it from the
repository root:

    python tools/make_demo.py
"""

import random
from pathlib import Path

STOPS = ("TAA", "TAG", "TGA")
COMPLEMENT = str.maketrans("ACGT", "TGCA")
SEED = 2026


def random_dna(rng, length, gc):
    at, cg = (1 - gc) / 2, gc / 2
    return "".join(rng.choices("ATGC", [at, at, cg, cg], k=length))


def random_gene(rng, aa, gc):
    """ATG, then aa-1 random sense codons, then a stop codon."""
    codons = ["ATG"]
    while len(codons) < aa:
        codon = random_dna(rng, 3, gc)
        if codon not in STOPS:
            codons.append(codon)
    codons.append(rng.choice(STOPS))
    return "".join(codons)


def reverse_complement(seq):
    return seq.translate(COMPLEMENT)[::-1]


def build(rng, layout):
    """layout is a list of (kind, size, gc); kind is gap, fwd or rev."""
    parts = []
    for kind, size, gc in layout:
        if kind == "gap":
            parts.append(random_dna(rng, size, gc))
        elif kind == "fwd":
            parts.append(random_gene(rng, size, gc))
        else:
            parts.append(reverse_complement(random_gene(rng, size, gc)))
    return "".join(parts)


def wrap(seq, width=70):
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))


def main():
    rng = random.Random(SEED)
    operon = build(rng, [
        ("gap", 150, 0.52), ("fwd", 420, 0.56), ("gap", 80, 0.50),
        ("fwd", 310, 0.56), ("gap", 60, 0.50), ("rev", 520, 0.57),
        ("gap", 120, 0.50), ("fwd", 260, 0.55), ("gap", 300, 0.31),
        ("rev", 340, 0.54), ("gap", 100, 0.50), ("fwd", 180, 0.55),
        ("gap", 80, 0.50), ("rev", 220, 0.54), ("gap", 150, 0.52),
    ])
    island = build(rng, [
        ("gap", 120, 0.38), ("fwd", 140, 0.40), ("gap", 90, 0.35),
        ("rev", 110, 0.40), ("gap", 100, 0.36),
    ])
    out = Path(__file__).resolve().parent.parent / "orf" / "data" / "demo.fasta"
    out.parent.mkdir(exist_ok=True)
    out.write_text(
        f">demo_operon synthetic: 7 planted genes on both strands, with an AT-rich insert (seed {SEED})\n"
        f"{wrap(operon)}\n"
        f">demo_island synthetic: 2 planted genes in AT-rich sequence (seed {SEED})\n"
        f"{wrap(island)}\n",
        encoding="utf-8", newline="\n")
    print(f"wrote {out} ({len(operon)} + {len(island)} nt)")


if __name__ == "__main__":
    main()
