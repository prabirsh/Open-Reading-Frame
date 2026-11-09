"""Genetic codes.

Tables are stored the way NCBI publishes them: 64 amino-acid letters, one per
codon, with the codons listed in TCAG order (TTT, TTC, TTA, TTG, TCT, ...).
A '*' marks a stop codon.
"""

from itertools import product

_CODON_ORDER = ["".join(c) for c in product("TCAG", repeat=3)]


class GeneticCode:
    def __init__(self, table_id, name, amino_acids, starts):
        if len(amino_acids) != 64:
            raise ValueError("a genetic code needs 64 amino-acid letters")
        self.id = table_id
        self.name = name
        self.table = dict(zip(_CODON_ORDER, amino_acids))
        self.stops = frozenset(c for c, aa in self.table.items() if aa == "*")
        self.starts = tuple(starts)

    def __repr__(self):
        return f"GeneticCode({self.id}, {self.name!r})"

    def translate(self, seq, offset=0, to_stop=False):
        """Translate whole codons of `seq` from `offset`. Unknown codons become X."""
        protein = []
        for i in range(offset, len(seq) - 2, 3):
            aa = self.table.get(seq[i:i + 3], "X")
            if aa == "*" and to_stop:
                break
            protein.append(aa)
        return "".join(protein)


_STANDARD = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
_VERT_MITO = "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG"

CODES = {
    1: GeneticCode(1, "Standard", _STANDARD, ("ATG", "TTG", "CTG")),
    2: GeneticCode(2, "Vertebrate mitochondrial", _VERT_MITO,
                   ("ATT", "ATC", "ATA", "ATG", "GTG")),
    11: GeneticCode(11, "Bacterial, archaeal and plant plastid", _STANDARD,
                    ("TTG", "CTG", "ATT", "ATC", "ATA", "ATG", "GTG")),
}


def get_code(table_id):
    try:
        return CODES[int(table_id)]
    except (KeyError, ValueError):
        known = ", ".join(str(i) for i in CODES)
        raise ValueError(f"unknown genetic code {table_id!r} (available: {known})") from None
