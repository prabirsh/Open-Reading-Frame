"""Numbers about sequences, proteins and ORF lengths."""

import math

# Average residue masses in Da (the values Expasy uses for average mass).
_RESIDUE_MASS = {
    "G": 57.0519, "A": 71.0788, "S": 87.0782, "P": 97.1167, "V": 99.1326,
    "T": 101.1051, "C": 103.1388, "L": 113.1594, "I": 113.1594, "N": 114.1038,
    "D": 115.0886, "Q": 128.1307, "K": 128.1741, "E": 129.1155, "M": 131.1926,
    "H": 137.1411, "F": 147.1766, "R": 156.1875, "Y": 163.1760, "W": 186.2132,
}
_WATER = 18.01524
_TYPICAL_RESIDUE = 110.0   # stand-in for X, which has no defined mass

AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"
AMINO_ACID_NAMES = {
    "A": "alanine", "C": "cysteine", "D": "aspartate", "E": "glutamate",
    "F": "phenylalanine", "G": "glycine", "H": "histidine", "I": "isoleucine",
    "K": "lysine", "L": "leucine", "M": "methionine", "N": "asparagine",
    "P": "proline", "Q": "glutamine", "R": "arginine", "S": "serine",
    "T": "threonine", "V": "valine", "W": "tryptophan", "Y": "tyrosine",
}

_COMPLEMENT_BASE = {"A": "T", "T": "A", "C": "G", "G": "C"}


def gc_percent(seq):
    """GC content as a percentage of unambiguous bases (A, C, G, T only)."""
    gc = seq.count("G") + seq.count("C")
    total = gc + seq.count("A") + seq.count("T")
    return 100.0 * gc / total if total else 0.0


def auto_window(length):
    """A GC window that gives roughly 120 points along the sequence."""
    return max(50, int(round(length / 120, -1)))


def gc_windows(seq, window):
    """Non-overlapping windows as (midpoint, gc%) pairs, 1-based midpoints."""
    n = len(seq)
    if n == 0:
        return []
    window = max(1, min(int(window), n))
    points = []
    for a in range(0, n, window):
        b = min(a + window, n)
        if b - a < window / 2 and points:
            break   # a short tail would only add noise
        chunk = seq[a:b]
        if any(chunk.count(base) for base in "ACGT"):
            points.append(((a + b + 1) / 2, gc_percent(chunk)))
    return points


def protein_mass_kda(protein):
    residues = sum(_RESIDUE_MASS.get(aa, _TYPICAL_RESIDUE) for aa in protein)
    return (residues + _WATER) / 1000.0


def residue_counts(proteins):
    counts = dict.fromkeys(AMINO_ACIDS, 0)
    for protein in proteins:
        for aa in AMINO_ACIDS:
            counts[aa] += protein.count(aa)
    return counts


# --- how surprising is an ORF? -------------------------------------------
#
# In random sequence a stop codon turns up every 1/p codons, so the stretch
# between two stops is geometric: it reaches L codons with probability
# (1 - p)^L. Multiply by the number of stops we expect to see in the search
# space and you get the number of stop-free stretches of that length, an upper
# bound on how many ORFs of that size chance alone could produce.

def base_frequencies(seq):
    counts = {b: seq.count(b) + 1 for b in "ACGT"}   # +1 keeps tiny inputs sane
    total = sum(counts.values())
    return {b: counts[b] / total for b in "ACGT"}


def flip_strand(freqs):
    return {b: freqs[_COMPLEMENT_BASE[b]] for b in "ACGT"}


def stop_probability(freqs, stops):
    total = 0.0
    for codon in stops:
        p = 1.0
        for base in codon:
            p *= freqs[base]
        total += p
    return total


def expected_by_chance(length, aa_length, stop_prob, frame_count):
    """Expected number of ORFs this long (or longer) on one strand by chance."""
    return length * (frame_count / 3.0) * stop_prob * (1.0 - stop_prob) ** aa_length


def support_label(expect):
    if expect < 1e-3:
        return "strong"
    if expect < 0.05:
        return "moderate"
    return "weak"


# --- histogram of ORF lengths ----------------------------------------------

def _nice_int_step(raw):
    raw = max(raw, 1.0)
    magnitude = 10 ** math.floor(math.log10(raw))
    for m in (1, 2, 5, 10):
        if raw <= m * magnitude:
            return int(m * magnitude)
    return int(10 * magnitude)


def length_histogram(forward, reverse, target_bins=12):
    """Bin ORF lengths (in aa) by strand.

    Returns a list of dicts: {lo, hi, label, forward, reverse}. When a long
    tail would squash everything else, it is folded into one final bin.
    """
    values = sorted(forward + reverse)
    if not values:
        return []
    lo, hi = values[0], values[-1]
    cap = hi
    folded = False
    p95 = values[int(0.95 * (len(values) - 1))]
    if len(values) >= 20 and hi > 2 * p95:
        cap, folded = p95, True
    step = _nice_int_step((cap - lo + 1) / target_bins)
    first = (lo // step) * step
    count = max(1, math.ceil((cap - first + 1) / step))
    bins = [{"lo": first + i * step, "hi": first + (i + 1) * step - 1,
             "forward": 0, "reverse": 0} for i in range(count)]
    if folded:
        bins.append({"lo": first + count * step, "hi": None, "forward": 0, "reverse": 0})
    for key, lengths in (("forward", forward), ("reverse", reverse)):
        for value in lengths:
            index = min((value - first) // step, len(bins) - 1)
            bins[index][key] += 1
    for b in bins:
        b["label"] = f"≥{b['lo']}" if b["hi"] is None else str(b["lo"])
    return bins


def median(values):
    ordered = sorted(values)
    if not ordered:
        return 0
    mid = len(ordered) // 2
    if len(ordered) % 2:
        return ordered[mid]
    return (ordered[mid - 1] + ordered[mid]) / 2
