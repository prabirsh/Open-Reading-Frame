"""The ORF scan itself."""

from dataclasses import dataclass

from .seqio import reverse_complement

STRANDS = ("both", "forward", "reverse")


@dataclass
class SearchParams:
    code: object
    strands: str = "both"
    frames: tuple = (1, 2, 3)
    min_aa: int = 30
    starts: tuple = ("ATG",)
    nested: bool = False
    partial: bool = False

    def __post_init__(self):
        if self.strands not in STRANDS:
            raise ValueError(f"strand must be one of: {', '.join(STRANDS)}")
        self.frames = tuple(sorted(set(self.frames)))
        if not self.frames or not set(self.frames) <= {1, 2, 3}:
            raise ValueError("frames must be chosen from 1, 2 and 3")
        if self.min_aa < 1:
            raise ValueError("minimum ORF length must be at least 1 amino acid")
        starts = tuple(dict.fromkeys(s.upper().replace("U", "T") for s in self.starts))
        for codon in starts:
            if len(codon) != 3 or set(codon) - set("ACGT"):
                raise ValueError(f"{codon!r} is not a valid start codon")
            if codon in self.code.stops:
                raise ValueError(f"{codon} is a stop codon in genetic code {self.code.id}")
        if not starts:
            raise ValueError("at least one start codon is required")
        self.starts = starts

    @property
    def strand_symbols(self):
        return {"both": "+-", "forward": "+", "reverse": "-"}[self.strands]


@dataclass
class Orf:
    seq_id: str
    strand: str
    frame: int
    start: int          # 1-based, inclusive, always on the forward strand
    end: int
    start_codon: str
    stop_codon: str     # empty when the ORF runs off the end of the sequence
    protein: str
    orf_id: str = ""
    gc: float = 0.0
    mass_kda: float = 0.0
    expect: float = 1.0
    support: str = ""

    @property
    def nt_length(self):
        return self.end - self.start + 1

    @property
    def aa_length(self):
        return len(self.protein)

    @property
    def partial(self):
        return not self.stop_codon

    @property
    def frame_label(self):
        return f"{self.strand}{self.frame}"

    def nucleotides(self, seq):
        piece = seq[self.start - 1:self.end]
        return piece if self.strand == "+" else reverse_complement(piece)


def _scan(seq, offset, params):
    """Yield (start, end, has_stop) for every ORF on one strand and frame.

    Positions are 0-based and `end` is exclusive, counted along `seq` itself.
    """
    if len(seq) - offset < 3:
        return
    stops = params.code.stops
    starts = set(params.starts)
    last = offset + 3 * ((len(seq) - offset) // 3)
    open_starts = []
    for i in range(offset, last, 3):
        codon = seq[i:i + 3]
        if codon in stops:
            for s in open_starts:
                if (i - s) // 3 >= params.min_aa:
                    yield s, i + 3, True
            open_starts = []
        elif codon in starts and (params.nested or not open_starts):
            open_starts.append(i)
    if params.partial:
        for s in open_starts:
            if (last - s) // 3 >= params.min_aa:
                yield s, last, False


def find_orfs(record, params):
    """All ORFs in one record, ordered by position and numbered."""
    seq = record.seq
    total = len(seq)
    found = []
    for strand in params.strand_symbols:
        strand_seq = seq if strand == "+" else reverse_complement(seq)
        for frame in params.frames:
            for s, e, stopped in _scan(strand_seq, frame - 1, params):
                body_end = e - 3 if stopped else e
                # whatever the start codon was, the initiator residue is Met
                protein = "M" + params.code.translate(strand_seq[s + 3:body_end])
                if strand == "+":
                    start, end = s + 1, e
                else:
                    start, end = total - e + 1, total - s
                found.append(Orf(
                    seq_id=record.id, strand=strand, frame=frame,
                    start=start, end=end,
                    start_codon=strand_seq[s:s + 3],
                    stop_codon=strand_seq[e - 3:e] if stopped else "",
                    protein=protein))
    found.sort(key=lambda o: (o.start, o.end, o.strand))
    for number, orf in enumerate(found, 1):
        orf.orf_id = f"{record.id}_ORF{number}"
    return found
