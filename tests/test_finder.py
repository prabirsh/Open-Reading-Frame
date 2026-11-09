import random
import unittest

from orf.codes import get_code
from orf.finder import SearchParams, find_orfs
from orf.seqio import Record, reverse_complement


def params(**overrides):
    base = dict(code=get_code(1), min_aa=1)
    base.update(overrides)
    return SearchParams(**base)


def find(seq, **overrides):
    return find_orfs(Record("s", "", seq), params(**overrides))


class BasicTests(unittest.TestCase):
    def test_single_forward_orf(self):
        (orf,) = find("ATGAAATAG")
        self.assertEqual((orf.strand, orf.frame, orf.start, orf.end), ("+", 1, 1, 9))
        self.assertEqual((orf.protein, orf.start_codon, orf.stop_codon), ("MK", "ATG", "TAG"))
        self.assertEqual((orf.aa_length, orf.nt_length, orf.partial), (2, 9, False))
        self.assertEqual(orf.orf_id, "s_ORF1")
        self.assertEqual(orf.frame_label, "+1")

    def test_frames_are_reported_from_one(self):
        (orf,) = find("CATGAAATAGC")
        self.assertEqual((orf.frame, orf.start, orf.end), (2, 2, 10))
        (orf,) = find("CCATGAAATAGC")
        self.assertEqual((orf.frame, orf.start, orf.end), (3, 3, 11))

    def test_reverse_strand_uses_forward_coordinates(self):
        seq = "GG" + reverse_complement("ATGAAATAG") + "CC"
        (orf,) = find(seq)
        self.assertEqual((orf.strand, orf.frame, orf.start, orf.end), ("-", 3, 3, 11))
        self.assertEqual(orf.protein, "MK")
        self.assertEqual(orf.nucleotides(seq), "ATGAAATAG")
        self.assertEqual(orf.frame_label, "-3")

    def test_strand_and_frame_selection(self):
        seq = "GG" + reverse_complement("ATGAAATAG") + "CC" + "ATGAAATAG"
        self.assertEqual({o.strand for o in find(seq)}, {"+", "-"})
        self.assertEqual({o.strand for o in find(seq, strands="forward")}, {"+"})
        self.assertEqual({o.strand for o in find(seq, strands="reverse")}, {"-"})
        self.assertEqual(find("CATGAAATAGC", frames=(1, 3)), [])
        self.assertEqual(len(find("CATGAAATAGC", frames=(2,))), 1)

    def test_minimum_length(self):
        self.assertEqual(len(find("ATGAAATAG", min_aa=2)), 1)
        self.assertEqual(find("ATGAAATAG", min_aa=3), [])

    def test_outermost_start_by_default_nested_on_request(self):
        (outer,) = find("ATGATGAAATAG")
        self.assertEqual((outer.start, outer.protein), (1, "MMK"))
        both = find("ATGATGAAATAG", nested=True)
        self.assertEqual([(o.start, o.protein) for o in both], [(1, "MMK"), (4, "MK")])

    def test_partial_orfs(self):
        self.assertEqual(find("ATGAAAGGG"), [])
        (orf,) = find("ATGAAAGGG", partial=True)
        self.assertTrue(orf.partial)
        self.assertEqual((orf.stop_codon, orf.protein, orf.end), ("", "MKG", 9))

    def test_partial_orf_ignores_an_incomplete_last_codon(self):
        (orf,) = find("ATGAAAGGGAC", partial=True)
        self.assertEqual((orf.end, orf.protein), (9, "MKG"))

    def test_alternative_start_codon_is_still_methionine(self):
        (orf,) = find("GTGAAATAG", starts=("GTG",))
        self.assertEqual((orf.start_codon, orf.protein), ("GTG", "MK"))
        self.assertEqual(find("GTGAAATAG"), [])

    def test_genetic_code_changes_the_stop_codons(self):
        seq = "ATGAGATAA"
        (standard,) = find(seq)
        self.assertEqual((standard.protein, standard.stop_codon), ("MR", "TAA"))
        (mito,) = find(seq, code=get_code(2))
        self.assertEqual((mito.protein, mito.stop_codon), ("M", "AGA"))

    def test_ambiguous_bases_never_start_or_stop_an_orf(self):
        (orf,) = find("ATGANATAG")
        self.assertEqual(orf.protein, "MX")
        self.assertEqual(find("ATGAAATNG"), [])

    def test_orfs_are_sorted_and_numbered(self):
        orfs = find("ATGAAATAG" + "CC" + "ATGCCCTAA")
        self.assertEqual([o.start for o in orfs], sorted(o.start for o in orfs))
        self.assertEqual([o.orf_id for o in orfs], [f"s_ORF{i}" for i in range(1, len(orfs) + 1)])

    def test_short_and_empty_sequences(self):
        for seq in ("A", "AT", "ATG", ""):
            self.assertEqual(find(seq), [])


class ParamValidationTests(unittest.TestCase):
    def test_bad_values(self):
        for overrides in ({"strands": "sideways"}, {"frames": ()}, {"frames": (4,)},
                          {"min_aa": 0}, {"starts": ()}, {"starts": ("AT",)},
                          {"starts": ("TAA",)}, {"starts": ("ANG",)}):
            with self.assertRaises(ValueError, msg=str(overrides)):
                params(**overrides)

    def test_starts_are_normalised(self):
        p = params(starts=("atg", "AUG", "gtg"))
        self.assertEqual(p.starts, ("ATG", "GTG"))

    def test_frames_are_sorted_and_deduplicated(self):
        self.assertEqual(params(frames=(3, 1, 3)).frames, (1, 3))


def oracle(seq, p):
    """A deliberately different way of finding ORFs: split codons at stops."""
    found = set()
    total = len(seq)
    for strand in p.strand_symbols:
        s = seq if strand == "+" else reverse_complement(seq)
        for frame in p.frames:
            codons = [s[i:i + 3] for i in range(frame - 1, total - 2, 3)]

            def emit(begin, end, stopped):
                starts = [j for j in range(begin, end) if codons[j] in p.starts]
                if not p.nested:
                    starts = starts[:1]
                for j in starts:
                    if end - j >= p.min_aa and (stopped or p.partial):
                        a = j * 3 + frame - 1
                        b = (end + (1 if stopped else 0)) * 3 + frame - 1
                        span = (a + 1, b) if strand == "+" else (total - b + 1, total - a)
                        found.add((strand, frame, span[0], span[1], stopped))

            begin = 0
            for k, codon in enumerate(codons):
                if codon in p.code.stops:
                    emit(begin, k, True)
                    begin = k + 1
            emit(begin, len(codons), False)
    return found


class AgainstOracleTests(unittest.TestCase):
    def test_random_sequences_agree_with_the_oracle(self):
        rng = random.Random(7)
        settings = [
            dict(),
            dict(min_aa=3),
            dict(nested=True, min_aa=2),
            dict(partial=True),
            dict(partial=True, nested=True, strands="reverse"),
            dict(frames=(2,), strands="forward", partial=True),
            dict(code=get_code(2), starts=get_code(2).starts, partial=True),
            dict(code=get_code(11), starts=get_code(11).starts, nested=True),
        ]
        for _ in range(120):
            alphabet = "ACGT" if rng.random() < 0.8 else "ACGTN"
            seq = "".join(rng.choice(alphabet) for _ in range(rng.randint(0, 260)))
            for overrides in settings:
                p = params(**overrides)
                got = {(o.strand, o.frame, o.start, o.end, not o.partial)
                       for o in find_orfs(Record("r", "", seq), p)}
                self.assertEqual(got, oracle(seq, p), (seq, overrides))

    def test_every_orf_is_internally_consistent(self):
        rng = random.Random(11)
        seq = "".join(rng.choice("ACGT") for _ in range(3000))
        p = params(partial=True, nested=True, min_aa=5)
        for orf in find_orfs(Record("r", "", seq), p):
            nt = orf.nucleotides(seq)
            self.assertEqual(len(nt) % 3, 0)
            self.assertEqual(nt[:3], orf.start_codon)
            self.assertIn(orf.start_codon, p.starts)
            body = nt if orf.partial else nt[:-3]
            if not orf.partial:
                self.assertEqual(nt[-3:], orf.stop_codon)
                self.assertIn(orf.stop_codon, p.code.stops)
            codons = [body[i:i + 3] for i in range(0, len(body), 3)]
            self.assertFalse(any(c in p.code.stops for c in codons))
            self.assertEqual(orf.protein, "M" + p.code.translate(body[3:]))
            self.assertEqual(orf.nt_length, len(nt))


if __name__ == "__main__":
    unittest.main()
