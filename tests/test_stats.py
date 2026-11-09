import random
import unittest

from orf import stats
from orf.codes import get_code
from orf.finder import SearchParams, find_orfs
from orf.seqio import Record

UBIQUITIN = "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"


class SequenceStatsTests(unittest.TestCase):
    def test_gc_percent(self):
        self.assertEqual(stats.gc_percent("GGCC"), 100.0)
        self.assertEqual(stats.gc_percent("ATAT"), 0.0)
        self.assertEqual(stats.gc_percent("ATGC"), 50.0)
        self.assertEqual(stats.gc_percent("GCNN"), 100.0)    # N does not dilute the result
        self.assertEqual(stats.gc_percent(""), 0.0)

    def test_gc_windows_cover_the_sequence(self):
        seq = "A" * 100 + "G" * 100
        points = stats.gc_windows(seq, 50)
        self.assertEqual([round(v) for _, v in points], [0, 0, 100, 100])
        self.assertEqual(points[0][0], 25.5)

    def test_gc_windows_edge_cases(self):
        self.assertEqual(stats.gc_windows("", 50), [])
        self.assertEqual(len(stats.gc_windows("ACGT", 1000)), 1)
        self.assertEqual(stats.gc_windows("NNNN", 2), [])
        # a short tail is dropped rather than plotted as noise
        self.assertEqual(len(stats.gc_windows("A" * 105, 50)), 2)

    def test_auto_window_scales_with_length(self):
        self.assertEqual(stats.auto_window(100), 50)
        self.assertEqual(stats.auto_window(12000), 100)
        self.assertEqual(stats.auto_window(1_200_000), 10000)


class ProteinTests(unittest.TestCase):
    def test_mass_of_human_ubiquitin(self):
        # Expasy lists 8564.84 Da for UBQ_HUMAN
        self.assertAlmostEqual(stats.protein_mass_kda(UBIQUITIN) * 1000, 8564.84, delta=0.5)

    def test_unknown_residues_still_get_a_mass(self):
        self.assertGreater(stats.protein_mass_kda("MXK"), stats.protein_mass_kda("MK"))

    def test_residue_counts(self):
        counts = stats.residue_counts(["MKK", "MX"])
        self.assertEqual((counts["M"], counts["K"]), (2, 2))
        self.assertEqual(sum(counts.values()), 4)      # X is not counted
        self.assertEqual(set(counts), set(stats.AMINO_ACIDS))


class ChanceExpectationTests(unittest.TestCase):
    def test_uniform_stop_probability(self):
        freqs = {b: 0.25 for b in "ACGT"}
        self.assertAlmostEqual(stats.stop_probability(freqs, get_code(1).stops), 3 / 64)
        self.assertAlmostEqual(stats.stop_probability(freqs, get_code(2).stops), 4 / 64)

    def test_at_rich_dna_has_more_stops(self):
        at_rich = stats.base_frequencies("AT" * 500 + "GC" * 50)
        gc_rich = stats.base_frequencies("GC" * 500 + "AT" * 50)
        stops = get_code(1).stops
        self.assertGreater(stats.stop_probability(at_rich, stops), stats.stop_probability(gc_rich, stops))

    def test_flipping_the_strand_swaps_complementary_bases(self):
        freqs = stats.base_frequencies("AAAAAC")
        flipped = stats.flip_strand(freqs)
        self.assertAlmostEqual(flipped["T"], freqs["A"])
        self.assertAlmostEqual(flipped["G"], freqs["C"])

    def test_frequencies_survive_degenerate_input(self):
        self.assertAlmostEqual(sum(stats.base_frequencies("").values()), 1.0)
        self.assertTrue(all(v > 0 for v in stats.base_frequencies("AAAA").values()))

    def test_longer_orfs_are_rarer(self):
        values = [stats.expected_by_chance(10_000, aa, 3 / 64, 3) for aa in (30, 100, 300)]
        self.assertGreater(values[0], values[1])
        self.assertGreater(values[1], values[2])

    def test_support_thresholds(self):
        self.assertEqual(stats.support_label(1e-6), "strong")
        self.assertEqual(stats.support_label(0.01), "moderate")
        self.assertEqual(stats.support_label(0.5), "weak")

    def test_prediction_matches_random_dna(self):
        """On random DNA, stop-to-stop stretches should appear as often as predicted."""
        rng = random.Random(2024)
        n = 150_000
        seq = "".join(rng.choice("ACGT") for _ in range(n))
        code = get_code(1)
        sense = tuple(c for c in code.table if c not in code.stops)   # any sense codon opens a stretch
        for aa in (25, 60):
            p = SearchParams(code, min_aa=aa, starts=sense)
            seen = len(find_orfs(Record("r", "", seq), p))
            freqs = stats.base_frequencies(seq)
            predicted = sum(
                stats.expected_by_chance(n, aa, stats.stop_probability(f, code.stops), 3)
                for f in (freqs, stats.flip_strand(freqs)))
            self.assertAlmostEqual(seen / predicted, 1.0, delta=0.15, msg=f"{aa} aa: {seen} vs {predicted:.0f}")


class HistogramTests(unittest.TestCase):
    def test_counts_add_up(self):
        forward = [30, 31, 45, 80, 120]
        reverse = [33, 90, 200]
        bins = stats.length_histogram(forward, reverse)
        self.assertEqual(sum(b["forward"] for b in bins), len(forward))
        self.assertEqual(sum(b["reverse"] for b in bins), len(reverse))
        self.assertLessEqual(bins[0]["lo"], 30)
        self.assertGreaterEqual(bins[-1]["hi"] if bins[-1]["hi"] is not None else 10**9, 200)

    def test_bins_do_not_overlap(self):
        bins = stats.length_histogram(list(range(30, 400, 7)), [100, 250])
        for a, b in zip(bins, bins[1:]):
            self.assertEqual(a["hi"] + 1, b["lo"])

    def test_empty_and_single_value(self):
        self.assertEqual(stats.length_histogram([], []), [])
        (only,) = stats.length_histogram([42], [])
        self.assertEqual((only["forward"], only["reverse"]), (1, 0))

    def test_long_tail_is_folded_into_one_bin(self):
        bins = stats.length_histogram([30] * 40 + [35] * 20 + [5000], [])
        self.assertIsNone(bins[-1]["hi"])
        self.assertTrue(bins[-1]["label"].startswith("≥"))
        self.assertEqual(bins[-1]["forward"], 1)
        self.assertEqual(sum(b["forward"] for b in bins), 61)
        self.assertLess(len(bins), 30)

    def test_median(self):
        self.assertEqual(stats.median([]), 0)
        self.assertEqual(stats.median([5]), 5)
        self.assertEqual(stats.median([1, 2, 3, 4]), 2.5)
        self.assertEqual(stats.median([9, 1, 5]), 5)


if __name__ == "__main__":
    unittest.main()
