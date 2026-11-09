import unittest

from orf.codes import CODES, get_code


class GeneticCodeTests(unittest.TestCase):
    def test_every_table_covers_all_64_codons(self):
        for code in CODES.values():
            self.assertEqual(len(code.table), 64, code)

    def test_standard_table_spot_checks(self):
        t = get_code(1).table
        self.assertEqual(t["ATG"], "M")
        self.assertEqual(t["TGG"], "W")
        self.assertEqual({c for c, aa in t.items() if aa == "*"}, {"TAA", "TAG", "TGA"})
        # fourfold-degenerate families
        for family, aa in (("GC", "A"), ("GG", "G"), ("CC", "P"), ("AC", "T"), ("GT", "V")):
            self.assertEqual({t[family + b] for b in "ACGT"}, {aa})
        # 20 amino acids plus stop
        self.assertEqual(len(set(t.values())), 21)

    def test_vertebrate_mito_differs_in_four_codons(self):
        std, mito = get_code(1).table, get_code(2).table
        different = {c for c in std if std[c] != mito[c]}
        self.assertEqual(different, {"AGA", "AGG", "ATA", "TGA"})
        self.assertEqual(mito["TGA"], "W")
        self.assertEqual(mito["ATA"], "M")
        self.assertEqual(get_code(2).stops, {"TAA", "TAG", "AGA", "AGG"})

    def test_bacterial_code_has_standard_amino_acids_but_more_starts(self):
        self.assertEqual(get_code(11).table, get_code(1).table)
        self.assertIn("GTG", get_code(11).starts)
        self.assertNotIn("GTG", get_code(1).starts)

    def test_translate(self):
        code = get_code(1)
        self.assertEqual(code.translate("ATGGCCTAA"), "MA*")
        self.assertEqual(code.translate("ATGGCCTAAGGG", to_stop=True), "MA")
        self.assertEqual(code.translate("ATGGCCTAA", offset=1), "WP")   # TGG CCT, then AA is left over
        self.assertEqual(code.translate("ATGAN"), "M")                  # trailing partial codon ignored
        self.assertEqual(code.translate("ANGATG"), "XM")                # ambiguous codon
        self.assertEqual(code.translate(""), "")

    def test_unknown_table(self):
        with self.assertRaises(ValueError) as caught:
            get_code(99)
        self.assertIn("available", str(caught.exception))
        with self.assertRaises(ValueError):
            get_code("abc")


if __name__ == "__main__":
    unittest.main()
