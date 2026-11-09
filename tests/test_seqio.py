import gzip
import io
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from orf.seqio import SequenceError, load, parse, reverse_complement


class ParseTests(unittest.TestCase):
    def test_multi_fasta(self):
        records = parse(">one first record\nATGC\nGGCC\n>two\nTTTT\n")
        self.assertEqual([r.id for r in records], ["one", "two"])
        self.assertEqual(records[0].seq, "ATGCGGCC")
        self.assertEqual(records[0].description, "first record")
        self.assertEqual(records[1].description, "")

    def test_plain_sequence_uses_default_id(self):
        (record,) = parse("atg aaa\n tag\n", "mine")
        self.assertEqual((record.id, record.seq), ("mine", "ATGAAATAG"))

    def test_rna_is_converted_and_case_is_ignored(self):
        (record,) = parse(">r\nAugcuu\n")
        self.assertEqual(record.seq, "ATGCTT")

    def test_iupac_codes_are_accepted(self):
        (record,) = parse(">x\nACGTRYSWKMBDHVN\n")
        self.assertEqual(record.length, 15)

    def test_bad_character_reports_record_and_position(self):
        with self.assertRaises(SequenceError) as caught:
            parse(">good\nACGT\n>bad\nACGT*ACGT\n")
        message = str(caught.exception)
        self.assertIn("'bad'", message)
        self.assertIn("'*'", message)
        self.assertIn("position 5", message)

    def test_digits_are_rejected(self):
        with self.assertRaises(SequenceError):
            parse(">x\n1 ATGC\n")

    def test_header_without_sequence(self):
        with self.assertRaises(SequenceError) as caught:
            parse(">empty\n>next\nACGT\n")
        self.assertIn("'empty'", str(caught.exception))

    def test_duplicate_ids_are_made_unique(self):
        records = parse(">a\nAC\n>a\nGT\n>a\nTT\n")
        self.assertEqual([r.id for r in records], ["a", "a_2", "a_3"])

    def test_blank_header_gets_a_number(self):
        (record,) = parse(">\nACGT\n")
        self.assertEqual(record.id, "seq1")

    def test_empty_text(self):
        self.assertEqual(parse(""), [])
        self.assertEqual(parse("\n  \n"), [])

    def test_windows_line_endings(self):
        (record,) = parse(">x\r\nATG\r\nAAA\r\n")
        self.assertEqual(record.seq, "ATGAAA")


class LoadTests(unittest.TestCase):
    def test_reads_plain_and_gzip_files(self):
        with tempfile.TemporaryDirectory() as tmp:
            plain = Path(tmp) / "sample.fa"
            plain.write_text(">s\nATGAAATAG\n")
            packed = Path(tmp) / "sample2.fa.gz"
            with gzip.open(packed, "wt") as handle:
                handle.write(">s\nATGAAATAG\n")
            self.assertEqual(load(str(plain))[0].seq, "ATGAAATAG")
            self.assertEqual(load(str(packed))[0].seq, "ATGAAATAG")

    def test_file_without_header_is_named_after_the_file(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "plasmid.txt"
            path.write_text("ATGAAATAG\n")
            self.assertEqual(load(str(path))[0].id, "plasmid")

    def test_utf8_bom_is_ignored(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "bom.fa"
            path.write_bytes(b"\xef\xbb\xbf>s\nATG\n")
            self.assertEqual(load(str(path))[0].id, "s")

    def test_stdin(self):
        with mock.patch("sys.stdin", io.StringIO(">piped\nATGAAATAG\n")):
            self.assertEqual(load("-")[0].id, "piped")

    def test_problems_become_sequence_errors(self):
        with tempfile.TemporaryDirectory() as tmp:
            empty = Path(tmp) / "empty.fa"
            empty.write_text("\n")
            binary = Path(tmp) / "blob.bin"
            binary.write_bytes(b"\xff\xfe\x00\x80\x81")
            for target, fragment in ((str(Path(tmp) / "nope.fa"), "not found"),
                                     (tmp, "directory"),
                                     (str(empty), "no sequence"),
                                     (str(binary), "not a text file")):
                with self.assertRaises(SequenceError) as caught:
                    load(target)
                self.assertIn(fragment, str(caught.exception))


class ReverseComplementTests(unittest.TestCase):
    def test_plain_bases(self):
        self.assertEqual(reverse_complement("ATGC"), "GCAT")
        self.assertEqual(reverse_complement("AAACCC"), "GGGTTT")

    def test_ambiguity_codes(self):
        self.assertEqual(reverse_complement("RYKMBDHVN"), "NBDHVKMRY")

    def test_involution(self):
        seq = "ACGTRYSWKMBDHVN" * 3
        self.assertEqual(reverse_complement(reverse_complement(seq)), seq)


if __name__ == "__main__":
    unittest.main()
