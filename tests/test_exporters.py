import json
import tempfile
import unittest
from pathlib import Path

from orf import exporters
from orf.analysis import analyze
from orf.codes import get_code
from orf.finder import SearchParams
from orf.seqio import load, reverse_complement

from .helpers import DEMO


def demo_analysis(**overrides):
    return analyze(load(str(DEMO)), SearchParams(get_code(1), **overrides))


class ExporterTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.dir = Path(self.tmp.name)
        self.analysis = demo_analysis()
        self.written = exporters.write_all(self.dir, self.analysis, "orf find demo", with_json=True)

    def test_reports_what_it_wrote(self):
        self.assertEqual(set(self.written), {"orfs.tsv", "proteins.faa", "orfs.fna", "summary.json"})
        for name in self.written:
            self.assertTrue((self.dir / name).is_file(), name)

    def test_tsv_has_one_row_per_orf(self):
        lines = (self.dir / "orfs.tsv").read_text(encoding="utf-8").splitlines()
        self.assertEqual(lines[0].split("\t"), list(exporters.TSV_COLUMNS))
        self.assertEqual(len(lines) - 1, len(self.analysis.orfs))
        first = dict(zip(exporters.TSV_COLUMNS, lines[1].split("\t")))
        orf = self.analysis.orfs[0]
        self.assertEqual(first["orf_id"], orf.orf_id)
        self.assertEqual(int(first["start"]), orf.start)
        self.assertEqual(int(first["aa_length"]), orf.aa_length)
        self.assertTrue(all(len(line.split("\t")) == len(exporters.TSV_COLUMNS) for line in lines))

    def test_protein_fasta_round_trips(self):
        text = (self.dir / "proteins.faa").read_text(encoding="utf-8")
        entries = [e for e in text.split(">") if e]
        self.assertEqual(len(entries), len(self.analysis.orfs))
        for entry, orf in zip(entries, self.analysis.orfs):
            header, *body = entry.splitlines()
            self.assertTrue(header.startswith(orf.orf_id))
            self.assertEqual("".join(body), orf.protein)
            self.assertTrue(all(len(line) <= 60 for line in body))

    def test_nucleotide_fasta_matches_the_source_sequence(self):
        sequences = {r.record.id: r.record.seq for r in self.analysis.results}
        text = (self.dir / "orfs.fna").read_text(encoding="utf-8")
        entries = [e for e in text.split(">") if e]
        self.assertEqual(len(entries), len(self.analysis.orfs))
        for entry, orf in zip(entries, self.analysis.orfs):
            body = "".join(entry.splitlines()[1:])
            piece = sequences[orf.seq_id][orf.start - 1:orf.end]
            self.assertEqual(body, piece if orf.strand == "+" else reverse_complement(piece))
            self.assertEqual(body[:3], orf.start_codon)

    def test_json_summary(self):
        data = json.loads((self.dir / "summary.json").read_text(encoding="utf-8"))
        self.assertEqual(data["command"], "orf find demo")
        self.assertEqual(data["parameters"]["min_aa"], 30)
        self.assertEqual(len(data["orfs"]), len(self.analysis.orfs))
        self.assertEqual(data["orfs"][0]["protein"], self.analysis.orfs[0].protein)
        self.assertEqual([s["id"] for s in data["sequences"]], ["demo_operon", "demo_island"])

    def test_files_use_unix_line_endings(self):
        for name in ("orfs.tsv", "proteins.faa", "orfs.fna"):
            self.assertNotIn(b"\r", (self.dir / name).read_bytes(), name)

    def test_empty_results_still_write_valid_files(self):
        empty = demo_analysis(min_aa=100_000)
        target = self.dir / "empty"
        exporters.write_all(target, empty, "cmd")
        self.assertEqual((target / "orfs.tsv").read_text().count("\n"), 1)    # header only
        self.assertEqual((target / "proteins.faa").read_text(), "")


if __name__ == "__main__":
    unittest.main()
