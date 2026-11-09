import io
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from orf import __version__

from .helpers import DEMO, run_cli


class HelpTests(unittest.TestCase):
    def test_no_arguments_shows_the_guide(self):
        code, out, _ = run_cli()
        self.assertEqual(code, 0)
        for text in ("usage: orf <command>", "find", "translate", "demo", "orf demo --open"):
            self.assertIn(text, out)

    def test_version(self):
        code, out, _ = run_cli("--version")
        self.assertEqual((code, out.strip()), (0, f"orf {__version__}"))

    def test_find_help_documents_every_option(self):
        code, out, _ = run_cli("find", "--help")
        self.assertEqual(code, 0)
        for flag in ("--seq", "--code", "--strand", "--frames", "--min-length", "--starts", "--nested",
                     "--partial", "--outdir", "--no-report", "--json", "--top", "--gc-window", "--title",
                     "--open", "--quiet", "--no-color", "examples:", "report.html"):
            self.assertIn(flag, out)

    def test_every_command_has_help(self):
        for command in ("find", "translate", "demo"):
            code, out, _ = run_cli(command, "-h")
            self.assertEqual(code, 0)
            self.assertIn(f"orf {command}", out)


class FindTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.out = Path(self.tmp.name) / "results"

    def test_find_on_a_file(self):
        code, out, err = run_cli("find", str(DEMO), "-o", str(self.out), "--no-color")
        self.assertEqual((code, err), (0, ""))
        for name in ("report.html", "orfs.tsv", "proteins.faa", "orfs.fna"):
            self.assertTrue((self.out / name).is_file(), name)
        self.assertFalse((self.out / "summary.json").exists())
        self.assertIn("demo_operon", out)
        self.assertIn("ORFs per frame", out)
        self.assertIn("Longest ORFs", out)
        self.assertNotIn("\033[", out)
        html = (self.out / "report.html").read_text(encoding="utf-8")
        self.assertIn("ORF report: demo.fasta", html)
        self.assertIn("orf find", html)       # the reproduce-this-run command

    def test_options_reach_the_search(self):
        run_cli("find", str(DEMO), "-o", str(self.out), "-m", "200", "-s", "forward", "-f", "1,3", "-q")
        rows = (self.out / "orfs.tsv").read_text().splitlines()[1:]
        self.assertTrue(rows)
        for row in rows:
            cells = row.split("\t")
            self.assertEqual(cells[2], "+")
            self.assertIn(cells[3], ("1", "3"))
            self.assertGreaterEqual(int(cells[7]), 200)

    def test_inline_sequence(self):
        code, out, _ = run_cli("find", "--seq", "ATGAAATAG", "-m", "1", "-o", str(self.out))
        self.assertEqual(code, 0)
        self.assertIn("1 sequence", out)
        self.assertIn("sequence_ORF1", (self.out / "orfs.tsv").read_text())

    def test_stdin(self):
        with mock.patch("sys.stdin", io.StringIO(">piped\nATGAAATAG\n")):
            code, out, _ = run_cli("find", "-", "-m", "1", "-o", str(self.out), "--no-report")
        self.assertEqual(code, 0)
        self.assertIn("piped_ORF1", (self.out / "orfs.tsv").read_text())

    def test_quiet_prints_nothing(self):
        code, out, err = run_cli("find", str(DEMO), "-o", str(self.out), "-q")
        self.assertEqual((code, out, err), (0, "", ""))
        self.assertTrue((self.out / "report.html").is_file())

    def test_no_report_and_json(self):
        run_cli("find", str(DEMO), "-o", str(self.out), "-q", "--no-report", "--json")
        self.assertFalse((self.out / "report.html").exists())
        self.assertTrue((self.out / "summary.json").is_file())

    def test_no_orfs_gives_advice_and_still_writes_files(self):
        code, out, _ = run_cli("find", "--seq", "ACGTACGTACGT", "-o", str(self.out), "--no-color")
        self.assertEqual(code, 0)
        self.assertIn("No ORFs found", out)
        self.assertIn("--partial", out)
        self.assertTrue((self.out / "report.html").is_file())

    def test_open_launches_the_browser(self):
        with mock.patch("webbrowser.open") as opener:
            run_cli("find", str(DEMO), "-o", str(self.out), "-q", "--open")
        opener.assert_called_once()
        self.assertTrue(opener.call_args[0][0].startswith("file:"))
        with mock.patch("webbrowser.open") as opener:
            run_cli("find", str(DEMO), "-o", str(self.out), "-q", "--open", "--no-report")
        opener.assert_not_called()

    def test_custom_title_and_top(self):
        run_cli("find", str(DEMO), "-o", str(self.out), "-q", "--title", "My genome", "--top", "3")
        html = (self.out / "report.html").read_text(encoding="utf-8")
        self.assertIn("<title>My genome</title>", html)
        self.assertIn("Showing the longest 3 of", html)

    def test_alternative_starts_and_codes(self):
        seq = "GTGAAATAG"
        run_cli("find", "--seq", seq, "-m", "1", "-o", str(self.out), "-q")
        self.assertEqual((self.out / "orfs.tsv").read_text().count("\n"), 1)
        run_cli("find", "--seq", seq, "-m", "1", "-c", "11", "--starts", "table", "-o", str(self.out), "-q")
        self.assertEqual((self.out / "orfs.tsv").read_text().count("\n"), 2)

    def test_partial_and_nested_flags(self):
        run_cli("find", "--seq", "ATGATGAAAGGG", "-m", "1", "--partial", "--nested", "-s", "forward",
                "-f", "1", "-o", str(self.out), "-q")
        rows = (self.out / "orfs.tsv").read_text().splitlines()[1:]
        self.assertEqual(len(rows), 2)
        self.assertTrue(all(row.split("\t")[10] == "yes" for row in rows))


class ErrorTests(unittest.TestCase):
    def test_missing_file(self):
        code, out, err = run_cli("find", "no_such_file.fa", "-o", "unused")
        self.assertEqual(code, 1)
        self.assertIn("file not found", err)
        self.assertEqual(out, "")

    def test_bad_sequence_is_explained(self):
        code, _, err = run_cli("find", "--seq", "ATGXX", "-o", "unused")
        self.assertEqual(code, 1)
        self.assertIn("unexpected character 'X'", err)

    def test_usage_errors_exit_with_two(self):
        cases = [
            ("find",),
            ("find", "--seq", "ATG", str(DEMO)),
            ("find", str(DEMO), "-f", "4"),
            ("find", str(DEMO), "-m", "0"),
            ("find", str(DEMO), "-c", "5"),
            ("find", str(DEMO), "--starts", "TAA"),
            ("find", str(DEMO), "--starts", "ATGC"),
            ("find", str(DEMO), "-s", "sideways"),
            ("translate",),
        ]
        for args in cases:
            code, _, err = run_cli(*args)
            self.assertEqual(code, 2, args)
            self.assertIn("error", err, args)

    def test_unwritable_output_folder(self):
        with tempfile.TemporaryDirectory() as tmp:
            blocker = Path(tmp) / "file"
            blocker.write_text("in the way")
            code, _, err = run_cli("find", "--seq", "ATGAAATAG", "-m", "1", "-o", str(blocker / "sub"))
        self.assertEqual(code, 1)
        self.assertIn("could not write", err)


class TranslateTests(unittest.TestCase):
    def test_all_frames_on_both_strands(self):
        code, out, _ = run_cli("translate", "--seq", "ATGGCCTAA")
        self.assertEqual(code, 0)
        self.assertEqual(out.count(">"), 6)
        self.assertIn(">sequence strand=+ frame=1 3aa\nMA*\n", out)

    def test_frame_strand_and_to_stop(self):
        _, out, _ = run_cli("translate", "--seq", "ATGGCCTAAGGG", "-s", "forward", "-f", "1", "--to-stop")
        self.assertEqual(out, ">sequence strand=+ frame=1 2aa\nMA\n")

    def test_reverse_strand(self):
        # reverse complement of ATGGCCTAA is TTAGGCCAT
        _, out, _ = run_cli("translate", "--seq", "ATGGCCTAA", "-s", "reverse", "-f", "1")
        self.assertEqual(out, ">sequence strand=- frame=1 3aa\nLGH\n")

    def test_rna_output(self):
        _, out, _ = run_cli("translate", "--seq", "ATGTTT", "--rna", "-s", "forward")
        self.assertEqual(out, ">sequence strand=+ rna\nAUGUUU\n")

    def test_wrapping_and_output_file(self):
        with tempfile.TemporaryDirectory() as tmp:
            target = Path(tmp) / "out.faa"
            code, out, _ = run_cli("translate", "--seq", "ATG" * 10, "-s", "forward", "-f", "1",
                                   "-w", "4", "-o", str(target))
            self.assertEqual((code, out), (0, ""))
            lines = target.read_text().splitlines()
        self.assertEqual(lines[1:], ["MMMM", "MMMM", "MM"])

    def test_too_short_frames_are_skipped(self):
        _, out, _ = run_cli("translate", "--seq", "AT", "-s", "forward")
        self.assertEqual(out, "")


class DemoTests(unittest.TestCase):
    def test_demo_runs_without_input(self):
        with tempfile.TemporaryDirectory() as tmp:
            target = Path(tmp) / "demo_out"
            code, out, _ = run_cli("demo", "-o", str(target), "--no-color")
            self.assertEqual(code, 0)
            self.assertTrue((target / "report.html").is_file())
        self.assertIn("Next steps", out)
        self.assertIn("orf find --help", out)

    def test_the_bundled_example_exists(self):
        self.assertTrue(DEMO.is_file())
        self.assertTrue(DEMO.read_text().startswith(">demo_operon"))


if __name__ == "__main__":
    unittest.main()
