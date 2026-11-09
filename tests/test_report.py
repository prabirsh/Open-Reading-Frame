import re
import unittest
import xml.etree.ElementTree as ET

from orf import charts, report
from orf.analysis import analyze
from orf.codes import get_code
from orf.finder import SearchParams
from orf.seqio import load, parse

from .helpers import DEMO


def build(records=None, **overrides):
    records = records or load(str(DEMO))
    analysis = analyze(records, SearchParams(get_code(1), **overrides))
    return analysis, report.render(analysis, title="Test report", command="orf find x")


class ReportTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.analysis, cls.html = build()

    def test_is_one_self_contained_page(self):
        self.assertTrue(self.html.startswith("<!doctype html>"))
        self.assertIn("<title>Test report</title>", self.html)
        for external in ("http://", "https://", "src=", "<link"):
            self.assertNotIn(external, self.html)

    def test_has_every_section(self):
        for marker in ("What was found", "Maps", "ORFs per reading frame", "How long are the ORFs?",
                       "Amino-acid composition", 'id="orf-table"', "Glossary", "Reproduce this run",
                       "orf find x", 'id="tip"'):
            self.assertIn(marker, self.html)

    def test_every_chart_is_well_formed_svg(self):
        svgs = re.findall(r"<svg.*?</svg>", self.html, re.S)
        self.assertGreaterEqual(len(svgs), 6)     # two maps, two GC tracks, three charts, hatch defs
        for svg in svgs:
            ET.fromstring(svg)

    def test_both_themes_are_defined(self):
        self.assertIn("prefers-color-scheme:dark", self.html)
        self.assertIn('data-theme="dark"', self.html)
        self.assertIn("#2a78d6", self.html)      # light forward
        self.assertIn("#3987e5", self.html)      # dark forward

    def test_findings_mention_the_longest_orf(self):
        longest = self.analysis.longest
        self.assertIn(longest.orf_id, self.html)
        self.assertIn(f"{longest.aa_length:,} amino acids", self.html)

    def test_table_lists_the_requested_number_of_orfs(self):
        rows = re.search(r'<table id="orf-table".*?</table>', self.html, re.S).group(0)
        self.assertEqual(rows.count("<tr>") - 1, min(50, len(self.analysis.orfs)))
        small = report.render(self.analysis, top=5)
        table = re.search(r'<table id="orf-table".*?</table>', small, re.S).group(0)
        self.assertEqual(table.count("<tr>") - 1, 5)
        self.assertIn("orfs.tsv", small)

    def test_nothing_unrendered_leaks_through(self):
        visible = re.sub(r"<(style|script)>.*?</\1>", "", self.html, flags=re.S)
        for leak in (r"\bNone\b", r"\bnan\b", r"\binf\b", r"[{}]", r"%s"):
            self.assertIsNone(re.search(leak, visible), leak)

    def test_user_text_is_escaped(self):
        records = parse(">evil<img_src=x> a <script>alert(1)</script>\nATGAAAGGGTAGATGAAATAGCCCATGAAATAG\n")
        _, html = build(records, min_aa=1)
        self.assertNotIn("<img_src", html)
        self.assertNotIn("<script>alert", html)
        self.assertIn("&lt;img_src=x&gt;", html)


class VariantTests(unittest.TestCase):
    def test_no_orfs(self):
        analysis, html = build(min_aa=10_000)
        self.assertIn("No ORFs", html)
        self.assertNotIn('id="orf-table"', html)
        self.assertNotIn('<svg class="cols"', html)
        self.assertIn("<svg class=\"gc\"", html)       # the GC track is still useful

    def test_single_strand_and_frame_skips_the_frame_chart(self):
        _, html = build(strands="forward", frames=(1,))
        self.assertNotIn("ORFs per reading frame", html)
        self.assertIn("How long are the ORFs?", html)
        self.assertNotIn("Reverse strand", html)

    def test_forward_only_has_no_reverse_lanes(self):
        _, html = build(strands="forward")
        self.assertNotIn(">−1<", html)
        self.assertIn(">+1<", html)

    def test_partial_orfs_are_explained(self):
        records = parse(">p\nATGAAAGGGCCCAAATTT\n")
        analysis, html = build(records, min_aa=1, partial=True)
        self.assertTrue(analysis.orfs[0].partial)
        self.assertIn("ran off the end", html)
        self.assertIn("without a stop codon", html)

    def test_many_sequences_get_a_cap_on_maps(self):
        text = "".join(f">s{i}\nATGAAATAGCCCATGAAATAG\n" for i in range(report.MAX_MAPS + 3))
        _, html = build(parse(text), min_aa=1)
        self.assertEqual(html.count('class="map"'), report.MAX_MAPS)
        self.assertIn(f"first {report.MAX_MAPS} of {report.MAX_MAPS + 3}", html)

    def test_ambiguous_bases_are_flagged(self):
        _, html = build(parse(">n\nATGAAANNNTAG\n"), min_aa=1)
        self.assertIn("ambiguous", html)


class ChartHelperTests(unittest.TestCase):
    def test_nice_ticks_start_at_zero_and_cover_the_top(self):
        for top in (1, 3, 7, 12, 99, 1234, 56_789):
            ticks = charts.nice_ticks(top)
            self.assertEqual(ticks[0], 0)
            self.assertGreaterEqual(ticks[-1], top)
            self.assertLessEqual(len(ticks), 8)

    def test_axis_labels(self):
        self.assertEqual(charts.fmt_axis(950), "950")
        self.assertEqual(charts.fmt_axis(12_000), "12K")
        self.assertEqual(charts.fmt_axis(2_000_000), "2M")

    def test_columns_never_divide_by_zero(self):
        svg = charts.columns([{"label": "a", "tip": "t", "parts": [("fwd", 0)]}], 0)
        ET.fromstring(svg)

    def test_tooltips_are_escaped(self):
        svg = charts.columns([{"label": "<b>", "tip": 'x"y<z', "parts": [("fwd", 3)]}], 3)
        ET.fromstring(svg)
        self.assertNotIn('x"y', svg)


if __name__ == "__main__":
    unittest.main()
