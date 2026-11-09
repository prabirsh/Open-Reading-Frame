"""The HTML report: one self-contained file, no network, no dependencies."""

from datetime import datetime

from . import __version__, charts, stats
from .charts import esc, fmt_count, strand_label

MAX_MAPS = 6

# Colours are the validated light/dark steps of the reference categorical palette:
# slot 1 (blue) = forward strand, slot 2 (orange) = reverse strand, slot 7 (violet) = residues.
_LIGHT = {
    "--page": "#f9f9f7", "--surface": "#fcfcfb", "--ink": "#0b0b0b", "--ink-2": "#52514e",
    "--grid": "#e1e0d9", "--axis": "#c3c2b7", "--border": "rgba(11,11,11,.10)",
    "--wash": "rgba(11,11,11,.06)", "--lane": "#f3f2ee",
    "--fwd": "#2a78d6", "--rev": "#eb6834", "--aa": "#6250d6", "--on-fill": "#ffffff",
}
_DARK = {
    "--page": "#0d0d0d", "--surface": "#1a1a19", "--ink": "#f0efec", "--ink-2": "#c3c2b7",
    "--grid": "#2c2c2a", "--axis": "#383835", "--border": "rgba(255,255,255,.10)",
    "--wash": "rgba(255,255,255,.07)", "--lane": "#232322",
    "--fwd": "#3987e5", "--rev": "#d95926", "--aa": "#9085e9", "--on-fill": "#ffffff",
}


def _vars(tokens):
    return ";".join(f"{k}:{v}" for k, v in tokens.items())


_THEMES = (
    f':root{{color-scheme:light;{_vars(_LIGHT)}}}'
    f'@media (prefers-color-scheme:dark){{:root:not([data-theme="light"]){{color-scheme:dark;{_vars(_DARK)}}}}}'
    f':root[data-theme="dark"]{{color-scheme:dark;{_vars(_DARK)}}}'
    f'@media print{{:root,:root[data-theme="dark"]{{color-scheme:light;{_vars(_LIGHT)}}}}}'
)

_STYLE = """
*{box-sizing:border-box}
html{-webkit-text-size-adjust:100%}
body{margin:0;background:var(--page);color:var(--ink);font:15px/1.55 system-ui,-apple-system,"Segoe UI",sans-serif}
a{color:var(--fwd)}
main{max-width:1120px;margin:0 auto;padding:0 20px 64px}
h1{font-size:26px;line-height:1.2;margin:0}
h2{font-size:19px;margin:44px 0 14px}
h3{font-size:16px;margin:0}
code,.mono{font-family:ui-monospace,"Cascadia Mono",Consolas,monospace;font-size:.9em}
.top{display:flex;gap:16px;align-items:flex-start;justify-content:space-between;padding:32px 0 8px}
.sub{color:var(--ink-2);margin:6px 0 0}
nav{display:flex;gap:18px;flex-wrap:wrap;margin:14px 0 4px;font-size:14px}
button{font:inherit;color:inherit}
.theme{border:1px solid var(--border);background:var(--surface);border-radius:8px;padding:6px 12px;cursor:pointer}
.tiles{display:grid;grid-template-columns:repeat(auto-fit,minmax(160px,1fr));gap:12px;margin-top:22px}
.tile{background:var(--surface);border:1px solid var(--border);border-radius:10px;padding:14px 16px}
.tile .k{color:var(--ink-2);font-size:13px}
.tile .v{font-size:26px;font-weight:600;line-height:1.25;margin-top:2px}
.tile .n{color:var(--ink-2);font-size:12.5px;margin-top:2px}
.card{background:var(--surface);border:1px solid var(--border);border-radius:12px;padding:20px;margin-bottom:16px}
.card-head{margin-bottom:10px}
.headline{margin:4px 0 0;color:var(--ink-2)}
.finds{margin:0;padding-left:20px}
.finds li{margin:6px 0}
.legend{display:flex;flex-wrap:wrap;gap:6px 18px;list-style:none;margin:8px 0 6px;padding:0;font-size:13px;color:var(--ink-2)}
.legend li{display:flex;align-items:center;gap:7px}
.sw{display:inline-block;width:12px;height:12px;border-radius:3px}
.sw.fwd{background:var(--fwd)}.sw.rev{background:var(--rev)}.sw.aa{background:var(--aa)}
.sw.partial{opacity:.5}
.how{margin:12px 0 0;font-size:13.5px;color:var(--ink-2);max-width:78ch}
.scroll{overflow-x:auto}
.gc-wrap{margin-top:16px}
.grid2{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:16px}
.grid2 .card{margin:0}
.span1{grid-column:span 1}.span2{grid-column:span 2}.span3{grid-column:1/-1}
@media (max-width:900px){.grid2{grid-template-columns:minmax(0,1fr)}.grid2 .card{grid-column:1/-1}}
svg{display:block;width:100%;height:auto}
.map,.gc{min-width:720px}
.cols{min-width:300px}
svg text{font:11px system-ui,-apple-system,"Segoe UI",sans-serif;fill:var(--ink-2)}
svg .value{fill:var(--ink);font-weight:600}
svg .orf-label{fill:var(--on-fill);font-size:10px;font-weight:600;pointer-events:none}
.grid{stroke:var(--grid);stroke-width:1}
.axis{stroke:var(--axis);stroke-width:1}
.lane{fill:var(--lane)}
.orf{pointer-events:all;stroke:transparent;stroke-width:8px;stroke-linejoin:round}
.orf.fwd,.bar.fwd{fill:var(--fwd)}
.orf.rev,.bar.rev{fill:var(--rev)}
.bar.aa{fill:var(--aa)}
.orf.partial{fill-opacity:.5}
.orf:hover{filter:brightness(1.15)}
.gc-line{fill:none;stroke:var(--ink-2);stroke-width:2;stroke-linejoin:round;stroke-linecap:round}
.gc-area{fill:var(--ink-2);fill-opacity:.1}
.gc-dot{fill:var(--ink-2);stroke:var(--surface);stroke-width:2}
.hit{fill:transparent;pointer-events:all}
.hit:hover,.col:hover .hit,.col:focus .hit{fill:var(--wash)}
.col{outline:none}
.col:focus-visible .hit{stroke:var(--ink-2)}
.col:hover .bar{filter:brightness(1.12)}
details{margin-top:12px;font-size:13.5px}
summary{cursor:pointer;color:var(--ink-2)}
.tablewrap{overflow-x:auto;margin-top:8px;max-height:420px;overflow-y:auto}
table{border-collapse:collapse;width:100%;font-size:13px}
th,td{text-align:left;padding:7px 10px;border-bottom:1px solid var(--grid);white-space:nowrap}
th{font-weight:600;color:var(--ink-2);position:sticky;top:0;background:var(--surface)}
th button{background:none;border:0;padding:0;cursor:pointer;font-weight:600}
th[aria-sort=ascending] button::after{content:" \\2191"}
th[aria-sort=descending] button::after{content:" \\2193"}
td.num,th.num{text-align:right}
tbody tr:hover{background:var(--wash)}
.dots{letter-spacing:1px;color:var(--ink-2)}
.tools{display:flex;gap:12px;align-items:center;flex-wrap:wrap;margin:6px 0}
input[type=search]{font:inherit;color:inherit;background:var(--page);border:1px solid var(--border);border-radius:8px;padding:6px 10px;min-width:240px}
.note{color:var(--ink-2);font-size:13.5px}
dl{display:grid;grid-template-columns:max-content 1fr;gap:6px 18px;margin:0}
dt{font-weight:600}
dd{margin:0;color:var(--ink-2)}
pre{background:var(--page);border:1px solid var(--border);border-radius:8px;padding:12px;overflow-x:auto;font-size:13px;margin:0}
footer{margin-top:40px;color:var(--ink-2);font-size:13px}
#tip{position:fixed;z-index:10;pointer-events:none;background:var(--ink);color:var(--surface);border-radius:6px;
padding:8px 10px;font-size:12px;line-height:1.45;white-space:pre-line;max-width:340px;box-shadow:0 4px 14px rgba(0,0,0,.25)}
#tip[hidden]{display:none}
@media print,(forced-colors:active){.orf.rev,.bar.rev{fill:url(#hatch)}.sw.rev{background:repeating-linear-gradient(135deg,var(--rev) 0 2px,transparent 2px 5px)}}
@media print{.theme,nav,.tools,#tip,details{display:none}.card{break-inside:avoid}body{background:#fff}}
"""

_SCRIPT = """
(function () {
  var root = document.documentElement, tip = document.getElementById('tip');
  var btn = document.getElementById('theme');
  if (btn) btn.addEventListener('click', function () {
    var dark = root.dataset.theme ? root.dataset.theme === 'dark'
      : window.matchMedia('(prefers-color-scheme: dark)').matches;
    root.dataset.theme = dark ? 'light' : 'dark';
  });

  function place(x, y) {
    var w = tip.offsetWidth, h = tip.offsetHeight;
    var left = x + 14, top = y + 16;
    if (left + w > window.innerWidth - 8) left = x - w - 14;
    if (top + h > window.innerHeight - 8) top = y - h - 12;
    tip.style.left = Math.max(8, left) + 'px';
    tip.style.top = Math.max(8, top) + 'px';
  }
  function owner(node) { return node && node.closest ? node.closest('[data-tip]') : null; }
  document.addEventListener('pointermove', function (e) {
    var t = owner(e.target);
    if (!t) { tip.hidden = true; return; }
    tip.textContent = t.getAttribute('data-tip');
    tip.hidden = false;
    place(e.clientX, e.clientY);
  });
  document.addEventListener('pointerout', function (e) { if (!e.relatedTarget) tip.hidden = true; });
  document.addEventListener('focusin', function (e) {
    var t = owner(e.target);
    if (!t) return;
    var r = t.getBoundingClientRect();
    tip.textContent = t.getAttribute('data-tip');
    tip.hidden = false;
    place(r.left + r.width / 2, r.top);
  });
  document.addEventListener('focusout', function () { tip.hidden = true; });

  document.querySelectorAll('table.sortable').forEach(function (table) {
    var body = table.tBodies[0];
    Array.prototype.forEach.call(table.tHead.rows[0].cells, function (th, col) {
      var b = th.querySelector('button');
      if (!b) return;
      b.addEventListener('click', function () {
        var asc = th.getAttribute('aria-sort') !== 'ascending';
        var numeric = th.classList.contains('num') || th.getAttribute('data-type') === 'num';
        Array.prototype.forEach.call(th.parentNode.cells, function (c) { c.removeAttribute('aria-sort'); });
        th.setAttribute('aria-sort', asc ? 'ascending' : 'descending');
        var rows = Array.prototype.slice.call(body.rows);
        rows.sort(function (a, c) {
          var x = a.cells[col].getAttribute('data-v'), y = c.cells[col].getAttribute('data-v');
          if (numeric) { x = parseFloat(x); y = parseFloat(y); return asc ? x - y : y - x; }
          var order = x.localeCompare(y, undefined, { numeric: true });
          return asc ? order : -order;
        });
        rows.forEach(function (r) { body.appendChild(r); });
      });
    });
  });

  var box = document.getElementById('filter');
  if (box) box.addEventListener('input', function () {
    var q = box.value.trim().toLowerCase();
    Array.prototype.forEach.call(document.querySelectorAll('#orf-table tbody tr'), function (r) {
      r.hidden = Boolean(q) && r.textContent.toLowerCase().indexOf(q) === -1;
    });
  });
})();
"""


def _plural(n, word, many=None):
    return f"{n:,} {word if n == 1 else (many or word + 's')}"


def _dots(support):
    return {"strong": "●●●", "moderate": "●●○", "weak": "●○○"}[support]


def _tile(label, value, note=""):
    extra = f'<div class="n">{esc(note)}</div>' if note else ""
    return (f'<div class="tile"><div class="k">{esc(label)}</div>'
            f'<div class="v">{esc(value)}</div>{extra}</div>')


def _legend(entries):
    items = "".join(f'<li><i class="sw {cls}"></i>{esc(text)}</li>' for cls, text in entries)
    return f'<ul class="legend">{items}</ul>'


def _table(headers, rows, numeric=(), sortable=False, table_id="", sort_as_number=()):
    """rows are lists of (display, sort_value) pairs.

    `numeric` columns are right-aligned and sort as numbers; `sort_as_number`
    adds columns that sort as numbers but stay left-aligned.
    """
    head = []
    for i, h in enumerate(headers):
        cls = ' class="num"' if i in numeric else ""
        kind = ' data-type="num"' if i in sort_as_number else ""
        label = f"<button>{esc(h)}</button>" if sortable else esc(h)
        head.append(f"<th{cls}{kind}>{label}</th>")
    body = []
    for row in rows:
        cells = []
        for i, (shown, value) in enumerate(row):
            cls = ' class="num"' if i in numeric else ""
            cells.append(f'<td{cls} data-v="{esc(value)}">{shown}</td>')
        body.append(f"<tr>{''.join(cells)}</tr>")
    ident = f' id="{table_id}"' if table_id else ""
    kind = ' class="sortable"' if sortable else ""
    return (f'<table{ident}{kind}><thead><tr>{"".join(head)}</tr></thead>'
            f'<tbody>{"".join(body)}</tbody></table>')


def _figure(title, headline, graphic, how, legend="", table_view=None, span=1):
    twin = ""
    if table_view:
        twin = (f'<details><summary>Show this chart as a table</summary>'
                f'<div class="tablewrap">{table_view}</div></details>')
    return (f'<section class="card span{span}"><div class="card-head"><h3>{esc(title)}</h3>'
            f'<p class="headline">{esc(headline)}</p></div>{legend}'
            f'<div class="scroll">{graphic}</div>'
            f'<p class="how"><strong>How to read this.</strong> {how}</p>{twin}</section>')


# --- sections --------------------------------------------------------------------

def _findings(analysis):
    orfs = analysis.orfs
    params = analysis.params
    n_seq = len(analysis.results)
    if not orfs:
        return ["No ORFs of at least "
                f"{_plural(params.min_aa, 'amino acid')} were found. Lower the minimum length, "
                "allow partial ORFs, or try alternative start codons to widen the search."]
    forward = sum(1 for o in orfs if o.strand == "+")
    reverse = len(orfs) - forward
    items = [
        f"Found <strong>{_plural(len(orfs), 'ORF')}</strong> of at least "
        f"{_plural(params.min_aa, 'amino acid')} in {_plural(n_seq, 'sequence')}: "
        f"{forward:,} on the forward strand, {reverse:,} on the reverse."]
    best = analysis.longest
    items.append(
        f"The longest is <span class='mono'>{esc(best.orf_id)}</span>: {best.aa_length:,} amino acids "
        f"({best.mass_kda:.1f} kDa), frame {strand_label(best.strand, best.frame)}, "
        f"positions {fmt_count(best.start)}–{fmt_count(best.end)}.")
    strong = [o for o in orfs if o.support == "strong"]
    if strong:
        items.append(
            f"<strong>{len(strong):,}</strong> {'ORF is' if len(strong) == 1 else 'ORFs are'} strongly supported: "
            "chance alone would produce one that long fewer than once in a thousand searches like this one. "
            "These are your best candidates for real genes.")
    else:
        items.append(
            "No ORF is strongly supported. Every ORF here is about as long as chance would produce in "
            "sequence with this base composition, so treat them as candidates at most.")
    partial = sum(1 for o in orfs if o.partial)
    if partial:
        items.append(f"{_plural(partial, 'ORF')} run off the end of the sequence without a stop codon (shown paler).")
    if n_seq == 1:
        items.append(f"ORFs cover {analysis.results[0].coverage:.0f}% of the sequence.")
    ambiguous = sum(r.ambiguous for r in analysis.results)
    if ambiguous:
        items.append(f"{_plural(ambiguous, 'base')} in the input are ambiguous (N or other IUPAC codes). "
                     "Codons containing them are translated as X and can never start or end an ORF.")
    return items


def _sequence_table(analysis):
    rows = []
    for r in analysis.results:
        best = r.longest
        longest = f"{best.aa_length:,} aa ({strand_label(best.strand, best.frame)})" if best else "–"
        rows.append([
            (f'<span class="mono">{esc(r.record.id)}</span>', r.record.id),
            (fmt_count(r.record.length), r.record.length),
            (f"{r.gc:.1f}", f"{r.gc:.3f}"),
            (fmt_count(len(r.orfs)), len(r.orfs)),
            (esc(longest), best.aa_length if best else 0),
            (f"{r.coverage:.0f}%", f"{r.coverage:.3f}"),
        ])
    return _table(["Sequence", "Length (nt)", "GC %", "ORFs", "Longest", "Covered"],
                  rows, numeric={1, 2, 3, 5}, sortable=True)


def _lanes(params):
    return [(s, f) for s in params.strand_symbols for f in params.frames]


def _map_legend(params, has_partial):
    entries = []
    if "+" in params.strand_symbols:
        entries.append(("fwd", "Forward strand (+)"))
    if "-" in params.strand_symbols:
        entries.append(("rev", f"Reverse strand ({charts.MINUS})"))
    if has_partial:
        entries.append(("fwd partial", "No stop codon (ran off the end)"))
    return _legend(entries)


def _record_section(result, params, index):
    rec = result.record
    dot = " · "
    meta = dot.join([f"{fmt_count(rec.length)} nt", f"GC {result.gc:.1f}%",
                     _plural(len(result.orfs), "ORF"), f"{result.coverage:.0f}% covered by ORFs"])
    if rec.description:
        meta += dot + rec.description
    head = f'<div class="card-head"><h3>{esc(rec.id)}</h3><p class="headline">{esc(meta)}</p></div>'
    best = result.longest
    how_html = ""
    if best:
        graphic = charts.orf_map(rec.length, result.orfs, _lanes(params))
        note = ""
        if len(result.orfs) > 4000:
            note = f'<p class="note">Showing the 4,000 longest of {len(result.orfs):,} ORFs.</p>'
        map_html = (f'{head}{_map_legend(params, any(o.partial for o in result.orfs))}'
                    f'<div class="scroll">{graphic}</div>{note}')
        how_html = (
            '<p class="how"><strong>How to read this.</strong> Each row is one reading frame. '
            "Rows marked + read the sequence as you gave it, left to right; "
            f"rows marked {charts.MINUS} read the opposite strand, right to left. "
            "A bar is one ORF: it starts at a start codon and runs to the next stop codon in that frame, "
            "and the arrow tip is the stop. Bar length is protein length, so the long bars matter most: "
            "stop codons turn up about every 20 codons in random DNA, so a long open stretch usually means "
            "something is encoded there. Hover a bar for its details; the full list is in the ORF table below.</p>")
    else:
        map_html = (f'{head}<p class="note">No ORFs of at least {_plural(params.min_aa, "amino acid")} '
                    "were found in this sequence.</p>")

    gc_html = ""
    if result.gc_points:
        gc_svg = charts.gc_track(result.gc_points, rec.length, result.gc, result.gc_window)
        rows = [[(f"{pos:,.0f}", f"{pos:.0f}"), (f"{value:.1f}", f"{value:.3f}")]
                for pos, value in result.gc_points]
        twin = _table(["Window midpoint (bp)", "GC %"], rows, numeric={0, 1})
        gc_html = (f'<div class="scroll gc-wrap">{gc_svg}</div>'
                   f'<p class="how"><strong>GC track.</strong> Share of G and C bases in each '
                   f'{result.gc_window:,} bp window, lined up with the map above. Genes often differ in GC '
                   "content from the DNA around them, so dips and plateaus can mark inserted DNA or gene "
                   f"boundaries. The flat line is the sequence average ({result.gc:.1f}%).</p>"
                   f'<details><summary>Show the GC track as a table</summary>'
                   f'<div class="tablewrap">{twin}</div></details>')
    return f'<section class="card" id="map-{index}">{map_html}{how_html}{gc_html}</section>'


def _frame_chart(analysis):
    params = analysis.params
    lanes = _lanes(params)
    if len(lanes) < 2:
        return ""
    counts = {lane: 0 for lane in lanes}
    longest = {lane: 0 for lane in lanes}
    for o in analysis.orfs:
        counts[(o.strand, o.frame)] += 1
        longest[(o.strand, o.frame)] = max(longest[(o.strand, o.frame)], o.aa_length)
    items = []
    for lane in lanes:
        label = strand_label(*lane)
        items.append({
            "label": label,
            "tip": f"Frame {label}\n{_plural(counts[lane], 'ORF')}\nlongest {longest[lane]} aa",
            "parts": [("fwd" if lane[0] == "+" else "rev", counts[lane])],
        })
    top = max(counts.values())
    if top == 0:
        return ""
    busiest = max(lanes, key=lambda lane: counts[lane])
    headline = (f"Frame {strand_label(*busiest)} holds the most ORFs "
                f"({counts[busiest]:,} of {len(analysis.orfs):,}).")
    graphic = charts.columns(items, top, width=340, height=240, x_title="Reading frame",
                             y_title="ORFs", cap_labels=True)
    table = _table(["Frame", "ORFs", "Longest (aa)"],
                   [[(strand_label(*l), strand_label(*l)), (fmt_count(counts[l]), counts[l]),
                     (fmt_count(longest[l]), longest[l])] for l in lanes], numeric={1, 2})
    return _figure(
        "ORFs per reading frame", headline, graphic,
        "Counts ORFs by the frame they sit in. Real coding DNA tends to pile long ORFs into a few frames, "
        "while chance ORFs spread more evenly across all of them.",
        legend=_map_legend(params, False), table_view=table)


def _length_chart(analysis, span):
    orfs = analysis.orfs
    fwd = [o.aa_length for o in orfs if o.strand == "+"]
    rev = [o.aa_length for o in orfs if o.strand == "-"]
    bins = stats.length_histogram(fwd, rev)
    if not bins:
        return ""
    items = []
    for b in bins:
        bin_range = f"≥{b['lo']} aa" if b["hi"] is None else f"{b['lo']}–{b['hi']} aa"
        total = b["forward"] + b["reverse"]
        items.append({
            "label": b["label"],
            "tip": f"{bin_range}\n{_plural(total, 'ORF')} ({b['forward']} forward, {b['reverse']} reverse)",
            "parts": [("fwd", b["forward"]), ("rev", b["reverse"])],
        })
    top = max(b["forward"] + b["reverse"] for b in bins)
    lengths = fwd + rev
    headline = (f"The median ORF is {stats.median(lengths):.0f} aa long; the longest is {max(lengths):,} aa.")
    graphic = charts.columns(items, top, width=700 if span == 2 else 1000, height=240,
                             x_title="ORF length in amino acids (start of each range)", y_title="ORFs")
    table = _table(["Length range (aa)", "Forward", "Reverse", "Total"],
                   [[(esc(f"≥{b['lo']}" if b["hi"] is None else f"{b['lo']}–{b['hi']}"), b["lo"]),
                     (str(b["forward"]), b["forward"]), (str(b["reverse"]), b["reverse"]),
                     (str(b["forward"] + b["reverse"]), b["forward"] + b["reverse"])] for b in bins],
                   numeric={1, 2, 3})
    return _figure(
        "How long are the ORFs?", headline, graphic,
        "Each bar counts the ORFs in a length range, split by strand. Most ORFs in any DNA are short and "
        "arise by chance, so a pile on the left is normal. The bars far to the right are the ones worth a closer look.",
        legend=_map_legend(analysis.params, False), table_view=table, span=span)


def _residue_chart(analysis):
    counts = stats.residue_counts(o.protein for o in analysis.orfs)
    total = sum(counts.values())
    if not total:
        return ""
    ranked = sorted(counts, key=lambda aa: (-counts[aa], aa))
    items = []
    for aa in ranked:
        pct = 100.0 * counts[aa] / total
        name = stats.AMINO_ACID_NAMES[aa]
        items.append({
            "label": aa,
            "tip": f"{name.capitalize()} ({aa})\n{pct:.1f}% of residues\n{counts[aa]:,} residues",
            "parts": [("aa", pct)],
        })
    top = 100.0 * counts[ranked[0]] / total
    low = ranked[-1]
    headline = (f"{stats.AMINO_ACID_NAMES[ranked[0]].capitalize()} is the most common residue "
                f"({top:.1f}%); {stats.AMINO_ACID_NAMES[low]} is the rarest "
                f"({100.0 * counts[low] / total:.1f}%).")
    graphic = charts.columns(items, top, width=1000, height=240, x_title="Amino acid (one-letter code)",
                             y_title="% of residues")
    table = _table(["Amino acid", "Residues", "Share (%)"],
                   [[(f"{stats.AMINO_ACID_NAMES[aa].capitalize()} ({aa})", aa),
                     (f"{counts[aa]:,}", counts[aa]), (f"{100.0 * counts[aa] / total:.1f}", f"{100.0 * counts[aa] / total:.3f}")]
                    for aa in ranked], numeric={1, 2})
    return _figure(
        "Amino-acid composition", headline, graphic,
        "The share of each amino acid across every translated ORF, most common first. Real proteins usually "
        "have leucine on top and tryptophan and cysteine near the bottom, so a very different profile suggests "
        "the ORFs are chance matches or the DNA has an unusual base composition.",
        table_view=table, span=3)


def _orf_table(analysis, top):
    ranked = sorted(analysis.orfs, key=lambda o: (-o.aa_length, o.seq_id, o.start))
    shown = ranked[:top]
    rows = []
    for o in shown:
        stop = o.stop_codon or "none"
        rows.append([
            (f'<span class="mono">{esc(o.orf_id)}</span>', o.orf_id),
            (strand_label(o.strand, o.frame), f"{o.strand}{o.frame}"),
            (fmt_count(o.start), o.start),
            (fmt_count(o.end), o.end),
            (fmt_count(o.aa_length), o.aa_length),
            (f"{o.mass_kda:.1f}", f"{o.mass_kda:.3f}"),
            (f"{o.gc:.0f}", f"{o.gc:.3f}"),
            (f"{o.start_codon} → {stop}", o.start_codon + stop),
            (f'<span class="dots" aria-hidden="true">{_dots(o.support)}</span> {o.support} '
             f'<span class="note">(E={o.expect:.2g})</span>', f"{o.expect:.6e}"),
        ])
    table = _table(["ORF", "Frame", "Start", "End", "Length (aa)", "kDa", "GC %", "Start → stop", "Support"],
                   rows, numeric={2, 3, 4, 5, 6}, sortable=True, table_id="orf-table",
                   sort_as_number={8})
    cut = ""
    if len(ranked) > len(shown):
        cut = (f' Showing the longest {len(shown):,} of {len(ranked):,}; the complete list is in '
               "<code>orfs.tsv</code>.")
    return (f'<section class="card" id="table"><div class="card-head"><h3>ORF table</h3>'
            f'<p class="headline">Longest first. Click a heading to sort.{cut}</p></div>'
            '<div class="tools"><input id="filter" type="search" placeholder="Filter by name, frame, codon..." '
            'aria-label="Filter the ORF table"></div>'
            f'<div class="tablewrap">{table}</div></section>')


def _about(analysis, command):
    p = analysis.params
    strands = {"both": "both strands", "forward": "forward strand only", "reverse": "reverse strand only"}[p.strands]
    settings = [
        ("Genetic code", f"{p.code.id} – {p.code.name}"),
        ("Strands", strands),
        ("Frames", ", ".join(str(f) for f in p.frames)),
        ("Minimum length", f"{p.min_aa} amino acids"),
        ("Start codons", ", ".join(p.starts)),
        ("Nested starts", "every start codon is reported" if p.nested else "only the outermost start per stop"),
        ("Partial ORFs", "kept" if p.partial else "dropped (a stop codon is required)"),
    ]
    dl = "".join(f"<dt>{esc(k)}</dt><dd>{esc(v)}</dd>" for k, v in settings)
    terms = [
        ("ORF", "open reading frame: a start codon followed by codons with no in-frame stop, ending at a stop codon."),
        ("Reading frame", "one of three ways to cut a strand into codons, depending on whether you start at base 1, 2 or 3. "
                          "Frames +1 to +3 are on the strand you gave; frames −1 to −3 are on its reverse complement."),
        ("Coordinates", "1-based and inclusive, always counted along the sequence you gave, even for reverse-strand ORFs."),
        ("Support", "how surprising the ORF's length is. E is the number of ORFs at least this long that chance would produce "
                    "in a search like this one, from random sequence with the same base composition. "
                    "Strong means E &lt; 0.001, moderate E &lt; 0.05, weak otherwise. "
                    "It is a screening aid, not a gene predictor: real short genes can look weak."),
        ("kDa", "average protein mass computed from the translated residues."),
    ]
    gl = "".join(f"<dt>{esc(k)}</dt><dd>{v if k == 'Support' else esc(v)}</dd>" for k, v in terms)
    return (f'<section class="card" id="about"><h3>Settings used</h3><div style="height:10px"></div><dl>{dl}</dl>'
            '<h3 style="margin-top:22px">Glossary</h3><div style="height:10px"></div>'
            f'<dl>{gl}</dl>'
            '<h3 style="margin-top:22px">Reproduce this run</h3><div style="height:10px"></div>'
            f'<pre>{esc(command)}</pre></section>')


def render(analysis, title="ORF report", top=50, command="", input_label=""):
    results = analysis.results
    orfs = analysis.orfs
    params = analysis.params
    strong = sum(1 for o in orfs if o.support == "strong")
    longest = analysis.longest
    tiles = [
        _tile("Sequences", f"{len(results):,}", input_label),
        _tile("Total length", f"{analysis.total_length:,} nt"),
        _tile("GC content", f"{analysis.gc:.1f}%"),
        _tile("ORFs found", f"{len(orfs):,}", f"≥ {params.min_aa} aa"),
        _tile("Longest ORF", f"{longest.aa_length:,} aa" if longest else "–",
              f"frame {strand_label(longest.strand, longest.frame)}" if longest else ""),
        _tile("Strongly supported", f"{strong:,}", f"of {len(orfs):,} ORFs" if orfs else ""),
    ]
    findings = "".join(f"<li>{item}</li>" for item in _findings(analysis))

    maps = [_record_section(r, params, i) for i, r in enumerate(results[:MAX_MAPS])]
    more = ""
    if len(results) > MAX_MAPS:
        more = (f'<p class="note">Maps are drawn for the first {MAX_MAPS} of {len(results):,} sequences. '
                "Every ORF is still in the table below and in <code>orfs.tsv</code>.</p>")
    seq_table = ""
    if len(results) > 1:
        seq_table = (f'<section class="card"><div class="card-head"><h3>Sequences</h3></div>'
                     f'<div class="tablewrap">{_sequence_table(analysis)}</div></section>')

    frames_html = _frame_chart(analysis)
    charts_html = "".join([frames_html, _length_chart(analysis, 2 if frames_html else 3),
                           _residue_chart(analysis)])
    charts_block = ""
    if charts_html:
        charts_block = f'<h2 id="charts">Across all sequences</h2><div class="grid2">{charts_html}</div>'
    table_block = _orf_table(analysis, top) if orfs else ""
    stamp = datetime.now().strftime("%Y-%m-%d %H:%M")

    hatch = ('<svg width="0" height="0" aria-hidden="true" style="position:absolute"><defs>'
             '<pattern id="hatch" width="6" height="6" patternUnits="userSpaceOnUse" patternTransform="rotate(135)">'
             '<rect width="6" height="6" style="fill:var(--rev)"/>'
             '<line x1="0" y1="0" x2="0" y2="6" style="stroke:var(--surface);stroke-width:2"/></pattern></defs></svg>')

    return f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{esc(title)}</title>
<style>{_THEMES}{_STYLE}</style>
</head>
<body>
{hatch}
<main>
<div class="top">
<div>
<h1>{esc(title)}</h1>
<p class="sub">Open reading frames · {esc(stamp)} · orf {__version__}</p>
<nav><a href="#findings">Findings</a><a href="#maps">Maps</a><a href="#charts">Charts</a>
<a href="#table">ORF table</a><a href="#about">About</a></nav>
</div>
<button class="theme" id="theme" type="button">Light / dark</button>
</div>
<div class="tiles">{"".join(tiles)}</div>
<h2 id="findings">What was found</h2>
<section class="card"><ul class="finds">{findings}</ul></section>
{seq_table}
<h2 id="maps">Maps</h2>
{"".join(maps)}{more}
{charts_block}
<h2>ORFs</h2>
{table_block or '<p class="note">No ORFs to list.</p>'}
<h2>About this report</h2>
{_about(analysis, command)}
<footer>Made with orf {__version__}. The report is a single file; it works offline and prints cleanly.</footer>
</main>
<div id="tip" role="tooltip" hidden></div>
<script>{_SCRIPT}</script>
</body>
</html>
"""


def write(path, analysis, **kwargs):
    with open(path, "w", encoding="utf-8", newline="\n") as out:
        out.write(render(analysis, **kwargs))
