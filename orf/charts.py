"""Inline SVG charts for the report.

Nothing here sets a colour directly. Marks carry CSS classes (fwd, rev, aa, ...)
and the report stylesheet maps those to the light and dark palettes.
"""

import html
import math

MINUS = "−"


def esc(value):
    return html.escape(str(value), quote=True)


def num(x):
    text = f"{x:.1f}"
    return text[:-2] if text.endswith(".0") else text


def fmt_count(value):
    return f"{value:,.0f}"


def fmt_axis(value):
    if value >= 1e6:
        return f"{value / 1e6:g}M"
    if value >= 1e4:
        return f"{value / 1e3:g}K"
    return f"{value:,.0f}"


def strand_label(strand, frame):
    return f"{'+' if strand == '+' else MINUS}{frame}"


def nice_ticks(top, count=4):
    """Round tick values from 0 up to at least `top`."""
    top = max(top, 1)
    raw = top / count
    magnitude = 10 ** math.floor(math.log10(raw))
    step = max(1, next(m * magnitude for m in (1, 2, 5, 10) if raw <= m * magnitude * 1.0000001))
    ticks = [0]
    while ticks[-1] < top - 1e-9:
        ticks.append(round(ticks[-1] + step, 10))
    return ticks


def _top_rounded(x, y, w, h, r):
    """A bar with rounded top corners and a square foot."""
    r = max(0.0, min(r, w / 2, h))
    if r == 0:
        return f"M{num(x)},{num(y)}h{num(w)}v{num(h)}h{num(-w)}z"
    return (f"M{num(x)},{num(y + h)}V{num(y + r)}a{num(r)},{num(r)} 0 0 1 {num(r)},{num(-r)}"
            f"h{num(w - 2 * r)}a{num(r)},{num(r)} 0 0 1 {num(r)},{num(r)}V{num(y + h)}z")


# --- six-frame ORF map -------------------------------------------------------

_LEFT, _RIGHT = 56, 16
_LANE_H, _LANE_GAP = 16, 3


def _orf_tip(orf):
    stop = orf.stop_codon or "no stop"
    lines = [
        f"{orf.orf_id}  ({strand_label(orf.strand, orf.frame)})",
        f"{fmt_count(orf.start)}–{fmt_count(orf.end)} · {orf.aa_length} aa · {orf.mass_kda:.1f} kDa",
        f"{orf.start_codon} → {stop} · GC {orf.gc:.0f}%",
        f"support: {orf.support} (chance E = {orf.expect:.2g})",
    ]
    return "\n".join(lines)


def _tick_positions(length):
    ticks = [t for t in nice_ticks(length, 6) if 0 < t <= length]
    return ticks or [length]


def orf_map(length, orfs, lanes, width=980, shown_limit=4000):
    """Six-frame map. `lanes` is a list of (strand, frame), forward lanes first."""
    plot_w = width - _LEFT - _RIGHT
    drawn = sorted(orfs, key=lambda o: -o.aa_length)[:shown_limit]
    top3 = {id(o) for o in drawn[:3]}

    def x_of(pos):
        return _LEFT + pos / length * plot_w

    forward = [l for l in lanes if l[0] == "+"]
    reverse = [l for l in lanes if l[0] == "-"]
    y = 6
    lane_y = {}
    for lane in forward:
        lane_y[lane] = y
        y += _LANE_H + _LANE_GAP
    axis_y = y + 2
    y = axis_y + 30
    for lane in reverse:
        lane_y[lane] = y
        y += _LANE_H + _LANE_GAP
    height = y + 4 if reverse else axis_y + 32

    out = [f'<svg class="map" viewBox="0 0 {width} {height}" role="img" '
           f'aria-label="Six-frame map of open reading frames">']
    for lane, ly in lane_y.items():
        out.append(f'<rect class="lane" x="{_LEFT}" y="{ly}" width="{plot_w}" height="{_LANE_H}" rx="2"/>')
        out.append(f'<text class="tick" x="{_LEFT - 10}" y="{ly + _LANE_H - 4}" text-anchor="end">'
                   f'{strand_label(*lane)}</text>')
    out.append(f'<line class="axis" x1="{_LEFT}" x2="{_LEFT + plot_w}" y1="{axis_y + 6}" y2="{axis_y + 6}"/>')
    for t in _tick_positions(length):
        tx = num(x_of(t))
        out.append(f'<line class="axis" x1="{tx}" x2="{tx}" y1="{axis_y + 6}" y2="{axis_y + 11}"/>')
        out.append(f'<text class="tick" x="{tx}" y="{axis_y + 23}" text-anchor="middle">{fmt_axis(t)}</text>')

    labels = []
    for orf in drawn:
        ly = lane_y.get((orf.strand, orf.frame))
        if ly is None:
            continue
        x0, x1 = x_of(orf.start - 1), x_of(orf.end)
        x1 = max(x1, x0 + 1.5)
        tip = 0 if orf.partial else min(7.0, (x1 - x0) * 0.45)
        if tip == 0:
            pts = [(x0, ly), (x1, ly), (x1, ly + _LANE_H), (x0, ly + _LANE_H)]
        elif orf.strand == "+":
            pts = [(x0, ly), (x1 - tip, ly), (x1, ly + _LANE_H / 2),
                   (x1 - tip, ly + _LANE_H), (x0, ly + _LANE_H)]
        else:
            pts = [(x1, ly), (x0 + tip, ly), (x0, ly + _LANE_H / 2),
                   (x0 + tip, ly + _LANE_H), (x1, ly + _LANE_H)]
        cls = "orf fwd" if orf.strand == "+" else "orf rev"
        if orf.partial:
            cls += " partial"
        points = " ".join(f"{num(px)},{num(py)}" for px, py in pts)
        out.append(f'<polygon class="{cls}" points="{points}" data-tip="{esc(_orf_tip(orf))}"/>')
        text = f"{orf.aa_length} aa"
        if id(orf) in top3 and (x1 - x0 - tip) >= len(text) * 6 + 14:
            body_mid = (x0 + x1 - tip) / 2 if orf.strand == "+" else (x0 + tip + x1) / 2
            labels.append(f'<text class="orf-label" x="{num(body_mid)}" y="{num(ly + _LANE_H - 4.5)}" '
                          f'text-anchor="middle">{text}</text>')
    out.extend(labels)
    out.append("</svg>")
    return "".join(out)


# --- GC content track ----------------------------------------------------------

def gc_track(points, length, mean, window, width=980, height=104):
    if not points:
        return ""
    plot_w = width - _LEFT - _RIGHT
    top_pad, bottom_pad = 10, 8
    plot_h = height - top_pad - bottom_pad
    values = [v for _, v in points]
    lo = max(0, math.floor((min(values) - 5) / 10) * 10)
    hi = min(100, math.ceil((max(values) + 5) / 10) * 10)
    if hi - lo < 20:
        hi = min(100, lo + 20)
        lo = max(0, hi - 20)

    def x_of(pos):
        return _LEFT + pos / length * plot_w

    def y_of(v):
        return top_pad + (hi - v) / (hi - lo) * plot_h

    out = [f'<svg class="gc" viewBox="0 0 {width} {height}" role="img" '
           f'aria-label="GC content along the sequence">']
    for tick in (lo, (lo + hi) / 2, hi):
        ty = num(y_of(tick))
        out.append(f'<line class="grid" x1="{_LEFT}" x2="{_LEFT + plot_w}" y1="{ty}" y2="{ty}"/>')
        out.append(f'<text class="tick" x="{_LEFT - 10}" y="{num(y_of(tick) + 4)}" text-anchor="end">{tick:g}%</text>')

    mean_y = num(y_of(min(max(mean, lo), hi)))
    out.append(f'<line class="axis" x1="{_LEFT}" x2="{_LEFT + plot_w}" y1="{mean_y}" y2="{mean_y}"/>')

    if len(points) > 1:
        line = " ".join(f"{num(x_of(p))},{num(y_of(v))}" for p, v in points)
        base = top_pad + plot_h
        first_x, last_x = num(x_of(points[0][0])), num(x_of(points[-1][0]))
        out.append(f'<polygon class="gc-area" points="{first_x},{num(base)} {line} {last_x},{num(base)}"/>')
        out.append(f'<polyline class="gc-line" points="{line}"/>')
    else:
        pos, value = points[0]
        out.append(f'<circle class="gc-dot" cx="{num(x_of(pos))}" cy="{num(y_of(value))}" r="4"/>')

    column = max(plot_w / len(points), 4)
    for pos, value in points:
        a = max(1, int(pos - window / 2 + 0.5))
        b = min(length, int(pos + window / 2))
        tip = f"GC {value:.1f}%\n{fmt_count(a)}–{fmt_count(b)} bp"
        out.append(f'<rect class="hit" x="{num(x_of(pos) - column / 2)}" y="{top_pad}" '
                   f'width="{num(column)}" height="{plot_h}" data-tip="{esc(tip)}"/>')
    out.append("</svg>")
    return "".join(out)


# --- column charts (bars and stacked histogram) -----------------------------------

def columns(items, top, width=520, height=230, x_title="", y_title="", cap_labels=False):
    """Vertical columns. Each item: {label, tip, parts: [(css class, value), ...]}.

    Parts are stacked bottom to top with a 2px gap between them.
    """
    left, right, top_pad, bottom = 44, 10, 26, 44
    plot_w, plot_h = width - left - right, height - top_pad - bottom
    ticks = nice_ticks(top)
    scale_top = ticks[-1]
    base = top_pad + plot_h
    band = plot_w / max(len(items), 1)
    bar_w = min(24.0, band * 0.7)

    out = [f'<svg class="cols" viewBox="0 0 {width} {height}" role="img">']
    if y_title:
        out.append(f'<text class="tick" x="{left - 34}" y="12">{esc(y_title)}</text>')
    for tick in ticks:
        ty = base - tick / scale_top * plot_h
        out.append(f'<line class="grid" x1="{left}" x2="{left + plot_w}" y1="{num(ty)}" y2="{num(ty)}"/>')
        out.append(f'<text class="tick" x="{left - 8}" y="{num(ty + 4)}" text-anchor="end">{fmt_axis(tick)}</text>')
    out.append(f'<line class="axis" x1="{left}" x2="{left + plot_w}" y1="{base}" y2="{base}"/>')

    skip = 2 if len(items) > 14 and band < 36 else 1
    for i, item in enumerate(items):
        x = left + band * i + (band - bar_w) / 2
        total = sum(v for _, v in item["parts"])
        out.append(f'<g class="col" tabindex="0" data-tip="{esc(item["tip"])}" '
                   f'aria-label="{esc(item["tip"].replace(chr(10), ", "))}">')
        out.append(f'<rect class="hit" x="{num(left + band * i)}" y="{top_pad}" '
                   f'width="{num(band)}" height="{num(plot_h + bottom - 14)}"/>')
        cursor = base
        filled = [(cls, v) for cls, v in item["parts"] if v > 0]
        for k, (cls, value) in enumerate(filled):
            h = value / scale_top * plot_h
            y = cursor - h
            # the 2px surface gap eats into the bottom of every stacked segment
            drawn_h = h if k == 0 else max(h - 2, 1)
            rounded = 4 if k == len(filled) - 1 else 0
            out.append(f'<path class="bar {cls}" d="{_top_rounded(x, y, bar_w, drawn_h, rounded)}"/>')
            cursor = y
        if cap_labels and total > 0:
            out.append(f'<text class="value" x="{num(x + bar_w / 2)}" y="{num(cursor - 6)}" '
                       f'text-anchor="middle">{fmt_count(total)}</text>')
        if i % skip == 0:
            out.append(f'<text class="tick" x="{num(x + bar_w / 2)}" y="{base + 16}" '
                       f'text-anchor="middle">{esc(item["label"])}</text>')
        out.append("</g>")
    if x_title:
        out.append(f'<text class="tick" x="{left + plot_w / 2}" y="{height - 6}" '
                   f'text-anchor="middle">{esc(x_title)}</text>')
    out.append("</svg>")
    return "".join(out)
