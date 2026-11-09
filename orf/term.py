"""Small terminal helpers: colour, tables and bars that degrade gracefully."""

import os
import shutil
import sys

_CODES = {"bold": "1", "dim": "2", "red": "31", "green": "32", "yellow": "33", "cyan": "36"}


def _enable_windows_ansi():
    if os.name != "nt":
        return True
    try:
        import ctypes
        kernel = ctypes.windll.kernel32
        handle = kernel.GetStdHandle(-11)
        mode = ctypes.c_ulong()
        if not kernel.GetConsoleMode(handle, ctypes.byref(mode)):
            return False
        return bool(kernel.SetConsoleMode(handle, mode.value | 0x0004))
    except Exception:
        return False


def _wants_color(stream):
    if os.environ.get("NO_COLOR"):
        return False
    if os.environ.get("FORCE_COLOR"):
        return True
    if not hasattr(stream, "isatty") or not stream.isatty():
        return False
    if os.environ.get("TERM") == "dumb":
        return False
    return _enable_windows_ansi()


def _can_print(stream, text):
    try:
        text.encode(getattr(stream, "encoding", None) or "ascii")
        return True
    except (LookupError, UnicodeEncodeError):
        return False


class Console:
    def __init__(self, color=None, quiet=False, stream=None):
        self.stream = stream or sys.stdout
        self.quiet = quiet
        self.color = _wants_color(self.stream) if color is None else color
        # box-drawing only on a real terminal that can show it; pipes and files get plain ASCII
        fancy = hasattr(self.stream, "isatty") and self.stream.isatty()
        self.block = "█" if fancy and _can_print(self.stream, "█") else "#"
        self.sep = "·" if fancy and _can_print(self.stream, "·") else "|"
        self.width = min(shutil.get_terminal_size((100, 24)).columns, 110)

    def paint(self, text, *styles):
        if not self.color or not styles:
            return text
        start = "".join(f"\033[{_CODES[s]}m" for s in styles)
        return f"{start}{text}\033[0m"

    def say(self, text=""):
        if not self.quiet:
            print(text, file=self.stream)

    def warn(self, text):
        print(f"{self.paint('warning:', 'yellow', 'bold')} {text}", file=sys.stderr)

    def error(self, text):
        print(f"{self.paint('orf: error:', 'red', 'bold')} {text}", file=sys.stderr)

    def table(self, headers, rows, align=None, indent=2):
        """Print rows of plain strings in aligned columns. align is e.g. 'lrr'."""
        align = align or "l" * len(headers)
        widths = [max(len(h), *(len(r[i]) for r in rows)) if rows else len(h)
                  for i, h in enumerate(headers)]

        def fit(cell, i):
            return cell.rjust(widths[i]) if align[i] == "r" else cell.ljust(widths[i])

        pad = " " * indent
        self.say(pad + self.paint("  ".join(fit(h, i) for i, h in enumerate(headers)).rstrip(), "dim"))
        for row in rows:
            self.say(pad + "  ".join(fit(c, i) for i, c in enumerate(row)).rstrip())

    def bar(self, value, biggest, width=24):
        if biggest <= 0 or value <= 0:
            return ""
        return self.block * max(1, round(width * value / biggest))
