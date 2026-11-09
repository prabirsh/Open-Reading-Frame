import contextlib
import io
from pathlib import Path

from orf.cli import main

ROOT = Path(__file__).resolve().parent.parent
DEMO = ROOT / "orf" / "data" / "demo.fasta"


def run_cli(*args):
    """Run the CLI in-process and return (exit code, stdout, stderr)."""
    out, err = io.StringIO(), io.StringIO()
    code = 0
    with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
        try:
            code = main(list(args))
        except SystemExit as exit_:
            code = exit_.code if isinstance(exit_.code, int) else 1
    return code, out.getvalue(), err.getvalue()
