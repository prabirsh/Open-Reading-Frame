"""Command line interface: `orf find`, `orf translate` and `orf demo`."""

import argparse
import os
import shlex
import shutil
import sys
import webbrowser
from pathlib import Path

from . import __version__, exporters, report
from .analysis import analyze
from .codes import CODES, get_code
from .finder import STRANDS, SearchParams
from .seqio import SequenceError, load, parse, reverse_complement
from .term import Console

DEMO_FASTA = Path(__file__).resolve().parent / "data" / "demo.fasta"

TOP_DESCRIPTION = """\
Find open reading frames (ORFs) in DNA and see what you've got: a table,
FASTA files with the proteins, and a self-contained HTML report with charts.
"""

TOP_EPILOG = """\
first time? run the demo, no input file needed:

  orf demo --open

then try your own sequence:

  orf find my_sequence.fasta --open

Run 'orf <command> --help' to see every option of a command.
"""

FIND_EPILOG = """\
examples:
  orf find genome.fasta                         scan all six frames, ORFs of 30 aa and up
  orf find genome.fasta -m 100 --open           long ORFs only, open the report when done
  orf find reads.fa -s forward -f 1             forward strand, frame 1 only
  orf find plasmid.fa -c 11 --starts table      bacterial code with its alternative starts
  cat seq.txt | orf find - --partial            read a pipe, keep ORFs cut off by the end
  orf find --seq ATGAAATAG -m 1                 try a short sequence straight from the prompt

what you get in the output folder:
  report.html    charts and tables, open this first
  orfs.tsv       one row per ORF
  proteins.faa   the translated proteins
  orfs.fna       the ORF nucleotide sequences

Positions are 1-based and always counted on the sequence you gave, even for
reverse-strand ORFs. Frames +1..+3 are on your strand, -1..-3 on the reverse.
"""

TRANSLATE_EPILOG = """\
examples:
  orf translate gene.fasta                      all six frames, FASTA on stdout
  orf translate gene.fasta -f 2 -s forward      just frame +2
  orf translate gene.fasta --to-stop -o out.faa stop at the first stop codon, write a file
  orf translate gene.fasta --rna                print the RNA instead of protein
"""


class HelpFormatter(argparse.RawDescriptionHelpFormatter):
    def __init__(self, prog):
        width = min(shutil.get_terminal_size((100, 24)).columns, 100)
        super().__init__(prog, max_help_position=30, width=width)

    def _format_action(self, action):
        # list the commands without argparse's "{find,translate,demo}" header line
        try:
            if isinstance(action, argparse._SubParsersAction):
                self._indent()
                lines = [super(HelpFormatter, self)._format_action(a) for a in action._get_subactions()]
                self._dedent()
                return "".join(lines)
        except AttributeError:
            pass
        return super()._format_action(action)


# --- argument types ----------------------------------------------------------------

def positive_int(text):
    try:
        value = int(text)
    except ValueError:
        raise argparse.ArgumentTypeError(f"{text!r} is not a whole number") from None
    if value < 1:
        raise argparse.ArgumentTypeError("must be 1 or more")
    return value


def frame_list(text):
    try:
        frames = tuple(int(part) for part in text.replace(" ", "").split(",") if part)
    except ValueError:
        frames = ()
    if not frames or any(f not in (1, 2, 3) for f in frames):
        raise argparse.ArgumentTypeError(f"{text!r}: use frames 1, 2 and/or 3, for example 1,3")
    return frames


def _code_help():
    return ", ".join(f"{c.id} {c.name.lower()}" for c in CODES.values())


# --- parser ------------------------------------------------------------------------

def _input_group(parser):
    group = parser.add_argument_group("input")
    group.add_argument("input", nargs="?", metavar="INPUT",
                       help="DNA file, FASTA or plain text (.gz is fine); use - for stdin")
    group.add_argument("--seq", metavar="DNA",
                       help="analyse a sequence given right here instead of a file")


def _scan_group(parser, with_starts):
    group = parser.add_argument_group("search")
    group.add_argument("-c", "--code", type=int, choices=sorted(CODES), default=1, metavar="TABLE",
                       help=f"genetic code (default 1): {_code_help()}")
    group.add_argument("-s", "--strand", choices=STRANDS, default="both", metavar="STRAND",
                       help="strand to scan: both, forward or reverse (default: both)")
    group.add_argument("-f", "--frames", type=frame_list, default=(1, 2, 3), metavar="LIST",
                       help="frames to scan, comma separated (default: 1,2,3)")
    if with_starts:
        group.add_argument("-m", "--min-length", type=positive_int, default=30, metavar="AA",
                           help="shortest ORF to keep, in amino acids (default: 30)")
        group.add_argument("--starts", default="ATG", metavar="LIST",
                           help="start codons, comma separated (default: ATG); 'table' allows every "
                                "start the genetic code knows, such as GTG and TTG in bacteria")
        group.add_argument("--nested", action="store_true",
                           help="report an ORF for every start codon, not just the outermost one")
        group.add_argument("--partial", action="store_true",
                           help="keep ORFs that run off the end of the sequence with no stop codon")


def build_parser():
    parser = argparse.ArgumentParser(
        prog="orf", usage="orf <command> [options]", formatter_class=HelpFormatter,
        description=TOP_DESCRIPTION, epilog=TOP_EPILOG)
    parser.add_argument("-V", "--version", action="version", version=f"orf {__version__}")
    commands = parser.add_subparsers(dest="command", metavar="<command>", title="commands")

    find = commands.add_parser(
        "find", formatter_class=HelpFormatter, epilog=FIND_EPILOG,
        usage="orf find INPUT [options]\n       orf find --seq DNA [options]",
        help="find ORFs and write tables, FASTA and an HTML report",
        description="Scan DNA for open reading frames on both strands and write the results.")
    _input_group(find)
    _scan_group(find, with_starts=True)
    out = find.add_argument_group("output")
    out.add_argument("-o", "--outdir", default="orf_results", metavar="DIR",
                     help="folder for the results, created if needed (default: orf_results)")
    out.add_argument("--no-report", action="store_true", help="skip the HTML report")
    out.add_argument("--json", action="store_true", help="also write summary.json")
    out.add_argument("--top", type=positive_int, default=50, metavar="N",
                     help="ORFs listed in the report table (default: 50)")
    out.add_argument("--gc-window", type=positive_int, metavar="BP",
                     help="window for the GC track (default: picked from the sequence length)")
    out.add_argument("--title", help="heading for the report")
    out.add_argument("--open", action="store_true", help="open the report in your browser when done")
    out.add_argument("-q", "--quiet", action="store_true", help="print nothing except errors")
    out.add_argument("--no-color", action="store_true",
                     help="plain output (the NO_COLOR variable works too)")
    find.set_defaults(func=cmd_find, _parser=find)

    translate = commands.add_parser(
        "translate", formatter_class=HelpFormatter, epilog=TRANSLATE_EPILOG,
        usage="orf translate INPUT [options]\n       orf translate --seq DNA [options]",
        help="translate DNA in the frames you pick and print FASTA",
        description="Translate whole frames, stop codons included, without looking for ORFs.")
    _input_group(translate)
    _scan_group(translate, with_starts=False)
    tout = translate.add_argument_group("output")
    tout.add_argument("--to-stop", action="store_true", help="end each protein at its first stop codon")
    tout.add_argument("--rna", action="store_true", help="print the RNA of each strand instead of protein")
    tout.add_argument("-w", "--width", type=int, default=60, metavar="N",
                      help="line width of the FASTA output, 0 for one line (default: 60)")
    tout.add_argument("-o", "--output", metavar="FILE", help="write to a file instead of stdout")
    translate.set_defaults(func=cmd_translate, _parser=translate)

    demo = commands.add_parser(
        "demo", formatter_class=HelpFormatter, usage="orf demo [options]",
        help="run on the bundled example, no input file needed",
        description="Run 'orf find' on a small synthetic sequence that ships with the tool.")
    demo.add_argument("-o", "--outdir", default="orf_demo", metavar="DIR",
                      help="folder for the results (default: orf_demo)")
    demo.add_argument("--open", action="store_true", help="open the report in your browser when done")
    demo.add_argument("--no-color", action="store_true", help="plain output")
    demo.set_defaults(func=cmd_demo, _parser=demo)
    return parser


# --- shared helpers ------------------------------------------------------------------

def _read_input(args):
    parser = args._parser
    if args.seq is not None and args.input is not None:
        parser.error("give either INPUT or --seq, not both")
    if args.seq is None and args.input is None:
        parser.error("no input given. Pass a FASTA or text file, '-' for stdin, or --seq DNA "
                     "(not sure where to start? try: orf demo)")
    if args.seq is not None:
        records = parse(args.seq, "sequence")
        if not records:
            raise SequenceError("--seq is empty")
        return records, "--seq"
    return load(args.input), args.input


def _search_params(args):
    code = get_code(args.code)
    if args.starts.strip().lower() == "table":
        starts = code.starts
    else:
        starts = tuple(s for s in args.starts.replace(" ", "").split(",") if s)
    try:
        return SearchParams(code, args.strand, args.frames, args.min_length, starts,
                            args.nested, args.partial)
    except ValueError as err:
        args._parser.error(str(err))


def _soften_encoding_errors():
    for stream in (sys.stdout, sys.stderr):
        try:
            stream.reconfigure(errors="replace")
        except (AttributeError, ValueError):
            pass


def _wrap(seq, width):
    if width <= 0:
        return seq
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))


# --- find --------------------------------------------------------------------------

def cmd_find(args, argv):
    console = Console(color=False if args.no_color else None, quiet=args.quiet)
    records, label = _read_input(args)
    params = _search_params(args)
    analysis = analyze(records, params, args.gc_window)
    command = "orf " + shlex.join(argv)

    outdir = Path(args.outdir)
    try:
        written = exporters.write_all(outdir, analysis, command, args.json)
        report_path = outdir / "report.html"
        if not args.no_report:
            title = args.title or ("ORF report" if label in ("-", "--seq") else f"ORF report: {Path(label).name}")
            report.write(report_path, analysis, title=title, top=args.top, command=command,
                         input_label="" if label in ("-", "--seq") else Path(label).name)
            written = {"report.html": "charts and tables, open this first", **written}
    except OSError as err:
        console.error(f"could not write to {outdir}: {err}")
        return 1

    _print_summary(console, analysis, label, outdir, written)
    if args.open and not args.no_report:
        webbrowser.open(report_path.resolve().as_uri())
    return 0


def _print_summary(console, analysis, label, outdir, written):
    c = console
    p = analysis.params
    results = analysis.results
    strands = {"both": "both strands", "forward": "forward strand", "reverse": "reverse strand"}[p.strands]
    c.say()
    c.say(c.paint(f"orf {__version__}", "bold") + c.paint(f"  {c.sep}  open reading frames", "dim"))
    c.say()
    plural = "" if len(results) == 1 else "s"
    c.say(f"  {c.paint('input'.ljust(9), 'dim')}{label}  ({len(results):,} sequence{plural}, "
          f"{analysis.total_length:,} nt, GC {analysis.gc:.1f}%)")
    c.say(f"  {c.paint('search'.ljust(9), 'dim')}code {p.code.id}  {c.sep}  {strands}  {c.sep}  "
          f"frames {','.join(map(str, p.frames))}  {c.sep}  min {p.min_aa} aa  {c.sep}  starts {','.join(p.starts)}")
    c.say()

    orfs = analysis.orfs
    if not orfs:
        c.say(f"  {c.paint('No ORFs found', 'yellow', 'bold')} at {p.min_aa} amino acids or longer.")
        hints = []
        if p.min_aa > 10:
            hints.append("lower --min-length")
        if not p.partial:
            hints.append("add --partial to keep ORFs cut off by the end of the sequence")
        if p.starts == ("ATG",):
            hints.append("try --starts table for alternative start codons")
        if hints:
            c.say("  Try: " + "; ".join(hints) + ".")
    else:
        shown = results[:10]
        rows = []
        for r in shown:
            best = r.longest
            longest = f"{best.aa_length} aa ({best.strand}{best.frame})" if best else "-"
            rows.append([r.record.id[:28], f"{r.record.length:,}", f"{r.gc:.1f}",
                         f"{len(r.orfs):,}", longest, f"{r.coverage:.0f}%"])
        c.table(["sequence", "length", "GC%", "ORFs", "longest", "covered"], rows, "lrrrrr")
        if len(results) > len(shown):
            c.say(c.paint(f"  ... and {len(results) - len(shown)} more (see orfs.tsv)", "dim"))
        c.say()

        lanes = [(s, f) for s in p.strand_symbols for f in p.frames]
        if len(lanes) > 1:
            counts = {lane: 0 for lane in lanes}
            for o in orfs:
                counts[(o.strand, o.frame)] += 1
            biggest = max(counts.values())
            c.say(f"  {c.paint('ORFs per frame', 'bold')}")
            for lane in lanes:
                style = "cyan" if lane[0] == "+" else "yellow"
                bar = c.bar(counts[lane], biggest).ljust(24)
                c.say(f"    {lane[0]}{lane[1]}  {c.paint(bar, style)} {counts[lane]}")
            c.say()

        c.say(f"  {c.paint('Longest ORFs', 'bold')}")
        ranked = sorted(orfs, key=lambda o: -o.aa_length)[:5]
        c.table(["ORF", "frame", "start", "end", "aa", "kDa", "support"],
                [[o.orf_id[:30], f"{o.strand}{o.frame}", f"{o.start:,}", f"{o.end:,}",
                  str(o.aa_length), f"{o.mass_kda:.1f}", o.support] for o in ranked],
                "lrrrrrl", indent=4)
        c.say()

    c.say(f"  {c.paint('Wrote to', 'bold')} {outdir}{os.sep}")
    width = max(len(name) for name in written)
    for name, what in written.items():
        c.say(f"    {name:<{width}}  {c.paint(what, 'dim')}")
    c.say()
    if "report.html" in written:
        c.say(f"  {c.paint('Open', 'green', 'bold')} {Path(outdir) / 'report.html'}")
        c.say()


# --- translate ---------------------------------------------------------------------

def cmd_translate(args, argv):
    records, _ = _read_input(args)
    code = get_code(args.code)
    symbols = {"both": "+-", "forward": "+", "reverse": "-"}[args.strand]
    chunks = []
    for rec in records:
        for strand in symbols:
            seq = rec.seq if strand == "+" else reverse_complement(rec.seq)
            if args.rna:
                chunks.append(f">{rec.id} strand={strand} rna\n{_wrap(seq.replace('T', 'U'), args.width)}")
                continue
            for frame in args.frames:
                protein = code.translate(seq, frame - 1, args.to_stop)
                if protein:
                    chunks.append(f">{rec.id} strand={strand} frame={frame} {len(protein)}aa\n"
                                  f"{_wrap(protein, args.width)}")
    text = "\n".join(chunks) + ("\n" if chunks else "")
    if args.output:
        try:
            with open(args.output, "w", encoding="utf-8", newline="\n") as out:
                out.write(text)
        except OSError as err:
            Console().error(f"could not write {args.output}: {err}")
            return 1
    else:
        sys.stdout.write(text)
    return 0


# --- demo --------------------------------------------------------------------------

def cmd_demo(args, argv):
    if not DEMO_FASTA.exists():
        Console().error(f"the bundled example is missing ({DEMO_FASTA}); reinstall the package")
        return 1
    forwarded = ["find", str(DEMO_FASTA), "--outdir", args.outdir]
    if args.open:
        forwarded.append("--open")
    if args.no_color:
        forwarded.append("--no-color")
    status = main(forwarded)
    if status == 0:
        c = Console(color=False if args.no_color else None)
        c.say(c.paint("  Next steps", "bold"))
        c.say(f"    orf find your_sequence.fasta --open   {c.paint('analyse your own DNA', 'dim')}")
        c.say(f"    orf find --help                       {c.paint('every option, with examples', 'dim')}")
        c.say(f"    {DEMO_FASTA}")
        c.say(f"    {c.paint('the demo input, a good template for the FASTA format', 'dim')}")
        c.say()
    return status


# --- entry point ---------------------------------------------------------------------

def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    _soften_encoding_errors()
    parser = build_parser()
    if not argv:
        parser.print_help()
        return 0
    args = parser.parse_args(argv)
    try:
        return args.func(args, argv)
    except SequenceError as err:
        Console().error(str(err))
        return 1
    except KeyboardInterrupt:
        print(file=sys.stderr)
        return 130
    except BrokenPipeError:
        # e.g. `orf translate x.fa | head`; silence the flush error at exit
        try:
            os.dup2(os.open(os.devnull, os.O_WRONLY), sys.stdout.fileno())
        except (OSError, ValueError):
            pass
        return 0
