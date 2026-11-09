"""Write results to disk: TSV table, FASTA files and a JSON summary."""

import json
from pathlib import Path

from . import __version__

TSV_COLUMNS = (
    "orf_id", "seq_id", "strand", "frame", "start", "end", "nt_length",
    "aa_length", "start_codon", "stop_codon", "partial", "gc_percent",
    "mass_kda", "expect", "support",
)


def _open(path):
    return open(path, "w", encoding="utf-8", newline="\n")


def tsv_row(orf):
    return (
        orf.orf_id, orf.seq_id, orf.strand, orf.frame, orf.start, orf.end,
        orf.nt_length, orf.aa_length, orf.start_codon, orf.stop_codon or "-",
        "yes" if orf.partial else "no", f"{orf.gc:.1f}", f"{orf.mass_kda:.2f}",
        format(orf.expect, ".3g"), orf.support,
    )


def write_tsv(path, analysis):
    with _open(path) as out:
        out.write("\t".join(TSV_COLUMNS) + "\n")
        for orf in analysis.orfs:
            out.write("\t".join(str(v) for v in tsv_row(orf)) + "\n")


def _fasta(out, header, seq, width=60):
    out.write(f">{header}\n")
    for i in range(0, len(seq), width):
        out.write(seq[i:i + width] + "\n")


def _header(orf, extra):
    return (f"{orf.orf_id} {orf.strand}{orf.frame} {orf.start}-{orf.end} "
            f"{extra}{' partial' if orf.partial else ''}")


def write_proteins(path, analysis):
    with _open(path) as out:
        for orf in analysis.orfs:
            _fasta(out, _header(orf, f"{orf.aa_length}aa"), orf.protein)


def write_nucleotides(path, analysis):
    with _open(path) as out:
        for result in analysis.results:
            for orf in result.orfs:
                _fasta(out, _header(orf, f"{orf.nt_length}nt"), orf.nucleotides(result.record.seq))


def write_json(path, analysis, command):
    params = analysis.params
    data = {
        "tool": "orf",
        "version": __version__,
        "command": command,
        "parameters": {
            "genetic_code": params.code.id,
            "strands": params.strands,
            "frames": list(params.frames),
            "min_aa": params.min_aa,
            "start_codons": list(params.starts),
            "nested": params.nested,
            "partial": params.partial,
        },
        "sequences": [
            {
                "id": r.record.id,
                "length": r.record.length,
                "gc_percent": round(r.gc, 2),
                "ambiguous_bases": r.ambiguous,
                "orf_count": len(r.orfs),
                "coverage_percent": round(r.coverage, 2),
            }
            for r in analysis.results
        ],
        "orfs": [
            dict(zip(TSV_COLUMNS, tsv_row(orf))) | {"protein": orf.protein}
            for orf in analysis.orfs
        ],
    }
    with _open(path) as out:
        json.dump(data, out, indent=2)
        out.write("\n")


def write_all(outdir, analysis, command, with_json=False):
    """Write every table/FASTA file and return {filename: description}."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    write_tsv(outdir / "orfs.tsv", analysis)
    write_proteins(outdir / "proteins.faa", analysis)
    write_nucleotides(outdir / "orfs.fna", analysis)
    written = {
        "orfs.tsv": "one row per ORF (open in Excel or pandas)",
        "proteins.faa": "translated proteins, FASTA",
        "orfs.fna": "ORF nucleotide sequences, FASTA",
    }
    if with_json:
        write_json(outdir / "summary.json", analysis, command)
        written["summary.json"] = "everything above in one JSON file"
    return written
