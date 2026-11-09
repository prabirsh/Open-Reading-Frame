# Changelog

## 1.0.0

The first proper release. The single-file script became a small command line tool.

- `orf find`: scan both strands and all frames for ORFs, with minimum length, start codons,
  nested starts, partial ORFs and three genetic codes.
- Outputs: `orfs.tsv`, `proteins.faa`, `orfs.fna`, optional `summary.json`.
- A self-contained HTML report: six-frame ORF map, GC track, ORFs per frame, length histogram,
  amino-acid composition and a sortable ORF table. Light and dark themes, hover details.
- A "support" score that compares each ORF's length with what chance would produce.
- `orf translate` for plain frame translation, `orf demo` for a first run with no input.
- Reads FASTA, plain text, gzip and stdin; handles lowercase, RNA and IUPAC codes.
- No dependencies. Python 3.9 or newer.

The 0.x script is still in the git history (`orf.py`, commit `d1bbc5b`).
