"""Reading DNA sequences from FASTA files, plain text files, gzip or stdin."""

import gzip
import re
import sys
from dataclasses import dataclass
from pathlib import Path

IUPAC = "ACGTRYSWKMBDHVN"
_INVALID = re.compile(f"[^{IUPAC}]")
_COMPLEMENT = str.maketrans(IUPAC, "TGCAYRSWMKVHDBN")


class SequenceError(ValueError):
    """The input could not be read or is not a DNA sequence."""


@dataclass
class Record:
    id: str
    description: str
    seq: str

    @property
    def length(self):
        return len(self.seq)


def reverse_complement(seq):
    return seq.translate(_COMPLEMENT)[::-1]


def load(source):
    """Read records from a path, or from stdin when `source` is '-'."""
    if source == "-":
        text = sys.stdin.read()
        default_id = "stdin"
    else:
        text = _read_text(Path(source))
        default_id = Path(source).name.split(".")[0] or "sequence"
    records = parse(text, default_id)
    if not records:
        raise SequenceError(f"no sequence found in {source}")
    return records


def _read_text(path):
    if path.is_dir():   # Windows raises PermissionError here, so look first
        raise SequenceError(f"{path} is a directory, not a sequence file")
    opener = gzip.open if path.suffix == ".gz" else open
    try:
        with opener(path, "rt", encoding="utf-8-sig") as handle:
            return handle.read()
    except FileNotFoundError:
        raise SequenceError(f"file not found: {path}") from None
    except UnicodeDecodeError:
        raise SequenceError(f"{path} is not a text file") from None
    except OSError as err:
        raise SequenceError(f"could not read {path}: {err}") from None


def parse(text, default_id="sequence"):
    """Parse FASTA, or a bare sequence with no header line."""
    records = []
    seen = set()
    header = None
    chunks = []

    def finish():
        if header is None and not chunks:
            return
        parts = (header or "").split(None, 1)
        if parts:
            name = parts[0]
        elif header is None:
            name = default_id
        else:
            name = f"seq{len(records) + 1}"
        description = parts[1] if len(parts) > 1 else ""
        seq = "".join("".join(chunks).split()).upper().replace("U", "T")
        if not seq:
            raise SequenceError(f"record '{name}' has no sequence")
        bad = _INVALID.search(seq)
        if bad:
            raise SequenceError(
                f"record '{name}' has an unexpected character {bad.group()!r} "
                f"at position {bad.start() + 1} (expected A, C, G, T, U or IUPAC codes)")
        unique, n = name, 1
        while unique in seen:
            n += 1
            unique = f"{name}_{n}"
        seen.add(unique)
        records.append(Record(unique, description, seq))

    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            finish()
            header, chunks = line[1:].strip(), []
        else:
            chunks.append(line)
    finish()
    return records
