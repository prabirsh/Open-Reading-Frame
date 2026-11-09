"""Run the finder over records and gather everything the outputs need."""

import time
from dataclasses import dataclass

from . import stats
from .finder import find_orfs


@dataclass
class RecordResult:
    record: object
    orfs: list
    gc: float
    gc_window: int
    gc_points: list
    ambiguous: int
    coverage: float      # % of the sequence covered by at least one ORF

    @property
    def longest(self):
        return max(self.orfs, key=lambda o: o.aa_length, default=None)


@dataclass
class Analysis:
    params: object
    results: list
    elapsed: float = 0.0

    @property
    def orfs(self):
        return [orf for result in self.results for orf in result.orfs]

    @property
    def total_length(self):
        return sum(r.record.length for r in self.results)

    @property
    def gc(self):
        weighted = sum(r.gc * r.record.length for r in self.results)
        return weighted / self.total_length if self.total_length else 0.0

    @property
    def longest(self):
        return max(self.orfs, key=lambda o: o.aa_length, default=None)


def _coverage(orfs, length):
    if not orfs or not length:
        return 0.0
    covered = 0
    run_start = run_end = None
    for orf in sorted(orfs, key=lambda o: o.start):
        if run_end is None or orf.start > run_end + 1:
            if run_end is not None:
                covered += run_end - run_start + 1
            run_start, run_end = orf.start, orf.end
        else:
            run_end = max(run_end, orf.end)
    covered += run_end - run_start + 1
    return 100.0 * covered / length


def _annotate(orfs, record, params):
    """Fill in GC, mass and the chance-expectation for each ORF."""
    forward = stats.base_frequencies(record.seq)
    stop_prob = {
        "+": stats.stop_probability(forward, params.code.stops),
        "-": stats.stop_probability(stats.flip_strand(forward), params.code.stops),
    }
    for orf in orfs:
        orf.gc = stats.gc_percent(orf.nucleotides(record.seq))
        orf.mass_kda = stats.protein_mass_kda(orf.protein)
        orf.expect = sum(
            stats.expected_by_chance(record.length, orf.aa_length, stop_prob[s], len(params.frames))
            for s in params.strand_symbols)
        orf.support = stats.support_label(orf.expect)


def analyze_record(record, params, gc_window=None):
    orfs = find_orfs(record, params)
    _annotate(orfs, record, params)
    window = gc_window or stats.auto_window(record.length)
    return RecordResult(
        record=record,
        orfs=orfs,
        gc=stats.gc_percent(record.seq),
        gc_window=window,
        gc_points=stats.gc_windows(record.seq, window),
        ambiguous=record.length - sum(record.seq.count(b) for b in "ACGT"),
        coverage=_coverage(orfs, record.length),
    )


def analyze(records, params, gc_window=None):
    began = time.perf_counter()
    results = [analyze_record(r, params, gc_window) for r in records]
    return Analysis(params, results, elapsed=time.perf_counter() - began)
