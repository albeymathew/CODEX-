"""Sequence validation and codon-pair statistics."""

from collections import Counter
from dataclasses import dataclass
from math import log2

from .genetic_code import ALL_CODONS, AMINO_ACID_TO_CODONS
from .io import SequenceRecord

VALID_BASES = frozenset("ATGC")


@dataclass(frozen=True)
class PairMetric:
    codon_1: str
    codon_2: str
    observed_count: int
    expected_count: float
    observed_frequency: float
    expected_probability: float
    observed_expected_ratio: float | None
    log2_observed_expected: float | None


@dataclass(frozen=True)
class SequenceAnalysis:
    name: str
    sequence_length: int
    frame: int
    trimmed_bases: int
    codons: tuple[str, ...] | None
    codon_count: int
    codon_counts: dict[str, int]
    codon_frequencies: dict[str, float]
    rscu: dict[str, float | None]
    pair_count: int
    pair_metrics: tuple[PairMetric, ...]


@dataclass(frozen=True)
class AnalysisReport:
    frame: int
    partial_codon_policy: str
    include_zero_pairs: bool
    records: tuple[SequenceAnalysis, ...]
    aggregate: SequenceAnalysis


def _normalize_sequence(record: SequenceRecord, frame: int, partial_codon: str) -> tuple[str, int]:
    sequence = "".join(record.sequence.split()).upper().replace("U", "T")
    invalid = sorted(set(sequence) - VALID_BASES)
    if invalid:
        raise ValueError(
            f"{record.name}: invalid nucleotide(s) {', '.join(invalid)}; "
            "only A, C, G, T, and U are supported"
        )

    sequence = sequence[frame:]
    remainder = len(sequence) % 3
    if remainder and partial_codon == "error":
        raise ValueError(
            f"{record.name}: {len(sequence)} nucleotide(s) remain after frame {frame}; "
            "sequence length must be divisible by 3 "
            "(use --partial-codon trim to discard the trailing bases)"
        )
    trimmed_bases = remainder
    if remainder:
        sequence = sequence[:-remainder]
    if not sequence:
        raise ValueError(f"{record.name}: no complete codons remain after applying the frame")
    return sequence, trimmed_bases


def _analyze_codons(
    name: str,
    sequence: str,
    frame: int,
    trimmed_bases: int,
    include_zero_pairs: bool,
) -> SequenceAnalysis:
    codons = tuple(sequence[index:index + 3] for index in range(0, len(sequence), 3))
    return _analysis_from_counts(
        name=name,
        codons=codons,
        codon_counts=Counter(codons),
        pair_counts=Counter(zip(codons, codons[1:])),
        sequence_length=len(sequence),
        frame=frame,
        trimmed_bases=trimmed_bases,
        include_zero_pairs=include_zero_pairs,
    )


def _analysis_from_counts(
    *,
    name: str,
    codons: tuple[str, ...] | None,
    codon_counts: Counter[str],
    pair_counts: Counter[tuple[str, str]],
    sequence_length: int,
    frame: int,
    trimmed_bases: int,
    include_zero_pairs: bool,
) -> SequenceAnalysis:
    codon_total = sum(codon_counts.values())
    codon_frequencies = {
        codon: count / codon_total for codon, count in sorted(codon_counts.items())
    }
    rscu: dict[str, float | None] = {}
    for synonymous_codons in AMINO_ACID_TO_CODONS.values():
        synonymous_total = sum(codon_counts.get(codon, 0) for codon in synonymous_codons)
        for codon in synonymous_codons:
            rscu[codon] = (
                codon_counts.get(codon, 0) * len(synonymous_codons) / synonymous_total
                if synonymous_total
                else None
            )

    pair_total = sum(pair_counts.values())
    pair_keys = (
        sorted(pair_counts)
        if not include_zero_pairs
        else [(first, second) for first in ALL_CODONS for second in ALL_CODONS]
    )
    metrics: list[PairMetric] = []
    for first, second in pair_keys:
        observed_count = pair_counts.get((first, second), 0)
        expected_probability = (
            codon_counts.get(first, 0) / codon_total
        ) * (
            codon_counts.get(second, 0) / codon_total
        )
        expected_count = pair_total * expected_probability
        if pair_total == 0 or expected_count == 0:
            ratio = None
            score = None
        else:
            ratio = observed_count / expected_count
            score = log2(ratio) if ratio > 0 else None
        metrics.append(
            PairMetric(
                codon_1=first,
                codon_2=second,
                observed_count=observed_count,
                expected_count=expected_count,
                observed_frequency=observed_count / pair_total if pair_total else 0.0,
                expected_probability=expected_probability,
                observed_expected_ratio=ratio,
                log2_observed_expected=score,
            )
        )

    return SequenceAnalysis(
        name=name,
        sequence_length=sequence_length,
        frame=frame,
        trimmed_bases=trimmed_bases,
        codons=codons,
        codon_count=codon_total,
        codon_counts=dict(sorted(codon_counts.items())),
        codon_frequencies=codon_frequencies,
        rscu=dict(sorted(rscu.items())),
        pair_count=pair_total,
        pair_metrics=tuple(metrics),
    )


def analyze_records(
    records: list[SequenceRecord],
    *,
    frame: int = 0,
    partial_codon: str = "error",
    include_zero_pairs: bool = False,
) -> AnalysisReport:
    if not records:
        raise ValueError("at least one sequence record is required")
    if frame not in (0, 1, 2):
        raise ValueError("frame must be 0, 1, or 2")
    if partial_codon not in ("error", "trim"):
        raise ValueError("partial_codon must be 'error' or 'trim'")

    analyses: list[SequenceAnalysis] = []
    normalized: list[tuple[SequenceRecord, str, int]] = []
    for record in records:
        sequence, trimmed_bases = _normalize_sequence(record, frame, partial_codon)
        normalized.append((record, sequence, trimmed_bases))
        analyses.append(
            _analyze_codons(record.name, sequence, frame, trimmed_bases, include_zero_pairs)
        )

    aggregate_trimmed = sum(trimmed for _, _, trimmed in normalized)
    codon_counts: Counter[str] = Counter()
    pair_counts: Counter[tuple[str, str]] = Counter()
    for analysis in analyses:
        codon_counts.update(analysis.codon_counts)
        codons = analysis.codons or ()
        pair_counts.update(zip(codons, codons[1:]))
    aggregate = _analysis_from_counts(
        name="all records",
        codons=None,
        codon_counts=codon_counts,
        pair_counts=pair_counts,
        sequence_length=sum(analysis.sequence_length for analysis in analyses),
        frame=frame,
        trimmed_bases=aggregate_trimmed,
        include_zero_pairs=include_zero_pairs,
    )

    return AnalysisReport(
        frame=frame,
        partial_codon_policy=partial_codon,
        include_zero_pairs=include_zero_pairs,
        records=tuple(analyses),
        aggregate=aggregate,
    )
