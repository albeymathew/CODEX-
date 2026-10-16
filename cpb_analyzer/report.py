"""Human-readable and machine-readable report rendering."""

import csv
import json
import sys
from dataclasses import asdict
from typing import TextIO

from .core import AnalysisReport, PairMetric, SequenceAnalysis

METHODOLOGY = (
    "Expected pair frequency assumes adjacent codons are independent, using "
    "the relevant sequence's codon frequencies (pooled frequencies for the "
    "aggregate). The observed/expected ratio is observed "
    "pair count divided by its expected count; log2 score is log2 of that "
    "ratio. A score of 0 indicates the model expectation, positive scores "
    "indicate enrichment, and negative scores indicate depletion. Zero-count "
    "pairs have no finite log2 score. These descriptive scores are not "
    "statistical significance tests."
)


def report_to_dict(report: AnalysisReport) -> dict[str, object]:
    def analysis_to_dict(analysis: SequenceAnalysis) -> dict[str, object]:
        return {
            "name": analysis.name,
            "sequence_length": analysis.sequence_length,
            "frame": analysis.frame,
            "trimmed_bases": analysis.trimmed_bases,
            "codons": list(analysis.codons) if analysis.codons is not None else None,
            "codon_count": analysis.codon_count,
            "codon_counts": analysis.codon_counts,
            "codon_frequencies": analysis.codon_frequencies,
            "rscu": analysis.rscu,
            "pair_count": analysis.pair_count,
            "pair_metrics": [asdict(metric) for metric in analysis.pair_metrics],
        }

    return {
        "tool": "CPB Analyzer",
        "methodology": METHODOLOGY,
        "options": {
            "frame": report.frame,
            "partial_codon_policy": report.partial_codon_policy,
            "include_zero_pairs": report.include_zero_pairs,
        },
        "records": [analysis_to_dict(analysis) for analysis in report.records],
        "aggregate": analysis_to_dict(report.aggregate),
    }


def render_json(report: AnalysisReport) -> str:
    return json.dumps(report_to_dict(report), indent=2, sort_keys=True) + "\n"


def _render_summary(analysis: SequenceAnalysis, top_pairs: int) -> list[str]:
    lines = [
        f"{analysis.name}:",
        f"  Sequence length: {analysis.sequence_length} nt",
        f"  Codons: {analysis.codon_count}",
        f"  Codon pairs: {analysis.pair_count}",
        f"  Trailing bases trimmed: {analysis.trimmed_bases}",
    ]
    if analysis.codon_counts:
        codon_text = ", ".join(
            f"{codon}={count}" for codon, count in analysis.codon_counts.items()
        )
        lines.append(f"  Codon counts: {codon_text}")

    scored_pairs = [
        metric for metric in analysis.pair_metrics
        if metric.log2_observed_expected is not None
    ]
    highest = sorted(
        (metric for metric in scored_pairs if metric.log2_observed_expected > 0),
        key=lambda metric: (metric.log2_observed_expected, metric.codon_1, metric.codon_2),
        reverse=True,
    )[:top_pairs]
    lowest = sorted(
        (metric for metric in scored_pairs if metric.log2_observed_expected < 0),
        key=lambda metric: (metric.log2_observed_expected, metric.codon_1, metric.codon_2),
    )[:top_pairs]

    def format_pair(metric: PairMetric) -> str:
        score = metric.log2_observed_expected
        return (
            f"{metric.codon_1}-{metric.codon_2} "
            f"(observed={metric.observed_count}, expected={metric.expected_count:.3f}, "
            f"log2 ratio={score:+.3f})"
        )

    if highest:
        lines.append("  Most enriched observed pairs:")
        lines.extend(f"    {format_pair(metric)}" for metric in highest)
    if lowest:
        lines.append("  Most depleted observed pairs:")
        lines.extend(f"    {format_pair(metric)}" for metric in lowest)
    if not highest and not lowest:
        if scored_pairs:
            lines.append("  No observed pair has a positive or negative score.")
        else:
            lines.append("  No codon-pair scores are available (at least two codons are needed).")
    return lines


def render_text(report: AnalysisReport, top_pairs: int = 10) -> str:
    lines = [
        "CPB Analyzer",
        f"Frame: {report.frame}; partial codon policy: {report.partial_codon_policy}",
        METHODOLOGY,
        "",
    ]
    for analysis in report.records:
        lines.extend(_render_summary(analysis, top_pairs))
        lines.append("")
    lines.extend(_render_summary(report.aggregate, top_pairs))
    return "\n".join(lines).rstrip() + "\n"


def write_pair_table(
    report: AnalysisReport,
    stream: TextIO | None = None,
    delimiter: str = ",",
) -> None:
    if stream is None:
        stream = sys.stdout
    writer = csv.writer(stream, delimiter=delimiter, lineterminator="\n")
    writer.writerow(
        [
            "record",
            "codon_1",
            "codon_2",
            "observed_count",
            "expected_count",
            "observed_frequency",
            "expected_probability",
            "observed_expected_ratio",
            "log2_observed_expected",
        ]
    )
    analyses = (*report.records, report.aggregate)
    for analysis in analyses:
        for metric in analysis.pair_metrics:
            writer.writerow(
                [
                    analysis.name,
                    metric.codon_1,
                    metric.codon_2,
                    metric.observed_count,
                    f"{metric.expected_count:.12g}",
                    f"{metric.observed_frequency:.12g}",
                    f"{metric.expected_probability:.12g}",
                    "" if metric.observed_expected_ratio is None else f"{metric.observed_expected_ratio:.12g}",
                    "" if metric.log2_observed_expected is None else f"{metric.log2_observed_expected:.12g}",
                ]
            )
