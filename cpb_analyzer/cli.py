"""Command-line interface for CPB Analyzer."""

import argparse
import sys
from pathlib import Path
from typing import Sequence

from .core import analyze_records
from .io import parse_sequences
from .report import render_json, render_text, write_pair_table


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="cpb-analyzer",
        description="Analyze codon usage and codon-pair bias in FASTA or raw sequences.",
        epilog="Use '-' as the input path to read a sequence from standard input.",
    )
    parser.add_argument(
        "input",
        nargs="?",
        default="-",
        help="input FASTA or plain sequence file (default: standard input)",
    )
    parser.add_argument(
        "-o", "--output",
        help="write the report to this path instead of standard output",
    )
    parser.add_argument(
        "--format",
        choices=("text", "json", "csv", "tsv"),
        default="text",
        help="report format; CSV/TSV contain pair metrics (default: text)",
    )
    parser.add_argument(
        "--frame",
        type=int,
        choices=(0, 1, 2),
        default=0,
        help="number of leading nucleotides to skip (default: 0)",
    )
    parser.add_argument(
        "--partial-codon",
        choices=("error", "trim"),
        default="error",
        help="reject or trim trailing incomplete codon(s) (default: error)",
    )
    parser.add_argument(
        "--include-zero-pairs",
        action="store_true",
        help="include all 4096 codon pairs, including pairs not observed",
    )
    parser.add_argument(
        "--top-pairs",
        type=int,
        default=10,
        metavar="N",
        help="number of enriched/depleted pairs shown in text reports (default: 10)",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.top_pairs < 0:
        parser.error("--top-pairs must be zero or greater")

    try:
        if args.input == "-":
            text = sys.stdin.read()
            source_name = "<stdin>"
        else:
            input_path = Path(args.input)
            text = input_path.read_text(encoding="utf-8")
            source_name = input_path.name
        records = parse_sequences(text, source_name)
        report = analyze_records(
            records,
            frame=args.frame,
            partial_codon=args.partial_codon,
            include_zero_pairs=args.include_zero_pairs,
        )
    except (OSError, UnicodeError, ValueError) as error:
        parser.error(str(error))

    if args.format == "text":
        output = render_text(report, args.top_pairs)
    elif args.format == "json":
        output = render_json(report)
    else:
        from io import StringIO

        stream = StringIO()
        write_pair_table(report, stream, delimiter="," if args.format == "csv" else "\t")
        output = stream.getvalue()

    try:
        if args.output:
            Path(args.output).write_text(output, encoding="utf-8")
        else:
            sys.stdout.write(output)
    except OSError as error:
        parser.error(f"could not write report: {error}")
    return 0
