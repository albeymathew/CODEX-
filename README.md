# CPB Analyzer

Codon Pair Bias Analyzer is a command-line tool for inspecting coding DNA or
RNA sequences. It validates and parses FASTA or plain sequence input, counts
codons and adjacent in-frame codon pairs, calculates relative synonymous codon
usage (RSCU), and reports pair enrichment scores.

## Project 3 — Final implementation

**Owners:** Albey Mathew and Arushi Tiwary

The final project is a complete, dependency-free Python application with:

- FASTA (single or multiple records) and plain sequence input
- DNA and RNA input (`U` is normalized to `T`)
- nucleotide validation, selectable reading frame, and explicit handling of
  trailing incomplete codons
- per-record and pooled codon usage, RSCU, and adjacent pair statistics
- human-readable, JSON, CSV, and TSV reports
- tests for input parsing, metrics, aggregation, and the CLI

## Install and run

Requires Python 3.10 or newer. Install from the repository:

```bash
python -m pip install .
```

Analyze a file:

```bash
cpb-analyzer sequences.fasta
```

Or run without installing the command:

```bash
python -m cpb_analyzer sequences.fasta --format json
```

Read a sequence from standard input:

```bash
printf 'ATG GCT TTA CGA\n' | cpb-analyzer -
```

Write a report to a file:

```bash
cpb-analyzer sequences.fasta --format csv --output codon-pairs.csv
```

Run the test suite:

```bash
python -m unittest discover -s tests -v
```

## Command options

```text
cpb-analyzer [input] [-o OUTPUT] [--format {text,json,csv,tsv}]
             [--frame {0,1,2}] [--partial-codon {error,trim}]
             [--include-zero-pairs] [--top-pairs N]
```

- `input`: FASTA or plain sequence path; defaults to standard input (`-`).
- `--frame`: number of leading nucleotides to skip before codonizing; defaults
  to `0`.
- `--partial-codon error|trim`: reject incomplete trailing codons by default,
  or discard and report the trailing bases.
- `--format`: text and JSON include sequence summaries and codon usage; CSV and
  TSV export pair metrics for each record and the pooled aggregate.
- `--include-zero-pairs`: include the complete 64 × 64 pair table in JSON,
  CSV, or TSV output. Without this option, only observed pairs are reported.
- `--top-pairs`: number of highest and lowest scored observed pairs shown in
  text output.

Only unambiguous `A`, `C`, `G`, `T`, and `U` bases are accepted. Ambiguous
IUPAC symbols, gaps, and non-nucleotide characters are rejected with an error.
FASTA records are analyzed independently; pairs never cross record boundaries.
The aggregate combines their counts and within-record pairs.

## Metrics and interpretation

For codons \(i\) and \(j\), the expected pair count is calculated under an
independence model:

```text
expected_count(i,j) = number_of_within-record_pairs × frequency(i) × frequency(j)
observed_expected_ratio(i,j) = observed_count(i,j) / expected_count(i,j)
log2_score(i,j) = log2(observed_expected_ratio(i,j))
```

The aggregate model uses codon frequencies pooled over all analyzed records.
A score of `0` matches the model expectation; positive values indicate
enrichment and negative values indicate depletion. The report includes the
expected count and ratio alongside the score. A pair with zero observed count
has no finite log2 score, and pairs whose expected count is zero have no ratio.
These descriptive scores are not significance tests and do not model
sequence-specific dependencies or biological covariates.

RSCU is the observed count of a codon divided by the average count of codons
encoding the same amino acid. The standard genetic code is used; stop codons
are counted as a synonymous group of one. RSCU is `null` in JSON when no
codon in that synonymous group occurs. The three stop codons are grouped as
termination codons for RSCU. The aggregate JSON object reports `codons` as
`null`; its pooled codon counts and frequencies are present, while ordered
codons remain available in each individual record.

## Repository layout

```text
cpb_analyzer/    validation, statistics, input parsing, and report rendering
tests/           unit and CLI tests
day-2-report.md  project planning and progress from day 2
day-3-report.md  final implementation and validation from day 3
```

## Project history

- [Second Day Report](day-2-report.md)
- [Third Day Report](day-3-report.md)
