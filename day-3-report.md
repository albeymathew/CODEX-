# Third Day Report

## Project
Codon Pair Bias Analyzer (CPB Analyzer)

## Date
2026-10-16

## Summary
The third day completed the CPB Analyzer as a runnable, dependency-free Python
command-line application. It now accepts FASTA and plain-sequence input,
validates and analyzes coding sequences, and produces human-readable and
machine-readable codon-pair reports.

## Work Completed
- Implemented FASTA parsing for one or more records and plain sequence input.
- Added DNA/RNA normalization (`U` is converted to `T`) and validation for
  unambiguous nucleotide bases.
- Added reading-frame selection and explicit reject-or-trim handling for
  incomplete trailing codons.
- Implemented codon counts, codon frequencies, relative synonymous codon usage
  (RSCU), adjacent codon-pair counts, and observed-versus-expected pair scores.
- Added per-record reports and a pooled aggregate that does not create
  artificial codon pairs between separate FASTA records.
- Added text, JSON, CSV, and TSV output, including an option to include all
  4,096 possible codon pairs.
- Added installation metadata, a `cpb-analyzer` command, and module execution
  through `python -m cpb_analyzer`.
- Documented installation, command options, metrics, assumptions, and example
  usage in the project README.

## Analysis Method
The expected count for a codon pair is calculated under an independence model
from the codon frequencies and the number of within-record adjacent pairs. The
observed-to-expected ratio and its base-2 logarithm describe enrichment or
depletion relative to that model. These are descriptive scores, not statistical
significance tests. RSCU is calculated relative to synonymous codons using the
standard genetic code.

## Validation
- Added unit and CLI tests covering input parsing, validation, DNA/RNA
  normalization, reading frames, incomplete codons, statistics, aggregation,
  and output formats.
- Ran `python -m unittest discover -s tests -v`: all 15 tests passed.
- Checked representative text output and multi-record JSON aggregation.

## How to Run

```bash
python -m pip install .
cpb-analyzer sequences.fasta
```

Alternatively, run directly from the repository:

```bash
python -m cpb_analyzer sequences.fasta --format json
```

## Challenges and Considerations
- Pair counts must not bridge separate FASTA records; aggregate statistics
  combine within-record counts while preserving record boundaries.
- Zero observed or expected counts cannot always produce a finite log2 score;
  those cases are represented without inventing a score.
- Codon-pair scores depend on an independence baseline and should not be
  interpreted as significance tests or as a complete biological model.

## Overall Status
The final CPB Analyzer project is implemented, documented, and tested. The
repository README contains complete setup and usage instructions, and this
report records the Day 3 implementation and validation results.
