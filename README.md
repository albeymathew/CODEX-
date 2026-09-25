# CPB Analyzer
Codon Pair Bias Analyzer

## Overview
Codon pair bias analysis studies how frequently adjacent codons occur together in a sequence compared with what would be expected from individual codon usage. This project is focused on building a practical analysis pipeline for DNA sequence data, starting with validation and counting logic and then expanding toward bias scoring and reporting.

## Project stage
Day 1 (Research & Planning)

## Owner
Albey Mathew

## Scope today
Understand codon pairs, decide how the analyzer will work, and write the basic Python logic: sequence -> codons -> codon pairs -> counts. Full statistical bias scoring, the interface, and integration are planned for later days.

## Notes: What is a codon and a codon pair?
- A CODON is a group of 3 consecutive nucleotides (A, T/U, G, C) in a DNA or mRNA sequence.
- A CODON PAIR is two consecutive, in-frame codons read together.
- Example sequence: ATG GCT TTA CGA
- Codon pairs: (ATG, GCT), (GCT, TTA), (TTA, CGA)

## Why codon pair analysis matters
- It can reveal over- or under-representation of codon combinations.
- It affects translation speed and accuracy.
- It is used in codon pair deoptimization and evolutionary studies.
- Full bias analysis compares observed versus expected codon pair frequency.

## Planned analyzer design
1. Take a raw DNA sequence as input.
2. Clean and validate it.
3. Split it into codons.
4. Generate codon pairs from consecutive codons.
5. Count how often each pair occurs.
6. Later compute a bias score and expose the results through an interface.

We chose Biopython's Seq object for validation and handling because it is the standard bioinformatics library and will make FASTA support easier later.

## Example implementation
```python
from collections import Counter

VALID_BASES = set("ATGC")


def clean_and_validate_sequence(raw_sequence: str) -> str:
    sequence = raw_sequence.strip().upper().replace("\n", "").replace(" ", "")

    invalid_chars = set(sequence) - VALID_BASES
    if invalid_chars:
        raise ValueError(f"Sequence contains invalid base(s): {sorted(invalid_chars)}")

    if len(sequence) % 3 != 0:
        raise ValueError(f"Sequence length ({len(sequence)}) is not a multiple of 3")

    return sequence


def split_into_codons(sequence: str) -> list[str]:
    return [sequence[i:i + 3] for i in range(0, len(sequence), 3)]


def get_codon_pairs(codons: list[str]) -> list[tuple[str, str]]:
    return [(codons[i], codons[i + 1)] for i in range(len(codons) - 1)]


def count_codon_pairs(codon_pairs: list[tuple[str, str]]) -> Counter:
    return Counter(codon_pairs)


def analyze_sequence(raw_sequence: str) -> dict:
    sequence = clean_and_validate_sequence(raw_sequence)
    codons = split_into_codons(sequence)
    pairs = get_codon_pairs(codons)
    pair_counts = count_codon_pairs(pairs)

    return {
        "sequence_length": len(sequence),
        "codon_count": len(codons),
        "codons": codons,
        "codon_pair_count": len(pairs),
        "codon_pair_frequencies": dict(pair_counts),
    }
```

## Daily Reports
- [Second Day Report](day-2-report.md)

## Project Code Overview

### Goal
Build a Codon Pair Bias Analyzer that reads sequence data, extracts codons, computes pairwise codon bias patterns, and reports the results in a structured and interpretable format.

### Core Code Components
- Sequence input loader
- Codon extraction and validation
- Codon frequency calculation
- Pair-bias metric computation
- Summary and reporting utilities

### Planned Implementation Structure
```python
# Example conceptual structure
class SequenceAnalyzer:
    def __init__(self, sequence):
        self.sequence = sequence

    def extract_codons(self):
        return [self.sequence[i:i+3] for i in range(0, len(self.sequence), 3)]

    def compute_codon_frequencies(self):
        pass

    def compute_pair_bias(self):
        pass

    def generate_report(self):
        pass
```

## Second Day Report

### Summary
The second day focused on clarifying the project scope, identifying the core analytical requirements, and outlining a practical implementation plan for the codon-pair bias analysis pipeline.

### Work Completed
- Reviewed the repository baseline and confirmed the project is currently in an early stage.
- Re-examined the project goal: analyze codon pair bias patterns in biological sequence data.
- Defined the key functional areas for the tool:
  - input sequence parsing
  - codon extraction and frequency analysis
  - pairwise codon bias calculations
  - output summaries and visualization support
- Identified likely output artifacts such as codon usage tables and bias score summaries.
- Drafted the implementation direction for a lightweight, modular design.

### Key Findings
- The project requires a clean separation between data ingestion, statistical calculations, and reporting.
- Codon pair bias analysis is best organized around deterministic metrics that can be validated.
- A modular structure will make the analyzer easier to test and scale.

### Next Steps
1. Define the exact input/output data model for sequence analysis.
2. Implement a prototype for codon extraction and frequency counting.
3. Add pair-bias scoring logic and sanity checks against sample inputs.
4. Build a basic reporting layer for summarizing the computed bias.
5. Validate the workflow with representative sequence examples.

### Overall Status
The project is progressing as an early-stage design and prototype effort. The core requirements are now clearer, and the next phase should focus on implementation and validation.
