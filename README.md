# CPB Analyzer
Codon Pair Bias Analyzer
"""
Codon Pair Bias Analyzer — Day 1 Draft
========================================
Project stage : Day 1 (Research & Planning)
Owner         : Albey Mathew
Scope today   : Understand codon pairs, decide how the analyzer will work,
                and write the *basic* Python logic (sequence -> codons ->
                codon pairs -> counts). Full statistical bias scoring,
                the interface, and integration are LATER days — not built
                here on purpose.

--------------------------------------------------------------------------
NOTES: What is a codon and a codon pair?
--------------------------------------------------------------------------
- A CODON is a group of 3 consecutive nucleotides (A, T/U, G, C) in a DNA
  or mRNA sequence. Each codon usually codes for one amino acid (or a
  stop signal). Example: "ATG" -> codes for Methionine (Met/Start).

- A CODON PAIR is two consecutive, in-frame codons read together.
  Example sequence:  ATG GCT TTA CGA
  Reading frame splits it into codons: [ATG, GCT, TTA, CGA]
  The codon PAIRS (overlapping, one step at a time) are:
      (ATG, GCT), (GCT, TTA), (TTA, CGA)
  i.e. for N codons there are (N - 1) codon pairs.

--------------------------------------------------------------------------
NOTES: Why codon pair analysis matters (background for the project)
--------------------------------------------------------------------------
- Not all codon pairs occur in a genome as often as random chance would
  predict. Some pairs are "over-represented" and some are
  "under-represented" — this non-random usage pattern is called
  CODON PAIR BIAS.
- It matters in bioinformatics/synthetic biology because:
    1. It affects translation speed and accuracy (ribosome efficiency).
    2. It's used in "codon pair deoptimization" — a technique to
       intentionally weaken viruses for vaccine design (e.g. SAVE method).
    3. It can hint at evolutionary pressure on a gene/organism.
- Full bias analysis usually compares OBSERVED codon pair frequency vs
  EXPECTED frequency (based on individual codon frequencies), producing
  a "Codon Pair Score" (CPS). That statistical scoring is planned for a
  LATER day — today we only build the counting foundation it depends on.

--------------------------------------------------------------------------
NOTES: Planned analyzer design (decided today, Day 1)
--------------------------------------------------------------------------
Pipeline we've agreed on:
    1. Take a raw DNA sequence as input (string, or later a FASTA file
       via Biopython's SeqIO).
    2. Clean/validate it (uppercase, check only A/T/G/C, check length is
       a multiple of 3).
    3. Split it into codons (non-overlapping triplets, respecting the
       reading frame).
    4. Generate codon PAIRS from consecutive codons (overlapping by one
       codon, as shown above).
    5. Count how often each codon pair occurs -> frequency table.
    6. (Day 2/3, NOT done here) Compare observed vs expected frequency,
       compute bias score, build the interface, connect to the website.

We chose Biopython's Seq object for validation/handling because it's the
standard bioinformatics library and will make FASTA file support (real
gene sequences) trivial to add later, instead of writing our own parser.

--------------------------------------------------------------------------
"""

from collections import Counter
from itertools import product

try:
    from Bio.Seq import Seq          # Biopython — standard bio sequence handling
    BIOPYTHON_AVAILABLE = True
except ImportError:
    BIOPYTHON_AVAILABLE = False       # Falls back to plain string logic if
                                       # Biopython isn't installed yet.


VALID_BASES = set("ATGC")


def clean_and_validate_sequence(raw_sequence: str) -> str:
    """
    Step 2 of the plan: clean + validate the input DNA sequence.

    - Uppercases the sequence.
    - Strips whitespace/newlines (useful later for pasted FASTA bodies).
    - Confirms only A, T, G, C are present (basic DNA alphabet check).
    - Confirms the length is a multiple of 3 (so it splits into whole
      codons with no leftover bases).

    Raises a ValueError with a clear message if anything is invalid —
    this keeps error-handling logic (Day 3 task) easy to extend later.
    """
    sequence = raw_sequence.strip().upper().replace("\n", "").replace(" ", "")

    invalid_chars = set(sequence) - VALID_BASES
    if invalid_chars:
        raise ValueError(
            f"Sequence contains invalid base(s): {sorted(invalid_chars)}. "
            f"Only A, T, G, C are allowed."
        )

    if len(sequence) % 3 != 0:
        raise ValueError(
            f"Sequence length ({len(sequence)}) is not a multiple of 3 — "
            f"it can't be split into whole codons."
        )

    if BIOPYTHON_AVAILABLE:
        # Round-trip through Biopython's Seq object as our validated form.
        sequence = str(Seq(sequence))

    return sequence


def split_into_codons(sequence: str) -> list[str]:
    """
    Step 3 of the plan: split a validated sequence into non-overlapping
    triplets (codons), respecting the reading frame from position 0.
    """
    return [sequence[i:i + 3] for i in range(0, len(sequence), 3)]


def get_codon_pairs(codons: list[str]) -> list[tuple[str, str]]:
    """
    Step 4 of the plan: build consecutive, overlapping codon pairs.
    For codons [c1, c2, c3, c4] this returns:
        [(c1, c2), (c2, c3), (c3, c4)]
    i.e. len(codons) - 1 pairs.
    """
    return [(codons[i], codons[i + 1]) for i in range(len(codons) - 1)]


def count_codon_pairs(codon_pairs: list[tuple[str, str]]) -> Counter:
    """
    Step 5 of the plan: frequency count of each codon pair.
    Returns a Counter mapping (codonA, codonB) -> occurrences.
    This table is the foundation the Day-2/3 bias-score math will use.
    """
    return Counter(codon_pairs)


def analyze_sequence(raw_sequence: str) -> dict:
    """
    Ties steps 2-5 together into one call, so this can be imported and
    reused as-is once the interface (Arushi's Day 2 task) is ready to
    call into it.
    """
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


# --------------------------------------------------------------------------
# Quick manual demo (Day 1 sanity check only — no real gene, no bias score)
# --------------------------------------------------------------------------
if __name__ == "__main__":
    sample_sequence = "ATGGCTTTACGAATGGCTTAAGGC"  # made-up test sequence

    result = analyze_sequence(sample_sequence)

    print("Input sequence :", sample_sequence)
    print("Sequence length:", result["sequence_length"])
    print("Codons         :", result["codons"])
    print("Codon pairs    :", result["codon_pair_count"])
    print("\nCodon pair frequency table:")
    for pair, count in result["codon_pair_frequencies"].items():
        print(f"  {pair[0]}-{pair[1]} : {count}")
