"""FASTA parsing and sequence record types."""

from dataclasses import dataclass


@dataclass(frozen=True)
class SequenceRecord:
    name: str
    sequence: str


def parse_fasta(text: str) -> list[SequenceRecord]:
    """Parse a FASTA-formatted string into a list of SequenceRecords.

    Also accepts bare (headerless) sequences — they are assigned the name
    'sequence_1', 'sequence_2', etc.
    """
    records: list[SequenceRecord] = []
    name: str | None = None
    parts: list[str] = []
    bare_index = 0

    lines = text.splitlines()

    # Detect headerless input: no line starts with '>'
    has_header = any(line.lstrip().startswith(">") for line in lines)

    if not has_header:
        # Treat the whole block as one or more whitespace-separated sequences
        sequence = "".join(lines).replace(" ", "").replace("\t", "")
        if sequence:
            bare_index += 1
            records.append(SequenceRecord(name=f"sequence_{bare_index}", sequence=sequence))
        return records

    for line in lines:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name is not None:
                records.append(SequenceRecord(name=name, sequence="".join(parts)))
            # Use only the first token as the name; rest is description
            name = line[1:].split()[0] if line[1:].strip() else f"sequence_{bare_index + 1}"
            parts = []
        else:
            parts.append(line)

    if name is not None:
        records.append(SequenceRecord(name=name, sequence="".join(parts)))

    return records
