"""Input parsing for FASTA files and plain nucleotide sequences."""

from dataclasses import dataclass


@dataclass(frozen=True)
class SequenceRecord:
    name: str
    sequence: str


def parse_sequences(text: str, source_name: str = "<input>") -> list[SequenceRecord]:
    """Parse FASTA records, or treat non-FASTA input as one sequence."""
    lines = text.splitlines()
    is_fasta = any(line.lstrip().startswith(">") for line in lines)
    records: list[SequenceRecord] = []

    if not is_fasta:
        sequence = "".join(text.split())
        if not sequence:
            raise ValueError(f"{source_name}: input contains no sequence")
        return [SequenceRecord(name=source_name, sequence=sequence)]

    current_name: str | None = None
    current_sequence: list[str] = []

    def finish_record(line_number: int) -> None:
        if current_name is None:
            return
        sequence = "".join("".join(current_sequence).split())
        if not sequence:
            raise ValueError(
                f"{source_name}:{line_number}: FASTA record {current_name!r} is empty"
            )
        records.append(SequenceRecord(current_name, sequence))

    for line_number, line in enumerate(lines, start=1):
        stripped = line.strip()
        if not stripped or stripped.startswith(";"):
            continue
        if stripped.startswith(">"):
            finish_record(line_number)
            current_name = stripped[1:].strip()
            if not current_name:
                raise ValueError(f"{source_name}:{line_number}: FASTA header is empty")
            current_sequence = []
        elif current_name is None:
            raise ValueError(
                f"{source_name}:{line_number}: sequence data appears before a FASTA header"
            )
        else:
            current_sequence.append(stripped)

    finish_record(len(lines) + 1)
    if not records:
        raise ValueError(f"{source_name}: FASTA input contains no records")
    return records
