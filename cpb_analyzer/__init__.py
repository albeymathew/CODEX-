"""Codon Pair Bias Analyzer."""

from .core import AnalysisReport, SequenceRecord, analyze_records
from .io import parse_sequences

__all__ = ["AnalysisReport", "SequenceRecord", "analyze_records", "parse_sequences"]
