import unittest

from cpb_analyzer.core import analyze_records
from cpb_analyzer.io import SequenceRecord, parse_sequences


class InputParsingTests(unittest.TestCase):
    def test_parses_multiline_fasta_records_and_comments(self) -> None:
        records = parse_sequences(
            "; generated file\n>gene one\nATG GCT\nTAA\n>gene-two\nUUU\n",
            "sample.fa",
        )
        self.assertEqual(
            records,
            [
                SequenceRecord("gene one", "ATGGCTTAA"),
                SequenceRecord("gene-two", "UUU"),
            ],
        )

    def test_plain_sequence_whitespace_is_removed(self) -> None:
        self.assertEqual(
            parse_sequences("atg gct\n", "sample.seq"),
            [SequenceRecord("sample.seq", "atggct")],
        )

    def test_rejects_sequence_before_fasta_header(self) -> None:
        with self.assertRaisesRegex(ValueError, "before a FASTA header"):
            parse_sequences("ATG\n>gene\nGCT\n", "sample.fa")

    def test_rejects_empty_fasta_record(self) -> None:
        with self.assertRaisesRegex(ValueError, "record 'gene' is empty"):
            parse_sequences(">gene\n>other\nATG\n", "sample.fa")


class AnalysisTests(unittest.TestCase):
    def test_calculates_codon_usage_pair_bias_and_rscu(self) -> None:
        report = analyze_records([SequenceRecord("repeat", "ATGATGATG")])
        analysis = report.aggregate

        self.assertEqual(analysis.codon_counts, {"ATG": 3})
        self.assertEqual(analysis.pair_count, 2)
        metric = analysis.pair_metrics[0]
        self.assertEqual((metric.codon_1, metric.codon_2), ("ATG", "ATG"))
        self.assertEqual(metric.observed_count, 2)
        self.assertAlmostEqual(metric.expected_count, 2)
        self.assertAlmostEqual(metric.log2_observed_expected, 0)
        self.assertEqual(analysis.rscu["ATG"], 1)

    def test_aggregate_does_not_create_pairs_across_records(self) -> None:
        report = analyze_records(
            [SequenceRecord("one", "ATG"), SequenceRecord("two", "GCT")]
        )
        self.assertEqual(report.aggregate.pair_count, 0)
        self.assertEqual(len(report.aggregate.pair_metrics), 0)

    def test_rna_is_normalized_to_dna(self) -> None:
        report = analyze_records([SequenceRecord("rna", "AUG")])
        self.assertEqual(report.aggregate.codon_counts, {"ATG": 1})

    def test_frame_and_partial_codon_trim_are_applied(self) -> None:
        report = analyze_records(
            [SequenceRecord("framed", "AATGAA")],
            frame=1,
            partial_codon="trim",
        )
        self.assertEqual(report.aggregate.codon_counts, {"ATG": 1})
        self.assertEqual(report.aggregate.trimmed_bases, 2)

    def test_rejects_invalid_bases(self) -> None:
        with self.assertRaisesRegex(ValueError, "invalid nucleotide"):
            analyze_records([SequenceRecord("bad", "ATN")])

    def test_rejects_incomplete_codon_by_default(self) -> None:
        with self.assertRaisesRegex(ValueError, "sequence length must be divisible by 3"):
            analyze_records([SequenceRecord("partial", "ATGA")])

    def test_zero_pairs_include_unobserved_pairs_without_inventing_scores(self) -> None:
        report = analyze_records(
            [SequenceRecord("repeated", "ATGATGATG")],
            include_zero_pairs=True,
        )
        self.assertEqual(len(report.aggregate.pair_metrics), 4096)
        metric = next(
            metric for metric in report.aggregate.pair_metrics
            if (metric.codon_1, metric.codon_2) == ("ATG", "ATG")
        )
        self.assertEqual(metric.log2_observed_expected, 0)
        impossible = next(
            metric for metric in report.aggregate.pair_metrics
            if (metric.codon_1, metric.codon_2) == ("AAA", "AAA")
        )
        self.assertIsNone(impossible.log2_observed_expected)


if __name__ == "__main__":
    unittest.main()
