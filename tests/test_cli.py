import contextlib
import io
import json
import unittest
from unittest.mock import patch

from cpb_analyzer.cli import main


class CliTests(unittest.TestCase):
    def test_json_output_contains_methodology_and_analysis(self) -> None:
        output = io.StringIO()
        with patch("sys.stdin", io.StringIO(">gene\nATGATG\n")), contextlib.redirect_stdout(output):
            result = main(["-", "--format", "json"])

        report = json.loads(output.getvalue())
        self.assertEqual(result, 0)
        self.assertEqual(report["aggregate"]["pair_count"], 1)
        self.assertIn("not statistical significance tests", report["methodology"])

    def test_csv_output_is_pair_level(self) -> None:
        output = io.StringIO()
        with patch("sys.stdin", io.StringIO("ATGATG")), contextlib.redirect_stdout(output):
            result = main(["-", "--format", "csv"])

        self.assertEqual(result, 0)
        self.assertIn("observed_expected_ratio", output.getvalue().splitlines()[0])
        self.assertIn("ATG,ATG,1", output.getvalue())

    def test_text_does_not_label_positive_scores_as_depleted(self) -> None:
        output = io.StringIO()
        with patch("sys.stdin", io.StringIO("ATG GCT TTA CGA\n")), contextlib.redirect_stdout(output):
            result = main(["-", "--top-pairs", "1"])

        self.assertEqual(result, 0)
        self.assertNotIn("Most depleted observed pairs:", output.getvalue())

    def test_rejects_negative_top_pairs(self) -> None:
        with self.assertRaises(SystemExit) as error:
            main(["-", "--top-pairs", "-1"])
        self.assertEqual(error.exception.code, 2)


if __name__ == "__main__":
    unittest.main()
