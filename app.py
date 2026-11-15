"""Flask web application for the CPB Analyzer."""

import csv
import io as builtin_io

from flask import Flask, make_response, render_template, request

from cpb.analysis import AnalysisReport, analyze_records
from cpb.genetic_code import GENETIC_CODE
from cpb.io import parse_fasta

app = Flask(__name__)

# Make GENETIC_CODE available in every template automatically
@app.context_processor
def inject_globals():
    return {"GENETIC_CODE": GENETIC_CODE}
app.secret_key = "cpb-analyzer-2024"


def _get_form_inputs():
    """Extract and validate common form inputs."""
    # Accept uploaded file or pasted text
    fasta_text = ""
    if "fasta_file" in request.files and request.files["fasta_file"].filename:
        fasta_text = request.files["fasta_file"].read().decode("utf-8", errors="replace")
    if not fasta_text:
        fasta_text = request.form.get("fasta_text", "").strip()

    frame = int(request.form.get("frame", 0))
    partial_codon = request.form.get("partial_codon", "trim")
    include_zero_pairs = request.form.get("include_zero_pairs") == "on"
    return fasta_text, frame, partial_codon, include_zero_pairs


def _run_analysis(fasta_text, frame, partial_codon, include_zero_pairs):
    """Parse FASTA and run analysis; returns (report, error_string)."""
    if not fasta_text:
        return None, "Please upload a FASTA file or paste a nucleotide sequence."
    try:
        records = parse_fasta(fasta_text)
        if not records:
            return None, "No valid FASTA records found in the input."
        report = analyze_records(
            records,
            frame=frame,
            partial_codon=partial_codon,
            include_zero_pairs=include_zero_pairs,
        )
        return report, None
    except ValueError as exc:
        return None, str(exc)


# ---------------------------------------------------------------------------
# Routes
# ---------------------------------------------------------------------------

@app.route("/", methods=["GET"])
def index():
    return render_template("index.html")


@app.route("/analyze", methods=["POST"])
def analyze():
    fasta_text, frame, partial_codon, include_zero_pairs = _get_form_inputs()
    report, error = _run_analysis(fasta_text, frame, partial_codon, include_zero_pairs)
    return render_template(
        "index.html",
        report=report,
        error=error,
        fasta_text=fasta_text,
        frame=frame,
        partial_codon=partial_codon,
        include_zero_pairs=include_zero_pairs,
    )


@app.route("/download/codon-pairs", methods=["POST"])
def download_codon_pairs():
    fasta_text, frame, partial_codon, include_zero_pairs = _get_form_inputs()
    report, error = _run_analysis(fasta_text, frame, partial_codon, include_zero_pairs)
    if error:
        return error, 400

    target = request.form.get("download_target", "aggregate")
    if target == "aggregate":
        analyses = [report.aggregate]
    else:
        analyses = list(report.records)

    output = builtin_io.StringIO()
    writer = csv.writer(output)
    writer.writerow([
        "sequence", "codon_1", "codon_2",
        "observed_count", "expected_count",
        "observed_frequency", "expected_probability",
        "observed_expected_ratio", "log2_observed_expected",
    ])
    for analysis in analyses:
        for m in analysis.pair_metrics:
            writer.writerow([
                analysis.name,
                m.codon_1, m.codon_2,
                m.observed_count,
                f"{m.expected_count:.4f}",
                f"{m.observed_frequency:.6f}",
                f"{m.expected_probability:.6f}",
                f"{m.observed_expected_ratio:.4f}" if m.observed_expected_ratio is not None else "",
                f"{m.log2_observed_expected:.4f}" if m.log2_observed_expected is not None else "",
            ])

    response = make_response(output.getvalue())
    response.headers["Content-Type"] = "text/csv"
    response.headers["Content-Disposition"] = "attachment; filename=codon_pairs.csv"
    return response


@app.route("/download/codon-counts", methods=["POST"])
def download_codon_counts():
    fasta_text, frame, partial_codon, include_zero_pairs = _get_form_inputs()
    report, error = _run_analysis(fasta_text, frame, partial_codon, include_zero_pairs)
    if error:
        return error, 400

    output = builtin_io.StringIO()
    writer = csv.writer(output)
    writer.writerow(["sequence", "codon", "amino_acid", "count", "frequency", "rscu"])

    from cpb.genetic_code import GENETIC_CODE
    for analysis in list(report.records) + [report.aggregate]:
        for codon, count in analysis.codon_counts.items():
            writer.writerow([
                analysis.name,
                codon,
                GENETIC_CODE.get(codon, "?"),
                count,
                f"{analysis.codon_frequencies.get(codon, 0):.6f}",
                f"{analysis.rscu.get(codon):.4f}" if analysis.rscu.get(codon) is not None else "",
            ])

    response = make_response(output.getvalue())
    response.headers["Content-Type"] = "text/csv"
    response.headers["Content-Disposition"] = "attachment; filename=codon_counts.csv"
    return response


if __name__ == "__main__":
    app.run(debug=True, host="0.0.0.0", port=8888)
