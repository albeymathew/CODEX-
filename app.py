"""CPB Analyzer — single-file Flask app (no subpackages needed)."""

import csv
import io as builtin_io
import os
from collections import Counter, defaultdict
from dataclasses import dataclass
from math import log2

from flask import Flask, make_response, render_template_string, request

# ── Genetic code ──────────────────────────────────────────────────────────────

GENETIC_CODE: dict[str, str] = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*", "TGA": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

ALL_CODONS: list[str] = sorted(GENETIC_CODE.keys())

_aa_map: dict[str, list[str]] = defaultdict(list)
for _c, _a in GENETIC_CODE.items():
    _aa_map[_a].append(_c)
AMINO_ACID_TO_CODONS: dict[str, list[str]] = {
    aa: sorted(codons) for aa, codons in sorted(_aa_map.items())
}

# ── Data classes ──────────────────────────────────────────────────────────────

@dataclass(frozen=True)
class SequenceRecord:
    name: str
    sequence: str

@dataclass(frozen=True)
class PairMetric:
    codon_1: str
    codon_2: str
    observed_count: int
    expected_count: float
    observed_frequency: float
    expected_probability: float
    observed_expected_ratio: float | None
    log2_observed_expected: float | None

@dataclass(frozen=True)
class SequenceAnalysis:
    name: str
    sequence_length: int
    frame: int
    trimmed_bases: int
    codons: tuple[str, ...] | None
    codon_count: int
    codon_counts: dict[str, int]
    codon_frequencies: dict[str, float]
    rscu: dict[str, float | None]
    pair_count: int
    pair_metrics: tuple[PairMetric, ...]

@dataclass(frozen=True)
class AnalysisReport:
    frame: int
    partial_codon_policy: str
    include_zero_pairs: bool
    records: tuple[SequenceAnalysis, ...]
    aggregate: SequenceAnalysis

# ── FASTA parser ──────────────────────────────────────────────────────────────

def parse_fasta(text: str) -> list[SequenceRecord]:
    records: list[SequenceRecord] = []
    name: str | None = None
    parts: list[str] = []
    has_header = any(l.lstrip().startswith(">") for l in text.splitlines())
    if not has_header:
        seq = "".join(text.split())
        if seq:
            records.append(SequenceRecord("sequence_1", seq))
        return records
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name is not None:
                records.append(SequenceRecord(name, "".join(parts)))
            name = line[1:].split()[0] if line[1:].strip() else "sequence"
            parts = []
        else:
            parts.append(line)
    if name is not None:
        records.append(SequenceRecord(name, "".join(parts)))
    return records

# ── Analysis ──────────────────────────────────────────────────────────────────

VALID_BASES = frozenset("ATGC")

def _normalize(record: SequenceRecord, frame: int, partial_codon: str):
    seq = "".join(record.sequence.split()).upper().replace("U", "T")
    invalid = sorted(set(seq) - VALID_BASES)
    if invalid:
        raise ValueError(f"{record.name}: invalid nucleotide(s) {', '.join(invalid)}")
    seq = seq[frame:]
    remainder = len(seq) % 3
    if remainder and partial_codon == "error":
        raise ValueError(f"{record.name}: sequence length not divisible by 3 after frame {frame}")
    trimmed = remainder
    if remainder:
        seq = seq[:-remainder]
    if not seq:
        raise ValueError(f"{record.name}: no complete codons remain")
    return seq, trimmed

def _from_counts(*, name, codons, codon_counts, pair_counts,
                 sequence_length, frame, trimmed_bases, include_zero_pairs):
    codon_total = sum(codon_counts.values())
    codon_frequencies = {c: n / codon_total for c, n in sorted(codon_counts.items())}
    rscu: dict[str, float | None] = {}
    for syn in AMINO_ACID_TO_CODONS.values():
        syn_total = sum(codon_counts.get(c, 0) for c in syn)
        for c in syn:
            rscu[c] = codon_counts.get(c, 0) * len(syn) / syn_total if syn_total else None

    pair_total = sum(pair_counts.values())
    pair_keys = (
        sorted(pair_counts) if not include_zero_pairs
        else [(a, b) for a in ALL_CODONS for b in ALL_CODONS]
    )
    metrics = []
    for first, second in pair_keys:
        obs = pair_counts.get((first, second), 0)
        exp_prob = (codon_counts.get(first, 0) / codon_total) * (codon_counts.get(second, 0) / codon_total)
        exp_count = pair_total * exp_prob
        if pair_total == 0 or exp_count == 0:
            ratio = score = None
        else:
            ratio = obs / exp_count
            score = log2(ratio) if ratio > 0 else None
        metrics.append(PairMetric(
            codon_1=first, codon_2=second,
            observed_count=obs, expected_count=exp_count,
            observed_frequency=obs / pair_total if pair_total else 0.0,
            expected_probability=exp_prob,
            observed_expected_ratio=ratio, log2_observed_expected=score,
        ))
    return SequenceAnalysis(
        name=name, sequence_length=sequence_length, frame=frame,
        trimmed_bases=trimmed_bases, codons=codons, codon_count=codon_total,
        codon_counts=dict(sorted(codon_counts.items())),
        codon_frequencies=codon_frequencies, rscu=dict(sorted(rscu.items())),
        pair_count=pair_total, pair_metrics=tuple(metrics),
    )

def analyze_records(records, *, frame=0, partial_codon="error", include_zero_pairs=False):
    if not records:
        raise ValueError("at least one sequence record is required")
    analyses, normalized = [], []
    for rec in records:
        seq, trimmed = _normalize(rec, frame, partial_codon)
        normalized.append((rec, seq, trimmed))
        codons = tuple(seq[i:i+3] for i in range(0, len(seq), 3))
        analyses.append(_from_counts(
            name=rec.name, codons=codons,
            codon_counts=Counter(codons),
            pair_counts=Counter(zip(codons, codons[1:])),
            sequence_length=len(seq), frame=frame,
            trimmed_bases=trimmed, include_zero_pairs=include_zero_pairs,
        ))
    agg_cc: Counter[str] = Counter()
    agg_pc: Counter[tuple[str, str]] = Counter()
    for a in analyses:
        agg_cc.update(a.codon_counts)
        c = a.codons or ()
        agg_pc.update(zip(c, c[1:]))
    aggregate = _from_counts(
        name="all records", codons=None, codon_counts=agg_cc, pair_counts=agg_pc,
        sequence_length=sum(a.sequence_length for a in analyses),
        frame=frame, trimmed_bases=sum(t for _, _, t in normalized),
        include_zero_pairs=include_zero_pairs,
    )
    return AnalysisReport(frame=frame, partial_codon_policy=partial_codon,
                          include_zero_pairs=include_zero_pairs,
                          records=tuple(analyses), aggregate=aggregate)

# ── Flask app ─────────────────────────────────────────────────────────────────

app = Flask(__name__)
app.secret_key = os.environ.get("SECRET_KEY", "cpb-analyzer-2024")

HTML = r"""
{% macro _pair_table(analysis, suffix) %}
<div class="table-responsive">
  <table class="table table-sm table-hover align-middle" id="tblPairs{{ suffix }}">
    <thead class="table-light"><tr>
      <th>Codon 1</th><th>AA 1</th><th>Codon 2</th><th>AA 2</th>
      <th>Observed</th><th>Expected</th><th>Obs Freq</th><th>Exp Prob</th>
      <th>O/E Ratio</th><th>log₂(O/E)</th>
    </tr></thead>
    <tbody>
    {% for m in analysis.pair_metrics %}
    <tr>
      <td><code>{{ m.codon_1 }}</code></td>
      <td><span class="badge-aa">{{ GENETIC_CODE.get(m.codon_1,'?') }}</span></td>
      <td><code>{{ m.codon_2 }}</code></td>
      <td><span class="badge-aa">{{ GENETIC_CODE.get(m.codon_2,'?') }}</span></td>
      <td>{{ m.observed_count }}</td>
      <td>{{ "%.3f"|format(m.expected_count) }}</td>
      <td>{{ "%.5f"|format(m.observed_frequency) }}</td>
      <td>{{ "%.5f"|format(m.expected_probability) }}</td>
      <td>{% if m.observed_expected_ratio is not none %}{{ "%.3f"|format(m.observed_expected_ratio) }}{% else %}&mdash;{% endif %}</td>
      <td>{% if m.log2_observed_expected is not none %}
        {% if m.log2_observed_expected > 0 %}<span class="log2-pos">+{{ "%.3f"|format(m.log2_observed_expected) }}</span>
        {% elif m.log2_observed_expected < 0 %}<span class="log2-neg">{{ "%.3f"|format(m.log2_observed_expected) }}</span>
        {% else %}<span class="log2-zero">0.000</span>{% endif %}
      {% else %}&mdash;{% endif %}</td>
    </tr>
    {% endfor %}
    </tbody>
  </table>
</div>
{% endmacro %}

{% macro _codon_table(analysis, suffix) %}
<div class="table-responsive">
  <table class="table table-sm table-hover align-middle" id="tblCodons{{ suffix }}">
    <thead class="table-light"><tr>
      <th>Codon</th><th>Amino Acid</th><th>Count</th><th>Frequency</th><th>RSCU</th>
    </tr></thead>
    <tbody>
    {% for codon, count in analysis.codon_counts.items() %}
    <tr>
      <td><code>{{ codon }}</code></td>
      <td><span class="badge-aa">{{ GENETIC_CODE.get(codon,'?') }}</span></td>
      <td>{{ count }}</td>
      <td>{{ "%.5f"|format(analysis.codon_frequencies.get(codon,0)) }}</td>
      <td>{% set r=analysis.rscu.get(codon) %}{% if r is not none %}{{ "%.3f"|format(r) }}{% else %}&mdash;{% endif %}</td>
    </tr>
    {% endfor %}
    </tbody>
  </table>
</div>
{% endmacro %}

<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8"/>
  <meta name="viewport" content="width=device-width,initial-scale=1"/>
  <title>CPB Analyzer — Codon Pair Bias</title>
  <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.3.3/dist/css/bootstrap.min.css" rel="stylesheet"/>
  <link href="https://cdn.datatables.net/1.13.8/css/dataTables.bootstrap5.min.css" rel="stylesheet"/>
  <link href="https://cdn.jsdelivr.net/npm/bootstrap-icons@1.11.3/font/bootstrap-icons.min.css" rel="stylesheet"/>
  <style>
    body{background:#f4f6fb}
    .hero{background:linear-gradient(135deg,#0d6efd,#6610f2);color:#fff;padding:2.5rem 0 2rem}
    .hero h1{font-size:2rem;font-weight:800}
    .card-panel{border:none;border-radius:1rem;box-shadow:0 2px 16px rgba(0,0,0,.08)}
    .stat-card{border-radius:.75rem;border:1px solid #dee2e6;background:#fff;padding:1.1rem 1.4rem}
    .stat-card .label{font-size:.78rem;color:#6c757d;text-transform:uppercase;letter-spacing:.05em}
    .stat-card .value{font-size:1.5rem;font-weight:700;color:#0d6efd}
    .badge-aa{font-family:monospace;font-size:.8rem;background:#e9ecef;color:#212529;border-radius:.3rem;padding:1px 5px}
    .upload-zone{border:2px dashed #adb5bd;border-radius:.75rem;padding:1.5rem;cursor:pointer;transition:border-color .2s}
    .upload-zone:hover{border-color:#0d6efd}
    textarea.form-control{font-family:monospace;font-size:.85rem}
    .log2-pos{color:#198754;font-weight:600}
    .log2-neg{color:#dc3545;font-weight:600}
    .log2-zero{color:#6c757d}
    #loadingOverlay{display:none;position:fixed;inset:0;background:rgba(255,255,255,.75);z-index:9999;align-items:center;justify-content:center;flex-direction:column}
    #loadingOverlay.active{display:flex}
  </style>
</head>
<body>
<div id="loadingOverlay">
  <div class="spinner-border text-primary mb-3" style="width:3rem;height:3rem;" role="status"></div>
  <p class="fw-semibold text-secondary">Analysing sequences…</p>
</div>
<nav class="navbar navbar-dark" style="background:#0d6efd">
  <div class="container">
    <span class="navbar-brand fw-bold"><i class="bi bi-dna me-2"></i>CPB Analyzer</span>
    <span class="text-white-50 small">Codon Pair Bias</span>
  </div>
</nav>
<div class="hero">
  <div class="container">
    <h1><i class="bi bi-graph-up-arrow me-2"></i>Codon Pair Bias Analyzer</h1>
    <p class="mb-0">Upload a FASTA file or paste nucleotide sequences to compute codon usage, RSCU, and codon-pair statistics.</p>
  </div>
</div>
<div class="container my-4">
  <div class="card card-panel mb-4">
    <div class="card-body p-4">
      <h5 class="fw-bold mb-3"><i class="bi bi-upload me-2 text-primary"></i>Input Sequence</h5>
      <form id="analyzeForm" method="POST" action="/analyze" enctype="multipart/form-data">
        <div class="mb-3">
          <label class="form-label fw-semibold">Upload FASTA file</label>
          <div class="upload-zone text-center" onclick="document.getElementById('fasta_file').click()">
            <i class="bi bi-file-earmark-text fs-2 text-secondary"></i>
            <p class="mb-0 text-secondary" id="fileLabel">Click to browse or drag & drop a <code>.fasta</code> / <code>.fa</code> / <code>.txt</code> file</p>
          </div>
          <input type="file" id="fasta_file" name="fasta_file" accept=".fasta,.fa,.txt,.seq" class="d-none" onchange="document.getElementById('fileLabel').textContent='📄 '+this.files[0].name"/>
        </div>
        <div class="text-center text-muted mb-3">— or —</div>
        <div class="mb-3">
          <label class="form-label fw-semibold">Paste sequence(s)</label>
          <textarea name="fasta_text" class="form-control" rows="5" placeholder=">sequence_1&#10;ATGAAACCC...">{{ fasta_text or '' }}</textarea>
        </div>
        <div class="row g-3 align-items-end mb-3">
          <div class="col-sm-4">
            <label class="form-label fw-semibold">Reading frame</label>
            <select name="frame" class="form-select">
              <option value="0" {% if frame==0 %}selected{% endif %}>Frame 0 (default)</option>
              <option value="1" {% if frame==1 %}selected{% endif %}>Frame 1</option>
              <option value="2" {% if frame==2 %}selected{% endif %}>Frame 2</option>
            </select>
          </div>
          <div class="col-sm-4">
            <label class="form-label fw-semibold">Partial codon</label>
            <select name="partial_codon" class="form-select">
              <option value="trim" {% if partial_codon=='trim' %}selected{% endif %}>Trim trailing bases</option>
              <option value="error" {% if partial_codon=='error' %}selected{% endif %}>Raise an error</option>
            </select>
          </div>
          <div class="col-sm-4 d-flex align-items-center pt-3">
            <div class="form-check">
              <input class="form-check-input" type="checkbox" name="include_zero_pairs" id="izp" {% if include_zero_pairs %}checked{% endif %}>
              <label class="form-check-label fw-semibold" for="izp">Include zero-count pairs</label>
            </div>
          </div>
        </div>
        <button type="submit" class="btn btn-primary px-4" onclick="document.getElementById('loadingOverlay').classList.add('active')">
          <i class="bi bi-play-circle me-1"></i>Analyse
        </button>
      </form>
    </div>
  </div>

  {% if error %}
  <div class="alert alert-danger d-flex align-items-center mb-4">
    <i class="bi bi-exclamation-triangle-fill me-2 fs-5"></i><div>{{ error }}</div>
  </div>
  {% endif %}

  {% if report %}
  {% set agg = report.aggregate %}
  <div class="row g-3 mb-4">
    <div class="col-6 col-md-3"><div class="stat-card"><div class="label">Sequences</div><div class="value">{{ report.records|length }}</div></div></div>
    <div class="col-6 col-md-3"><div class="stat-card"><div class="label">Total length (nt)</div><div class="value">{{ "{:,}".format(agg.sequence_length) }}</div></div></div>
    <div class="col-6 col-md-3"><div class="stat-card"><div class="label">Total codons</div><div class="value">{{ "{:,}".format(agg.codon_count) }}</div></div></div>
    <div class="col-6 col-md-3"><div class="stat-card"><div class="label">Codon pairs</div><div class="value">{{ "{:,}".format(agg.pair_count) }}</div></div></div>
  </div>

  {% if report.records|length > 1 %}
  <div class="card card-panel mb-3">
    <div class="card-body py-2 px-4 d-flex align-items-center gap-3 flex-wrap">
      <span class="fw-semibold text-secondary small">Viewing:</span>
      <div class="btn-group btn-group-sm">
        <input type="radio" class="btn-check" name="sv" id="svAgg" checked onchange="switchView('aggregate')">
        <label class="btn btn-outline-primary" for="svAgg">All (aggregate)</label>
        {% for rec in report.records %}
        <input type="radio" class="btn-check" name="sv" id="svR{{ loop.index }}" onchange="switchView('rec{{ loop.index }}')">
        <label class="btn btn-outline-secondary" for="svR{{ loop.index }}">{{ rec.name }}</label>
        {% endfor %}
      </div>
    </div>
  </div>
  {% endif %}

  <div class="card card-panel mb-5">
    <div class="card-body p-0">
      <ul class="nav nav-tabs px-3 pt-3">
        <li class="nav-item"><a class="nav-link active" data-bs-toggle="tab" href="#tabPairs"><i class="bi bi-grid-3x3 me-1"></i>Codon Pair Metrics</a></li>
        <li class="nav-item"><a class="nav-link" data-bs-toggle="tab" href="#tabCodons"><i class="bi bi-table me-1"></i>Codon Counts &amp; RSCU</a></li>
        {% if report.records|length > 1 %}<li class="nav-item"><a class="nav-link" data-bs-toggle="tab" href="#tabSeq"><i class="bi bi-list-ul me-1"></i>Per-Sequence</a></li>{% endif %}
      </ul>
      <div class="tab-content p-3">
        <div class="tab-pane fade show active" id="tabPairs">
          <div class="d-flex justify-content-between align-items-center mb-2 flex-wrap gap-2">
            <p class="mb-0 text-secondary small">log<sub>2</sub>(O/E) > 0 = over-represented | &lt; 0 = under-represented</p>
            <form method="POST" action="/download/codon-pairs" enctype="multipart/form-data" class="d-inline">
              <input type="hidden" name="fasta_text" value="{{ fasta_text or '' }}">
              <input type="hidden" name="frame" value="{{ frame }}">
              <input type="hidden" name="partial_codon" value="{{ partial_codon }}">
              {% if include_zero_pairs %}<input type="hidden" name="include_zero_pairs" value="on">{% endif %}
              <button class="btn btn-sm btn-outline-success"><i class="bi bi-download me-1"></i>Download CSV</button>
            </form>
          </div>
          <div id="pairsAggregate">{{ _pair_table(agg,'agg') }}</div>
          {% for rec in report.records %}<div id="pairsRec{{ loop.index }}" class="d-none">{{ _pair_table(rec,'rec'~loop.index) }}</div>{% endfor %}
        </div>
        <div class="tab-pane fade" id="tabCodons">
          <div id="codonsAggregate">{{ _codon_table(agg,'cagg') }}</div>
          {% for rec in report.records %}<div id="codonsRec{{ loop.index }}" class="d-none">{{ _codon_table(rec,'crec'~loop.index) }}</div>{% endfor %}
          <form method="POST" action="/download/codon-counts" enctype="multipart/form-data" class="mt-2">
            <input type="hidden" name="fasta_text" value="{{ fasta_text or '' }}">
            <input type="hidden" name="frame" value="{{ frame }}">
            <input type="hidden" name="partial_codon" value="{{ partial_codon }}">
            {% if include_zero_pairs %}<input type="hidden" name="include_zero_pairs" value="on">{% endif %}
            <button class="btn btn-sm btn-outline-success"><i class="bi bi-download me-1"></i>Download codon counts CSV</button>
          </form>
        </div>
        {% if report.records|length > 1 %}
        <div class="tab-pane fade" id="tabSeq">
          <table class="table table-sm table-hover" id="tblPerSeq">
            <thead class="table-light"><tr><th>Sequence</th><th>Length (nt)</th><th>Codons</th><th>Trimmed</th><th>Pairs</th></tr></thead>
            <tbody>
            {% for rec in report.records %}
            <tr><td><code>{{ rec.name }}</code></td><td>{{ "{:,}".format(rec.sequence_length) }}</td><td>{{ "{:,}".format(rec.codon_count) }}</td><td>{{ rec.trimmed_bases }}</td><td>{{ "{:,}".format(rec.pair_count) }}</td></tr>
            {% endfor %}
            </tbody>
          </table>
        </div>
        {% endif %}
      </div>
    </div>
  </div>
  {% endif %}
</div>
<footer class="text-center py-3 border-top" style="font-size:.82rem;color:#6c757d">
  CPB Analyzer &mdash; Codon Pair Bias &nbsp;|&nbsp; Built with Flask &amp; Bootstrap
</footer>
<script src="https://cdn.jsdelivr.net/npm/bootstrap@5.3.3/dist/js/bootstrap.bundle.min.js"></script>
<script src="https://code.jquery.com/jquery-3.7.1.min.js"></script>
<script src="https://cdn.datatables.net/1.13.8/js/jquery.dataTables.min.js"></script>
<script src="https://cdn.datatables.net/1.13.8/js/dataTables.bootstrap5.min.js"></script>
<script>
function switchView(key){
  document.querySelectorAll('[id^="pairsAggregate"],[id^="pairsRec"],[id^="codonsAggregate"],[id^="codonsRec"]').forEach(el=>el.classList.add('d-none'));
  if(key==='aggregate'){
    document.getElementById('pairsAggregate')?.classList.remove('d-none');
    document.getElementById('codonsAggregate')?.classList.remove('d-none');
  } else {
    const i=key.replace('rec','');
    document.getElementById('pairsRec'+i)?.classList.remove('d-none');
    document.getElementById('codonsRec'+i)?.classList.remove('d-none');
  }
}
$(document).ready(function(){
  $('table[id^="tblPairs"],table[id^="tblCodons"],#tblPerSeq').each(function(){
    if(!$.fn.DataTable.isDataTable(this)){
      $(this).DataTable({pageLength:25,order:[],language:{search:'Filter:'},dom:'<"d-flex justify-content-between align-items-center mb-2"lf>rtip'});
    }
  });
});
</script>
</body>
</html>
"""

def _get_inputs():
    fasta_text = ""
    if "fasta_file" in request.files and request.files["fasta_file"].filename:
        fasta_text = request.files["fasta_file"].read().decode("utf-8", errors="replace")
    if not fasta_text:
        fasta_text = request.form.get("fasta_text", "").strip()
    frame = int(request.form.get("frame", 0))
    partial_codon = request.form.get("partial_codon", "trim")
    include_zero_pairs = request.form.get("include_zero_pairs") == "on"
    return fasta_text, frame, partial_codon, include_zero_pairs

def _run(fasta_text, frame, partial_codon, include_zero_pairs):
    if not fasta_text:
        return None, "Please upload a FASTA file or paste a nucleotide sequence."
    try:
        records = parse_fasta(fasta_text)
        if not records:
            return None, "No valid FASTA records found."
        return analyze_records(records, frame=frame, partial_codon=partial_codon,
                               include_zero_pairs=include_zero_pairs), None
    except ValueError as e:
        return None, str(e)

def _render(**kwargs):
    return render_template_string(HTML, GENETIC_CODE=GENETIC_CODE, **kwargs)

@app.route("/", methods=["GET"])
def index():
    return _render()

@app.route("/analyze", methods=["POST"])
def analyze():
    fasta_text, frame, partial_codon, include_zero_pairs = _get_inputs()
    report, error = _run(fasta_text, frame, partial_codon, include_zero_pairs)
    return _render(report=report, error=error, fasta_text=fasta_text,
                   frame=frame, partial_codon=partial_codon,
                   include_zero_pairs=include_zero_pairs)

@app.route("/download/codon-pairs", methods=["POST"])
def download_codon_pairs():
    fasta_text, frame, partial_codon, include_zero_pairs = _get_inputs()
    report, error = _run(fasta_text, frame, partial_codon, include_zero_pairs)
    if error:
        return error, 400
    out = builtin_io.StringIO()
    w = csv.writer(out)
    w.writerow(["sequence","codon_1","codon_2","observed_count","expected_count",
                "observed_frequency","expected_probability","observed_expected_ratio","log2_observed_expected"])
    for a in list(report.records) + [report.aggregate]:
        for m in a.pair_metrics:
            w.writerow([a.name, m.codon_1, m.codon_2, m.observed_count,
                        f"{m.expected_count:.4f}", f"{m.observed_frequency:.6f}",
                        f"{m.expected_probability:.6f}",
                        f"{m.observed_expected_ratio:.4f}" if m.observed_expected_ratio is not None else "",
                        f"{m.log2_observed_expected:.4f}" if m.log2_observed_expected is not None else ""])
    resp = make_response(out.getvalue())
    resp.headers["Content-Type"] = "text/csv"
    resp.headers["Content-Disposition"] = "attachment; filename=codon_pairs.csv"
    return resp

@app.route("/download/codon-counts", methods=["POST"])
def download_codon_counts():
    fasta_text, frame, partial_codon, include_zero_pairs = _get_inputs()
    report, error = _run(fasta_text, frame, partial_codon, include_zero_pairs)
    if error:
        return error, 400
    out = builtin_io.StringIO()
    w = csv.writer(out)
    w.writerow(["sequence","codon","amino_acid","count","frequency","rscu"])
    for a in list(report.records) + [report.aggregate]:
        for codon, count in a.codon_counts.items():
            w.writerow([a.name, codon, GENETIC_CODE.get(codon,"?"), count,
                        f"{a.codon_frequencies.get(codon,0):.6f}",
                        f"{a.rscu.get(codon):.4f}" if a.rscu.get(codon) is not None else ""])
    resp = make_response(out.getvalue())
    resp.headers["Content-Type"] = "text/csv"
    resp.headers["Content-Disposition"] = "attachment; filename=codon_counts.csv"
    return resp

if __name__ == "__main__":
    app.run()
