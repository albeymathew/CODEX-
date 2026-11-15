# 🧬 CPB Analyzer — Codon Pair Bias

A local web application for analysing **codon usage** and **codon-pair bias (CPB)** in nucleotide sequences. Upload a FASTA file (or paste a sequence) and instantly get codon counts, RSCU values, and codon-pair O/E statistics — all in a clean, sortable browser interface.

---

## ✨ Features

| Feature | Detail |
|---|---|
| **FASTA upload** | Drag & drop or browse `.fasta` / `.fa` / `.txt` files |
| **Paste sequences** | Bare sequence or multi-record FASTA |
| **Reading frame** | Frame 0, 1, or 2 |
| **Partial codon policy** | Trim trailing bases or raise an error |
| **Codon Pair Metrics** | Observed count, Expected count, O/E ratio, log₂(O/E) — colour-coded |
| **Codon Counts & RSCU** | Per-codon frequency and Relative Synonymous Codon Usage |
| **Per-sequence view** | Switch between individual records and the aggregate |
| **CSV export** | Download codon pair metrics or codon counts as CSV |

---

## 🚀 Quick Start

### 1. Clone the repo

```bash
git clone https://github.com/YOUR_USERNAME/cpb-analyzer.git
cd cpb-analyzer
```

### 2. Install dependencies

```bash
pip install -r requirements.txt
```

### 3. Run the app

```bash
python app.py
```

Then open **http://localhost:8888** in your browser.

> **Windows shortcut:** double-click `run.bat` — it installs deps and starts the server automatically.

---

## 📁 Project Structure

```
cpb-analyzer/
├── cpb/
│   ├── __init__.py
│   ├── genetic_code.py   # Standard genetic code (64 codons, AMINO_ACID_TO_CODONS)
│   ├── io.py             # FASTA parser → SequenceRecord
│   └── analysis.py       # Core analysis: codon counts, RSCU, pair metrics
├── templates/
│   └── index.html        # Jinja2 template (Bootstrap 5 + DataTables)
├── app.py                # Flask app + download routes
├── requirements.txt
├── run.bat               # Windows one-click launcher
└── README.md
```

---

## 📊 Output Explained

### Codon Pair Metrics

| Column | Meaning |
|---|---|
| `Observed` | Raw count of this codon pair in the sequence |
| `Expected` | Expected count assuming codon independence |
| `Obs Freq` | Observed pair frequency (observed / total pairs) |
| `Exp Prob` | Expected probability (freq(codon1) × freq(codon2)) |
| `O/E Ratio` | Observed / Expected — > 1 means over-represented |
| `log₂(O/E)` | Log-odds score — **green** = favoured, **red** = avoided |

### RSCU (Relative Synonymous Codon Usage)

$$\text{RSCU}_{ij} = \frac{X_{ij}}{\frac{1}{n_i} \sum_{j=1}^{n_i} X_{ij}}$$

Where $X_{ij}$ is the count of the $j$-th codon for amino acid $i$, and $n_i$ is the number of synonymous codons for that amino acid. RSCU = 1.0 means no bias; > 1.0 means preferred; < 1.0 means avoided.

---

## 🧪 Example

Paste this into the textarea and click **Analyse**:

```
>example_gene
ATGAAACCCGGGTTTATGAAACCCGGGTTT
```

---

## 🛠 Requirements

- Python 3.10+
- Flask 3.0+

---

## 📄 License

MIT License — see [LICENSE](LICENSE) for details.
