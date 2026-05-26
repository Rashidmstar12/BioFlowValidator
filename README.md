# BioFlowValidator

> **A transparent, rule-based validator for RNA-seq differential expression analysis workflows.**

BioFlowValidator catches common scientific and computational errors in RNA-seq data *before* expensive analysis begins — acting as a pre-analysis guard rail for wet-lab biologists, students, and clinical researchers.

---

## Features

- ✅ **32 validation rules** across 5 categories (format, sample, gene ID, normalization, biology)
- 🔬 Detects: sample mismatches, mixed gene ID namespaces, pre-normalized counts, too few replicates, library size outliers, and more
- 📊 Human-readable HTML report + machine-readable JSON
- 🚀 REST API (FastAPI) + React/TypeScript frontend
- 🐳 Single-command Docker startup

---

## Quick Start

### Docker (recommended)

```bash
git clone https://github.com/Rashidmstar12/BioFlowValidator.git
cd BioFlowValidator
docker compose up --build
```

Open [http://localhost:3000](http://localhost:3000) in your browser.

### Local Development

**Backend:**
```bash
cd backend
pip install -r requirements.txt
uvicorn app.main:app --reload --port 8000
```

**Frontend:**
```bash
cd frontend
npm install
npm run dev
```

Open [http://localhost:5173](http://localhost:5173).

---

## Inputs

| File | Format | Required |
|---|---|---|
| Count matrix | TSV / CSV / XLSX (genes × samples or samples × genes) | ✅ |
| Sample metadata | TSV / CSV (sample IDs + condition column) | Optional |

---

## Validation Rule Categories

| Category | Rules | Description |
|---|---|---|
| **Format** | FMT-001 – FMT-008 | Encoding, delimiters, headers, duplicates, non-negatives, matrix orientation |
| **Sample** | SMP-001 – SMP-005 | Sample ID matching, duplicates, replicates, near-identical replicate diagnostics |
| **Gene ID** | GEN-001 – GEN-005 | Namespace consistency, duplicates, version suffixes, organism detection |
| **Normalization** | NRM-001 – NRM-006 | Integer counts, library size ratios, zero genes, duplicate count profiles |
| **Biology** | BIO-001 – BIO-008 | Single condition, MT fraction, label sanity, batch confounding, ERCC spike-ins |

See [`docs/validation_rules.md`](docs/validation_rules.md) for the full rule reference.

---

## Running Tests

```bash
cd backend
python -m pytest tests/ -v
```

Run the dataset benchmark:
```bash
python datasets/benchmark.py
```

---

## API Reference

See [`docs/api_spec.md`](docs/api_spec.md) or browse the interactive docs at `http://localhost:8000/docs`.

---

## Repository Structure

```
BioFlowValidator/
├── backend/           # Python FastAPI application
│   ├── app/
│   │   ├── engine/    # FileParser, RuleRegistry, RuleRunner
│   │   ├── models/    # RuleResult, ValidationReport, ValidationContext
│   │   ├── rules/     # format/, sample/, gene/, normalization/, biology/
│   │   ├── report/    # JSONExporter, HTMLExporter
│   │   └── routers/   # FastAPI route handlers
│   └── tests/         # Unit + integration tests
├── frontend/          # React + TypeScript + Vite SPA
├── datasets/          # Valid + faulty example datasets + benchmark
├── docs/              # API spec, validation rules reference
├── Dockerfile.backend
├── Dockerfile.frontend
└── docker-compose.yml
```

---

## Design Principles

- **Validation only** — no analysis, no statistical computation
- **Transparent** — every rule has a documented ID, description, and suggestion
- **Auditable** — JSON report includes file SHA-256 hash and timestamp
- **Scientifically conservative** — ambiguous cases produce WARNING not ERROR
- **Reproducible** — same inputs always produce identical outputs

---

## License

MIT
