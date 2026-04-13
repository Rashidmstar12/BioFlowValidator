# BioFlowValidator API Specification

Base URL: `http://localhost:8000`

---

## `GET /health`

Health check endpoint.

**Response 200:**
```json
{ "status": "ok" }
```

---

## `POST /upload`

Upload count matrix and optional metadata files.

**Request:** `multipart/form-data`
- `count_matrix` (required): TSV/CSV/XLSX count matrix file
- `metadata` (optional): TSV/CSV sample metadata file

**Constraints:**
- Max file size: 50 MB per file
- Allowed extensions: `.tsv`, `.csv`, `.txt`, `.xlsx`

**Response 200:**
```json
{
  "job_id": "550e8400-e29b-41d4-a716-446655440000",
  "count_matrix_filename": "counts.tsv",
  "count_matrix_size": 204800,
  "metadata_filename": "metadata.tsv",
  "metadata_size": 1024
}
```

**Response 400:** Unsupported file extension.  
**Response 413:** File exceeds size limit.

---

## `POST /validate`

Run validation engine on previously uploaded files.

**Request:** `application/json`
```json
{
  "job_id": "550e8400-e29b-41d4-a716-446655440000",
  "count_matrix_filename": "counts.tsv",
  "metadata_filename": "metadata.tsv"
}
```

**Response 200:**
```json
{
  "job_id": "550e8400-e29b-41d4-a716-446655440000",
  "status": "complete",
  "summary": {
    "error_count": 0,
    "warning_count": 2,
    "pass_count": 20,
    "skip_count": 3,
    "total_rules": 25
  }
}
```

**Response 404:** Job ID not found. Upload files first.

---

## `GET /report/results/{job_id}`

Retrieve full validation report as JSON.

**Response 200:**
```json
{
  "job_id": "...",
  "timestamp": "2025-01-01T12:00:00+00:00",
  "files": [
    { "filename": "counts.tsv", "size_bytes": 204800, "sha256": "abc123..." }
  ],
  "summary": {
    "error_count": 0,
    "warning_count": 2,
    "pass_count": 20,
    "skip_count": 3,
    "total_rules": 25
  },
  "results": [
    {
      "rule_id": "FMT-001",
      "category": "format",
      "severity": "ERROR",
      "status": "PASS",
      "message": "File encoding detected as 'utf-8' (UTF-8 compatible).",
      "affected_items": [],
      "suggestion": "",
      "details": {}
    }
  ]
}
```

**Response 404:** No report found for given job ID.

---

## `GET /report/{job_id}`

Retrieve full validation report as rendered HTML page.

**Response 200:** `text/html` — standalone HTML report with all rule results.  
**Response 404:** No report found for given job ID.
