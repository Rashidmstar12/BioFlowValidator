import tempfile
from pathlib import Path
import pytest
from app import store
from app.models.rule_result import FileMeta, RuleResult, ValidationReport

def test_job_dir_and_report_persistence():
    # Setup temporary directory and files
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_path = Path(tmpdir)
        # Override upload directory for testing
        store._UPLOAD_DIR = tmp_path
        store._REGISTRY_FILE = tmp_path / "registry.json"
        
        # Clear caches
        store._job_dirs.clear()
        store._reports.clear()
        
        job_id = "test-job-uuid-1234"
        job_dir = tmp_path / job_id
        job_dir.mkdir()
        
        # 1. Register job directory
        store.register_job_dir(job_id, job_dir)
        
        # Clear memory cache to force reading from disk
        store._job_dirs.clear()
        
        # Retrieve it
        retrieved_dir = store.get_job_dir(job_id)
        assert retrieved_dir is not None
        assert retrieved_dir.resolve() == job_dir.resolve()
        
        # 2. Save validation report
        files = [FileMeta("counts.tsv", 500, "dummy-sha256")]
        results = [
            RuleResult(
                rule_id="FMT-001",
                category="format",
                severity="WARNING",
                status="PASS",
                message="File encoding is UTF-8"
            )
        ]
        report = ValidationReport.build(job_id, files, results)
        store.save(job_id, report)
        
        # Verify file is written to disk
        report_file = job_dir / "report.json"
        assert report_file.exists()
        
        # Clear memory cache
        store._reports.clear()
        
        # Retrieve report from disk
        retrieved_report = store.get(job_id)
        assert retrieved_report is not None
        assert retrieved_report.job_id == job_id
        assert retrieved_report.files[0].filename == "counts.tsv"
        assert retrieved_report.results[0].rule_id == "FMT-001"
        assert retrieved_report.results[0].status == "PASS"
