"""Report router — returns JSON or HTML reports for completed jobs."""
from __future__ import annotations

from fastapi import APIRouter, HTTPException
from fastapi.responses import HTMLResponse, JSONResponse

from app import store
from app.report.html_exporter import to_html
from app.report.json_exporter import to_json

router = APIRouter()


@router.get("/results/{job_id}")
async def get_results(job_id: str) -> JSONResponse:
    """Return the full validation report as JSON."""
    report = store.get(job_id)
    if report is None:
        raise HTTPException(status_code=404, detail=f"No report found for job '{job_id}'.")
    return JSONResponse(content=report.to_dict())


@router.get("/{job_id}", response_class=HTMLResponse)
async def get_html_report(job_id: str) -> HTMLResponse:
    """Return the validation report as a rendered HTML page."""
    report = store.get(job_id)
    if report is None:
        raise HTTPException(status_code=404, detail=f"No report found for job '{job_id}'.")
    html = to_html(report)
    return HTMLResponse(content=html)
