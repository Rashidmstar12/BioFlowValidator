"""FastAPI application entry point."""
from __future__ import annotations

from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware

from app.routers import upload, validate, report

app = FastAPI(
    title="BioFlowValidator API",
    description="Validates RNA-seq bioinformatics workflows and detects common errors.",
    version="0.1.0",
)

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

app.include_router(upload.router, prefix="/upload", tags=["upload"])
app.include_router(validate.router, prefix="/validate", tags=["validate"])
app.include_router(report.router, prefix="/report", tags=["report"])


@app.get("/health", tags=["health"])
async def health() -> dict:
    return {"status": "ok"}
