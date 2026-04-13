"""HTMLExporter — renders a ValidationReport as a standalone HTML file."""
from __future__ import annotations

from pathlib import Path

from jinja2 import Environment, FileSystemLoader, select_autoescape

from app.models.rule_result import ValidationReport

_TEMPLATE_DIR = Path(__file__).parent.parent / "templates"


def to_html(report: ValidationReport) -> str:
    env = Environment(
        loader=FileSystemLoader(str(_TEMPLATE_DIR)),
        autoescape=select_autoescape(["html"]),
    )
    template = env.get_template("report.html.j2")
    return template.render(report=report, data=report.to_dict())
