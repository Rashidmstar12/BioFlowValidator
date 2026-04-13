"""RuleRunner — executes all registered rules against a ValidationContext."""
from __future__ import annotations

import traceback

from app.engine.registry import get_rules
from app.models.context import ValidationContext
from app.models.rule_result import RuleResult


def run_all(context: ValidationContext) -> list[RuleResult]:
    """Run every registered rule and return their results.

    If a rule raises an unexpected exception it is caught and reported as an
    ERROR result so the overall report is still produced.
    """
    rules = get_rules()
    results: list[RuleResult] = []
    for rule in rules:
        try:
            result = rule.run(context)
        except Exception as exc:  # noqa: BLE001
            result = RuleResult(
                rule_id=getattr(rule, "rule_id", "UNKNOWN"),
                category=getattr(rule, "category", "unknown"),
                severity="ERROR",
                status="FAIL",
                message=f"Rule raised an unexpected exception: {exc}",
                details={"traceback": traceback.format_exc()},
            )
        results.append(result)
    return results
