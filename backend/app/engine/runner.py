"""RuleRunner — executes all registered rules against a ValidationContext."""
from __future__ import annotations

import traceback

from app.engine.registry import get_rules
from app.models.context import ValidationContext
from app.models.rule_result import RuleResult


GATE_ON_FAILURE = {
    "SMP-001": ["SMP-002", "SMP-004", "BIO-001", "BIO-002", "BIO-003", "BIO-007"],
    "NRM-001": ["NRM-002", "NRM-003", "NRM-004", "NRM-005", "NRM-006"],
}


def run_all(context: ValidationContext) -> list[RuleResult]:
    """Run every registered rule and return their results.

    If a rule raises an unexpected exception it is caught and reported as an
    ERROR result so the overall report is still produced.
    """
    rules = get_rules()
    results: list[RuleResult] = []
    blocked_rules: set[str] = set()
    
    for rule in rules:
        rule_id = getattr(rule, "rule_id", "UNKNOWN")
        category = getattr(rule, "category", "unknown")
        severity = getattr(rule, "severity", "WARNING")
        
        if rule_id in blocked_rules:
            # Find the blocking rule ID
            blocking_rule_id = None
            for gate_id, dependent_ids in GATE_ON_FAILURE.items():
                if rule_id in dependent_ids:
                    blocking_rule_id = gate_id
                    break
            
            result = RuleResult(
                rule_id=rule_id,
                category=category,
                severity=severity,
                status="SKIP",
                message=f"Skipped: {blocking_rule_id} detected a fatal input error. This rule cannot produce valid results on the current data."
            )
            results.append(result)
            continue

        try:
            result = rule.run(context)
        except Exception as exc:  # noqa: BLE001
            result = RuleResult(
                rule_id=rule_id,
                category=category,
                severity=severity,
                status="FAIL",
                message=f"Rule raised an unexpected exception: {exc}",
                details={"traceback": traceback.format_exc()},
            )
            
        if result.status == "FAIL" and result.severity == "ERROR" and rule_id in GATE_ON_FAILURE:
            for dep_id in GATE_ON_FAILURE[rule_id]:
                blocked_rules.add(dep_id)
            
        results.append(result)
        
    return results


