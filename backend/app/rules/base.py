"""Abstract base class for all validation rules."""
from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Literal

from app.models.context import ValidationContext
from app.models.rule_result import RuleResult, RuleSeverity, RuleStatus


class BaseRule(ABC):
    rule_id: str
    category: str
    severity: RuleSeverity
    description: str

    @abstractmethod
    def run(self, context: ValidationContext) -> RuleResult:
        ...

    def _pass(self, message: str = "Check passed.") -> RuleResult:
        return RuleResult(
            rule_id=self.rule_id,
            category=self.category,
            severity=self.severity,
            status="PASS",
            message=message,
        )

    def _fail(
        self,
        message: str,
        affected_items: list[str] | None = None,
        suggestion: str = "",
        details: dict | None = None,
    ) -> RuleResult:
        return RuleResult(
            rule_id=self.rule_id,
            category=self.category,
            severity=self.severity,
            status="FAIL",
            message=message,
            affected_items=affected_items or [],
            suggestion=suggestion,
            details=details or {},
        )

    def _skip(self, reason: str = "Prerequisites not met; rule skipped.") -> RuleResult:
        return RuleResult(
            rule_id=self.rule_id,
            category=self.category,
            severity=self.severity,
            status="SKIP",
            message=reason,
        )
