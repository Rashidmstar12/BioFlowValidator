/** TypeScript types matching the backend ValidationReport model. */

export type RuleSeverity = 'ERROR' | 'WARNING' | 'INFO';
export type RuleStatus = 'PASS' | 'FAIL' | 'SKIP';

export interface RuleResult {
  rule_id: string;
  category: string;
  severity: RuleSeverity;
  status: RuleStatus;
  message: string;
  affected_items: string[];
  suggestion: string;
  details: Record<string, unknown>;
}

export interface FileMeta {
  filename: string;
  size_bytes: number;
  sha256: string;
}

export interface ReportSummary {
  error_count: number;
  warning_count: number;
  pass_count: number;
  skip_count: number;
  total_rules: number;
}

export interface ValidationReport {
  job_id: string;
  timestamp: string;
  files: FileMeta[];
  summary: ReportSummary;
  results: RuleResult[];
}

export interface UploadResponse {
  job_id: string;
  count_matrix_filename: string;
  count_matrix_size: number;
  metadata_filename?: string;
  metadata_size?: number;
}

export interface ValidateResponse {
  job_id: string;
  status: string;
  summary: ReportSummary;
}
