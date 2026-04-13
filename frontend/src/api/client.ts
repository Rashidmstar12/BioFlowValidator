/** API client wrappers for the BioFlowValidator backend. */
import type {
  UploadResponse,
  ValidateResponse,
  ValidationReport,
} from '../types/report';

const BASE = '';

export async function uploadFiles(
  countMatrix: File,
  metadata?: File
): Promise<UploadResponse> {
  const form = new FormData();
  form.append('count_matrix', countMatrix);
  if (metadata) form.append('metadata', metadata);

  const res = await fetch(`${BASE}/upload`, { method: 'POST', body: form });
  if (!res.ok) {
    const err = await res.json().catch(() => ({ detail: res.statusText }));
    throw new Error(err.detail ?? 'Upload failed');
  }
  return res.json();
}

export async function validateJob(jobId: string): Promise<ValidateResponse> {
  const res = await fetch(`${BASE}/validate`, {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify({ job_id: jobId }),
  });
  if (!res.ok) {
    const err = await res.json().catch(() => ({ detail: res.statusText }));
    throw new Error(err.detail ?? 'Validation failed');
  }
  return res.json();
}

export async function getResults(jobId: string): Promise<ValidationReport> {
  const res = await fetch(`${BASE}/report/results/${jobId}`);
  if (!res.ok) throw new Error('Report not found');
  return res.json();
}

export function getHtmlReportUrl(jobId: string): string {
  return `${BASE}/report/${jobId}`;
}
