import { useState } from 'react';
import type { ValidationReport } from '../types/report';
import RuleCard from './RuleCard';
import DownloadButton from './DownloadButton';

interface Props {
  report: ValidationReport;
}

type Tab = 'all' | 'errors' | 'warnings' | 'passed' | 'skipped';

const TABS: { id: Tab; label: string }[] = [
  { id: 'all', label: 'All' },
  { id: 'errors', label: 'Errors' },
  { id: 'warnings', label: 'Warnings' },
  { id: 'passed', label: 'Passed' },
  { id: 'skipped', label: 'Skipped' },
];

export default function ReportViewer({ report }: Props) {
  const [activeTab, setActiveTab] = useState<Tab>('all');

  const filtered = report.results.filter((r) => {
    if (activeTab === 'all') return true;
    if (activeTab === 'errors') return r.status === 'FAIL' && r.severity === 'ERROR';
    if (activeTab === 'warnings') return r.status === 'FAIL' && r.severity === 'WARNING';
    if (activeTab === 'passed') return r.status === 'PASS';
    if (activeTab === 'skipped') return r.status === 'SKIP';
    return true;
  });

  const categories = [...new Set(filtered.map((r) => r.category))];
  const { error_count, warning_count, pass_count, skip_count } = report.summary;

  const overallStatus =
    error_count > 0
      ? { label: '❌ Validation Failed', cls: 'bg-red-50 border-red-200 text-red-800' }
      : warning_count > 0
      ? { label: '⚠️ Passed with Warnings', cls: 'bg-amber-50 border-amber-200 text-amber-800' }
      : { label: '✅ All Checks Passed', cls: 'bg-green-50 border-green-200 text-green-800' };

  return (
    <div className="space-y-6">
      {/* Overall status banner */}
      <div className={`border rounded-xl px-5 py-4 font-semibold text-lg ${overallStatus.cls}`}>
        {overallStatus.label}
      </div>

      {/* Summary stats */}
      <div className="grid grid-cols-2 sm:grid-cols-4 gap-3">
        {[
          { label: 'Errors', count: error_count, cls: 'bg-red-50 text-red-700' },
          { label: 'Warnings', count: warning_count, cls: 'bg-amber-50 text-amber-700' },
          { label: 'Passed', count: pass_count, cls: 'bg-green-50 text-green-700' },
          { label: 'Skipped', count: skip_count, cls: 'bg-indigo-50 text-indigo-700' },
        ].map(({ label, count, cls }) => (
          <div key={label} className={`rounded-xl p-4 text-center ${cls}`}>
            <div className="text-3xl font-bold">{count}</div>
            <div className="text-xs font-semibold uppercase tracking-wide mt-1">{label}</div>
          </div>
        ))}
      </div>

      {/* File metadata */}
      <div className="bg-white border border-slate-200 rounded-xl overflow-hidden">
        <div className="px-4 py-3 bg-slate-50 border-b border-slate-200">
          <span className="text-sm font-semibold text-slate-600">Submitted Files</span>
        </div>
        <table className="w-full text-sm">
          <thead className="bg-slate-50 text-xs text-slate-500 uppercase">
            <tr>
              <th className="px-4 py-2 text-left">Filename</th>
              <th className="px-4 py-2 text-right">Size</th>
              <th className="px-4 py-2 text-left hidden sm:table-cell">SHA-256</th>
            </tr>
          </thead>
          <tbody>
            {report.files.map((f) => (
              <tr key={f.filename} className="border-t border-slate-100">
                <td className="px-4 py-2 font-medium">{f.filename}</td>
                <td className="px-4 py-2 text-right text-slate-500">{(f.size_bytes / 1024).toFixed(1)} KB</td>
                <td className="px-4 py-2 font-mono text-xs text-slate-400 hidden sm:table-cell">{f.sha256.slice(0, 16)}…</td>
              </tr>
            ))}
          </tbody>
        </table>
      </div>

      {/* Download buttons */}
      <div className="flex gap-3">
        <DownloadButton jobId={report.job_id} report={report} />
      </div>

      {/* Tabs */}
      <div className="flex gap-1 border-b border-slate-200">
        {TABS.map((t) => (
          <button
            key={t.id}
            onClick={() => setActiveTab(t.id)}
            className={`px-4 py-2 text-sm font-medium rounded-t transition-colors ${
              activeTab === t.id
                ? 'border-b-2 border-blue-600 text-blue-600'
                : 'text-slate-500 hover:text-slate-800'
            }`}
          >
            {t.label}
          </button>
        ))}
      </div>

      {/* Results by category */}
      {categories.length === 0 ? (
        <p className="text-slate-400 text-sm text-center py-8">No results in this category.</p>
      ) : (
        categories.map((cat) => (
          <div key={cat}>
            <h3 className="text-xs font-bold text-slate-400 uppercase tracking-widest mb-2">
              {cat}
            </h3>
            <div className="space-y-2">
              {filtered
                .filter((r) => r.category === cat)
                .map((r) => (
                  <RuleCard key={r.rule_id} result={r} />
                ))}
            </div>
          </div>
        ))
      )}
    </div>
  );
}
