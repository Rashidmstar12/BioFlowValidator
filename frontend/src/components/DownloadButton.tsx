import type { ValidationReport } from '../types/report';
import { getHtmlReportUrl } from '../api/client';

interface Props {
  jobId: string;
  report: ValidationReport;
}

export default function DownloadButton({ jobId, report }: Props) {
  function downloadJson() {
    const blob = new Blob([JSON.stringify(report, null, 2)], { type: 'application/json' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = `bioflowvalidator-report-${jobId}.json`;
    a.click();
    URL.revokeObjectURL(url);
  }

  return (
    <div className="flex gap-2 flex-wrap">
      <button
        onClick={downloadJson}
        className="px-4 py-2 bg-slate-700 hover:bg-slate-800 text-white text-sm font-medium rounded-lg transition-colors"
      >
        ⬇ Download JSON
      </button>
      <a
        href={getHtmlReportUrl(jobId)}
        target="_blank"
        rel="noopener noreferrer"
        className="px-4 py-2 bg-blue-600 hover:bg-blue-700 text-white text-sm font-medium rounded-lg transition-colors inline-block"
      >
        🌐 Open HTML Report
      </a>
    </div>
  );
}
