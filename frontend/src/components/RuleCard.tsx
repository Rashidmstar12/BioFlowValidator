import type { RuleResult } from '../types/report';

interface Props {
  result: RuleResult;
}

const SEVERITY_STYLES: Record<string, string> = {
  ERROR: 'bg-red-100 text-red-800',
  WARNING: 'bg-amber-100 text-amber-800',
  INFO: 'bg-sky-100 text-sky-800',
};

const STATUS_STYLES: Record<string, string> = {
  PASS: 'bg-green-100 text-green-800',
  FAIL: '',
  SKIP: 'bg-indigo-100 text-indigo-800',
};

export default function RuleCard({ result }: Props) {
  const badgeStyle =
    result.status === 'PASS'
      ? STATUS_STYLES.PASS
      : result.status === 'SKIP'
      ? STATUS_STYLES.SKIP
      : SEVERITY_STYLES[result.severity] ?? 'bg-slate-100 text-slate-700';

  const badgeLabel =
    result.status === 'PASS' ? 'PASS' : result.status === 'SKIP' ? 'SKIP' : result.severity;

  const isFail = result.status === 'FAIL';

  return (
    <details
      className="bg-white border border-slate-200 rounded-lg overflow-hidden group"
      open={isFail}
    >
      <summary className="flex items-center gap-3 px-4 py-3 cursor-pointer list-none hover:bg-slate-50">
        <span className="font-mono text-xs text-slate-400 min-w-[5rem]">{result.rule_id}</span>
        <span className={`text-xs font-bold px-2 py-0.5 rounded-full ${badgeStyle}`}>
          {badgeLabel}
        </span>
        <span className="flex-1 text-sm text-slate-700">{result.message}</span>
        <span className="text-slate-400 text-xs group-open:rotate-180 transition-transform">▼</span>
      </summary>

      <div className="px-4 pb-4 border-t border-slate-100 space-y-2 pt-3">
        {result.suggestion && (
          <div className="bg-green-50 border-l-4 border-green-400 px-3 py-2 text-sm text-green-800 rounded">
            💡 {result.suggestion}
          </div>
        )}
        {result.affected_items.length > 0 && (
          <div>
            <p className="text-xs font-semibold text-slate-500 mb-1">Affected items:</p>
            <ul className="list-disc list-inside space-y-0.5">
              {result.affected_items.map((item, i) => (
                <li key={i} className="font-mono text-xs text-slate-600">{item}</li>
              ))}
            </ul>
          </div>
        )}
        {result.details && Object.keys(result.details).length > 0 && (
          <div className="mt-2">
            <details className="text-xs">
              <summary className="cursor-pointer text-indigo-600 hover:text-indigo-800 font-semibold underline mb-1">
                Show Rule Details
              </summary>
              <pre className="bg-slate-50 p-2 rounded border border-slate-200 overflow-x-auto text-[10px] font-mono text-slate-700 whitespace-pre-wrap">
                {JSON.stringify(result.details, null, 2)}
              </pre>
            </details>
          </div>
        )}
      </div>
    </details>
  );
}
