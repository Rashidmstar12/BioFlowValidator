import { useState } from 'react';
import FileUploadZone from '../components/FileUploadZone';
import ValidationProgress from '../components/ValidationProgress';
import ReportViewer from '../components/ReportViewer';
import { uploadFiles, validateJob, getResults } from '../api/client';
import type { ValidationReport } from '../types/report';

type AppState = 'idle' | 'loading' | 'done' | 'error';

export default function HomePage() {
  const [state, setState] = useState<AppState>('idle');
  const [report, setReport] = useState<ValidationReport | null>(null);
  const [errorMsg, setErrorMsg] = useState<string>('');

  async function handleFilesSelected(countMatrix: File, metadata?: File) {
    setState('loading');
    setErrorMsg('');
    try {
      const upload = await uploadFiles(countMatrix, metadata);
      await validateJob(upload.job_id);
      const result = await getResults(upload.job_id);
      setReport(result);
      setState('done');
    } catch (e) {
      setErrorMsg(e instanceof Error ? e.message : String(e));
      setState('error');
    }
  }

  function handleReset() {
    setState('idle');
    setReport(null);
    setErrorMsg('');
  }

  return (
    <div className="min-h-screen bg-gradient-to-br from-slate-50 to-blue-50">
      {/* Header */}
      <header className="bg-slate-900 text-white py-6 px-4 shadow-lg">
        <div className="max-w-4xl mx-auto">
          <h1 className="text-2xl font-bold tracking-tight">🧬 BioFlowValidator</h1>
          <p className="text-slate-400 text-sm mt-1">
            Validate RNA-seq differential expression workflows before analysis
          </p>
        </div>
      </header>

      <main className="max-w-4xl mx-auto px-4 py-8 space-y-8">
        {state === 'idle' && (
          <div className="bg-white rounded-2xl shadow-sm border border-slate-200 p-6">
            <h2 className="text-lg font-semibold text-slate-800 mb-1">Upload Your Files</h2>
            <p className="text-sm text-slate-500 mb-4">
              Upload your count matrix and optional sample metadata to detect common RNA-seq workflow errors.
            </p>
            <FileUploadZone onFilesSelected={handleFilesSelected} />
          </div>
        )}

        {state === 'loading' && <ValidationProgress />}

        {state === 'error' && (
          <div className="bg-red-50 border border-red-200 rounded-xl p-5 space-y-3">
            <p className="text-red-700 font-semibold">⚠️ Validation Error</p>
            <p className="text-red-600 text-sm">{errorMsg}</p>
            <button
              onClick={handleReset}
              className="px-4 py-2 bg-red-600 text-white text-sm rounded-lg hover:bg-red-700"
            >
              Try Again
            </button>
          </div>
        )}

        {state === 'done' && report && (
          <div className="space-y-4">
            <div className="flex items-center justify-between">
              <h2 className="text-lg font-semibold text-slate-800">Validation Report</h2>
              <button
                onClick={handleReset}
                className="text-sm text-blue-600 hover:text-blue-800 underline"
              >
                ← Validate another file
              </button>
            </div>
            <ReportViewer report={report} />
          </div>
        )}

        {/* About section */}
        {state === 'idle' && (
          <div className="bg-white rounded-2xl shadow-sm border border-slate-200 p-6">
            <h2 className="text-lg font-semibold text-slate-800 mb-3">What does BioFlowValidator check?</h2>
            <div className="grid grid-cols-1 sm:grid-cols-2 gap-3 text-sm text-slate-600">
              {[
                ['📋 Format', 'Encoding, delimiters, headers, non-numeric values, negative counts'],
                ['🔗 Sample Matching', 'Count matrix vs metadata sample ID agreement'],
                ['🧪 Gene IDs', 'Namespace consistency, duplicates, version suffixes'],
                ['📊 Normalization', 'Raw vs normalized counts, library size issues'],
                ['🔁 Replicates', 'Minimum biological replicates per condition'],
                ['🔬 Biology', 'Single-condition check, MT fraction, label sanity'],
              ].map(([title, desc]) => (
                <div key={title} className="flex gap-2">
                  <span className="text-base">{title.split(' ')[0]}</span>
                  <div>
                    <p className="font-medium text-slate-700">{title.split(' ').slice(1).join(' ')}</p>
                    <p className="text-slate-500 text-xs">{desc}</p>
                  </div>
                </div>
              ))}
            </div>
          </div>
        )}
      </main>

      <footer className="text-center text-xs text-slate-400 py-8">
        BioFlowValidator — transparent, rule-based bioinformatics workflow validation.
      </footer>
    </div>
  );
}
