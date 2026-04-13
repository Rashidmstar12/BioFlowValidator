import React, { useRef, useState } from 'react';

interface Props {
  onFilesSelected: (countMatrix: File, metadata?: File) => void;
  disabled?: boolean;
}

const ALLOWED = new Set(['.tsv', '.csv', '.txt', '.xlsx']);

function extOf(f: File): string {
  return f.name.slice(f.name.lastIndexOf('.')).toLowerCase();
}

export default function FileUploadZone({ onFilesSelected, disabled }: Props) {
  const [dragging, setDragging] = useState(false);
  const [countFile, setCountFile] = useState<File | null>(null);
  const [metaFile, setMetaFile] = useState<File | null>(null);
  const [error, setError] = useState<string>('');

  const countRef = useRef<HTMLInputElement>(null);
  const metaRef = useRef<HTMLInputElement>(null);

  function validate(f: File): boolean {
    if (!ALLOWED.has(extOf(f))) {
      setError(`Unsupported file type: ${extOf(f)}. Allowed: ${[...ALLOWED].join(', ')}`);
      return false;
    }
    if (f.size > 50 * 1024 * 1024) {
      setError('File exceeds 50 MB limit.');
      return false;
    }
    setError('');
    return true;
  }

  function handleCountChange(e: React.ChangeEvent<HTMLInputElement>) {
    const f = e.target.files?.[0];
    if (f && validate(f)) setCountFile(f);
  }

  function handleMetaChange(e: React.ChangeEvent<HTMLInputElement>) {
    const f = e.target.files?.[0];
    if (f && validate(f)) setMetaFile(f);
  }

  function handleDrop(e: React.DragEvent) {
    e.preventDefault();
    setDragging(false);
    const files = Array.from(e.dataTransfer.files);
    files.forEach((f) => {
      if (!validate(f)) return;
      const ext = extOf(f);
      if (!countFile && ALLOWED.has(ext)) setCountFile(f);
      else if (!metaFile && ALLOWED.has(ext)) setMetaFile(f);
    });
  }

  function handleSubmit() {
    if (!countFile) { setError('Please select a count matrix file.'); return; }
    onFilesSelected(countFile, metaFile ?? undefined);
  }

  return (
    <div className="space-y-4">
      <div
        className={`border-2 border-dashed rounded-xl p-8 text-center transition-colors ${dragging ? 'border-blue-500 bg-blue-50' : 'border-slate-300 bg-white'}`}
        onDragOver={(e) => { e.preventDefault(); setDragging(true); }}
        onDragLeave={() => setDragging(false)}
        onDrop={handleDrop}
      >
        <div className="text-4xl mb-2">📂</div>
        <p className="text-slate-600 text-sm">Drag &amp; drop files here, or use the buttons below.</p>
        <p className="text-slate-400 text-xs mt-1">Supported: TSV, CSV, TXT, XLSX · Max 50 MB</p>
      </div>

      <div className="grid grid-cols-1 sm:grid-cols-2 gap-3">
        <div>
          <label className="block text-xs font-semibold text-slate-500 mb-1 uppercase tracking-wide">
            Count Matrix <span className="text-red-500">*</span>
          </label>
          <button
            onClick={() => countRef.current?.click()}
            disabled={disabled}
            className="w-full text-left px-3 py-2 border border-slate-300 rounded-lg text-sm hover:bg-slate-50 truncate disabled:opacity-50"
          >
            {countFile ? `✅ ${countFile.name}` : 'Select count matrix file…'}
          </button>
          <input ref={countRef} type="file" accept=".tsv,.csv,.txt,.xlsx" className="hidden" onChange={handleCountChange} />
        </div>

        <div>
          <label className="block text-xs font-semibold text-slate-500 mb-1 uppercase tracking-wide">
            Sample Metadata <span className="text-slate-400">(optional)</span>
          </label>
          <button
            onClick={() => metaRef.current?.click()}
            disabled={disabled}
            className="w-full text-left px-3 py-2 border border-slate-300 rounded-lg text-sm hover:bg-slate-50 truncate disabled:opacity-50"
          >
            {metaFile ? `✅ ${metaFile.name}` : 'Select metadata file…'}
          </button>
          <input ref={metaRef} type="file" accept=".tsv,.csv,.txt,.xlsx" className="hidden" onChange={handleMetaChange} />
        </div>
      </div>

      {error && (
        <p className="text-red-600 text-sm bg-red-50 border border-red-200 rounded px-3 py-2">{error}</p>
      )}

      <button
        onClick={handleSubmit}
        disabled={disabled || !countFile}
        className="w-full py-3 bg-blue-600 hover:bg-blue-700 text-white font-semibold rounded-xl transition-colors disabled:opacity-50 disabled:cursor-not-allowed"
      >
        {disabled ? 'Validating…' : '🔬 Validate Files'}
      </button>
    </div>
  );
}
