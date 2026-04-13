export default function ValidationProgress() {
  return (
    <div className="flex flex-col items-center justify-center py-16 space-y-4">
      <div className="relative w-16 h-16">
        <div className="absolute inset-0 rounded-full border-4 border-blue-100"></div>
        <div className="absolute inset-0 rounded-full border-4 border-blue-600 border-t-transparent animate-spin"></div>
      </div>
      <p className="text-slate-600 font-medium">Running validation checks…</p>
      <p className="text-slate-400 text-sm">This usually takes a few seconds.</p>
    </div>
  );
}
