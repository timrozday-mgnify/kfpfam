// Candidates and recall as the diagonal bar rises, against the two-seed rule: real reads, GA.
import { loadCommon, loadReal, loadShape, withSet } from '../web/load.js';
import { evaluate } from '../web/policy.js';

const common = await loadCommon();
for (const [shape, value] of [['k5', 9], ['k4', 7]]) {
  const data = await loadShape(shape, common);
  data.real = await loadReal(shape, common);
  const d = withSet(common, data, 'real');
  const base = { shape, value, method: 'lr', target: 'fixed', truth: 'ga', set: 'real' };
  const show = (name, q) => {
    const r = evaluate({ ...base, ...q }, d);
    console.log(shape, 'T', value, name.padEnd(16), 'recall', r.rec.recall.toFixed(3), 'cand', r.cand.mean.toFixed(1), 'decoy', r.decoyCand.toFixed(1));
  };
  for (const bar of [20, 22, 24, 26, 28, 30, 34]) show(`diag >= ${bar}`, { rule: 'diag', bar });
  if (shape === 'k5') for (const w of [0, 3, 6]) show(`two anchors +-${w}`, { rule: 'anchors', anchorN: 2, anchorWindow: w, anchorSum: 0, bar: 22 });
}
