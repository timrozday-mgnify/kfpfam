// Policy engine against the separate Python scan (docs/plan.md, Seeds and read scores). GA, as its
// E-value truth used a much smaller Z. Candidates are compared loosely: different samples of reads.
import assert from 'node:assert';
import { loadCommon, loadShape } from '../web/load.js';
import { evaluate } from '../web/policy.js';

const common = await loadCommon();
const check = async (shape, value, bar, want) => {
  const d = { ...common, ...(await loadShape(shape, common)), meta: common.meta };
  const r = evaluate({ shape, method: 'lr', value, target: 'fixed', rule: 'diag', bar, truth: 'ga', set: 'reads' }, d);
  console.log(shape, value, bar, 'recall', r.rec.recall.toFixed(3), 'candidates', r.cand.mean.toFixed(1), 'entries', r.cost.entries.toExponential(2), 'hits', r.hits.toExponential(2));
  assert(Math.abs(r.rec.recall - want.recall) < 0.05, 'recall');
  assert(Math.abs(Math.log2(r.cand.mean / want.candidates)) < 1, 'candidates');
  assert(Math.abs(Math.log10(r.cost.entries / want.entries)) < 0.15, 'entries');
};
await check('k5', 10, 22, { recall: 0.923, candidates: 8.3, entries: 4.7e8 });
await check('k5', 8, 22, { recall: 0.991, candidates: 13.5, entries: 3.2e9 });
await check('k4', 8, 22, { recall: 0.955, candidates: 10.5, entries: 1.14e8 });
await check('k6', 12, 22, { recall: 0.877, candidates: 6.5, entries: 1.9e9 });
console.log('ok');
