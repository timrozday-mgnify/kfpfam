// Read hit rules side by side at one k-mer set: the numbers behind the plan's table.
import { loadCommon, loadShape } from '../web/load.js';
import { evaluate } from '../web/policy.js';

const common = await loadCommon();
const row = (name, r) => console.log(name.padEnd(46), 'recall', r.rec.recall.toFixed(3), 'p10', r.rec.familyP10.toFixed(2), 'under50', r.rec.below50.toFixed(3), 'cand', r.cand.mean.toFixed(1));
for (const [shape, value] of [['k4', 9], ['k5', 9.5], ['k5', 10], ['k6', 10.5], ['k6', 12]]) {
  const d = { ...common, ...(await loadShape(shape, common)), meta: common.meta };
  for (const set of ['all', 'reads']) {
    console.log(`\n${shape} T${value}, ${set}, GA`);
    const base = { shape, method: 'lr', value, target: 'fixed', truth: 'ga', set, bar: 22, anchorN: 2, anchorWindow: 3, anchorSum: 0 };
    row('diagonal >= 22', evaluate({ ...base, rule: 'diag' }, d));
    row('diagonal >= 25', evaluate({ ...base, rule: 'diag', bar: 25 }, d));
    for (const w of [0, 1, 3, 6]) row(`anchors, window ±${w}`, evaluate({ ...base, rule: 'anchors', anchorWindow: w }, d));
    row('anchors ±3, chain >= 24', evaluate({ ...base, rule: 'anchors', anchorSum: 24 }, d));
    row('3 anchors ±3', evaluate({ ...base, rule: 'anchors', anchorN: 3 }, d));
    row('diagonal >= 25 or anchors ±3', evaluate({ ...base, rule: 'either', bar: 25 }, d));
    row('diagonal >= 25 or anchors ±3 chain >= 24', evaluate({ ...base, rule: 'either', bar: 25, anchorSum: 24 }, d));
    row('diagonal >= 22 and anchors ±3', evaluate({ ...base, rule: 'both' }, d));
  }
}
