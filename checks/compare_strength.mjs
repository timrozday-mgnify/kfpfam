// Recall and candidates by match strength, for policies and nail.
import { loadCommon, loadReal, loadShape, withSet } from '../web/load.js';
import { binsFor, byStrength, evaluate } from '../web/policy.js';

const common = await loadCommon();
const policies = [
  ['k5 T9, d22', { shape: 'k5', value: 9, bar: 22 }],
  ['k5 T8, d22', { shape: 'k5', value: 8, bar: 22 }],
  ['k4 T7, d22', { shape: 'k4', value: 7, bar: 22 }],
  ['k4 T6, d18', { shape: 'k4', value: 6, bar: 18 }],
  ['k5 T6, d18', { shape: 'k5', value: 6, bar: 18 }],
  ['k5 T8, d18', { shape: 'k5', value: 8, bar: 18 }],
  ['k5 T8, d20', { shape: 'k5', value: 8, bar: 20 }],
];
const shapes = {};
const fmt = (a) => a.map((x) => (x == null ? '  – ' : x.toFixed(2).padStart(5))).join(' ');
const scheme = process.argv[2] || 'ga';
console.log(binsFor(scheme).map((b) => b.label).join(' | '));
for (const [set, key] of [['real', 'real'], ['reads', 'reads'], ['f50', 'f50'], ['f33', 'f33']]) {
  console.log(`\n== ${set}`);
  for (const t of common.byBin.tools) {
    const r = t.recall[key]?.[scheme];
    if (r) console.log(t.tool.padEnd(18), fmt(r));
  }
  for (const [name, p] of policies) {
    shapes[p.shape] ??= await loadShape(p.shape, common);
    if (set === 'real') shapes[p.shape].real ??= await loadReal(p.shape, common);
    const d = withSet(common, shapes[p.shape], set);
    const q = { method: 'lr', target: 'fixed', rule: 'diag', truth: 'all', set, ...p };
    const res = evaluate(q, d);
    const b = byStrength(res, q, d, scheme);
    console.log(name.padEnd(18), fmt(b.map((x) => (x.weight ? x.ok / x.weight : null))), `cand ${res.cand.mean.toFixed(1)} decoy ${res.decoyCand.toFixed(1)} GB ${(res.cost.entries * 4 / 1e9).toFixed(1)}`);
    if (set === 'real' || set === 'reads') console.log('  real/decoy cand  ', fmt(b.map((x) => x.real)), '\n                   ', fmt(b.map((x) => x.decoy)));
  }
}
for (const t of common.byBin.tools) {
  for (const [k, c] of Object.entries(t.candidates[scheme])) console.log(t.tool.padEnd(18), k.padEnd(10), fmt(c));
}
