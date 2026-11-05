// What a loose filter passes on the real reads, for measure/gapped_test.py. Usage: SHAPE T BAR OUT_DIR
import { writeFileSync } from 'node:fs';
import { loadCommon, loadReal, loadShape, withSet } from '../web/load.js';
import { evaluate } from '../web/policy.js';

const [shape, value, bar, out] = process.argv.slice(2);
const common = await loadCommon();
const data = await loadShape(shape, common);
data.real = await loadReal(shape, common);
const d = withSet(common, data, 'real');
const r = evaluate({ shape, value: +value, method: 'lr', target: 'fixed', truth: 'all', set: 'real', rule: 'diag', bar: +bar }, d);
const { pairs, models } = d;
const lines = ['g\tset\tmodel\tweight\tscore\tga\tdiag'];
for (let g = 0; g < pairs.n; g++) {
  if (pairs.set[g].startsWith('real')) lines.push([g, pairs.set[g], pairs.model[g], pairs.weight[g], pairs.score[g], models.ga[pairs.model[g]], r.truthEv.best[g]].join('\t'));
}
writeFileSync(`${out}/pairs.tsv`, `${lines.join('\n')}\n`);
for (const [name, set, ev] of [['bg', d.bg, r.bgEv], ['decoy', d.decoy, r.decoyEv]]) {
  const rows = ['read\tmodel\tdiag\tscore'];
  for (let g = 0; g < set.rows.groups; g++) {
    if (ev.best[g] >= +bar) rows.push([set.read[g], set.model[g], ev.best[g], set.score[g]].join('\t'));
  }
  writeFileSync(`${out}/${name}.tsv`, `${rows.join('\n')}\n`);
}
console.log('reads', d.meta.reads, 'truth reads', d.meta.truthReads, 'candidates', r.cand.mean.toFixed(1), 'decoys', r.decoyCand.toFixed(1));
