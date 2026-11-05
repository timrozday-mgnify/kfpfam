// @ts-check
/* global Plotly */
// One family: the score that chooses its k-mers, along its columns, and where the anchors fall.
import { loadProfile } from './load.js';
import { selector } from './policy.js';

const A = 'ACDEFGHIKLMNPQRSTVWY';
const Q = [0.0787945, 0.01516, 0.0535222, 0.0668298, 0.0397062, 0.0695071, 0.0229198, 0.0590092, 0.0594422, 0.0963728, 0.0237718, 0.0414386, 0.0482904, 0.0395639, 0.0540978, 0.0683364, 0.0540687, 0.0673417, 0.0114135, 0.0304133];
const total = Q.reduce((a, b) => a + b, 0);
for (let i = 0; i < 20; i++) Q[i] /= total;
const el = (id) => /** @type {HTMLInputElement} */ (document.getElementById(id));
const STEP = 0.5; // bits per bin of a window's score distribution
const compact = (x) => (x >= 1e4 || (x > 0 && x < 0.01) ? x.toExponential(1).replace('e+', 'e') : x >= 100 ? x.toFixed(0) : x.toPrecision(3));

/**
 * Exact k-mer counts of each window by score, for the axis the method chooses by: log-odds, or
 * -log2 emission probability. Returns per window a Float64Array of counts by bin.
 */
function windowDistributions(lo, length, cols, byProbability) {
  const span = cols[cols.length - 1] + 1;
  const wins = Math.max(0, length - span + 1);
  const NB = byProbability ? 240 : 320;
  const zero = byProbability ? 0 : 120; // bin of score 0
  const out = [];
  for (let w = 0; w < wins; w++) {
    let cur = new Float64Array(NB);
    cur[zero] = 1;
    for (const o of cols) {
      const next = new Float64Array(NB);
      for (let a = 0; a < 20; a++) {
        const s = lo[(w + o) * 20 + a] / 4;
        const v = byProbability ? -Math.log2(Q[a] * 2 ** s) : s;
        const sh = Math.round(v / STEP);
        for (let b = 0; b < NB; b++) {
          const c = cur[b];
          if (c) next[Math.max(0, Math.min(NB - 1, b + sh))] += c;
        }
      }
      cur = next;
    }
    out.push(cur);
  }
  return { out, wins, NB, zero, span };
}

function columnInfo(lo, length) {
  const info = new Float64Array(length);
  for (let c = 0; c < length; c++) {
    for (let a = 0; a < 20; a++) {
      const s = lo[c * 20 + a] / 4;
      info[c] += Q[a] * 2 ** s * s;
    }
  }
  return info;
}

/** k-mers of one window at or above a score, best first (at most `limit`). */
function topKmers(lo, w, cols, limit) {
  const best = cols.map((o) => Math.max(...Array.from({ length: 20 }, (_, a) => lo[(w + o) * 20 + a] / 4)));
  const rest = best.map((_, i) => best.slice(i + 1).reduce((a, b) => a + b, 0));
  const found = [];
  let floor = -Infinity;
  const walk = (i, score, logp, letters) => {
    if (score + (i < cols.length ? best[i] + rest[i] : 0) < floor) return;
    if (i === cols.length) {
      found.push({ kmer: letters, score, logp });
      if (found.length > limit * 3) {
        found.sort((x, y) => y.score - x.score);
        found.length = limit;
        floor = found[limit - 1].score;
      }
      return;
    }
    const order = Array.from({ length: 20 }, (_, a) => a).sort((x, y) => lo[(w + cols[i]) * 20 + y] - lo[(w + cols[i]) * 20 + x]);
    for (const a of order) {
      const s = lo[(w + cols[i]) * 20 + a] / 4;
      walk(i + 1, score + s, logp - Math.log2(Q[a] * 2 ** s), letters + A[a]);
    }
  };
  walk(0, 0, 0, '');
  found.sort((x, y) => y.score - x.score);
  return found.slice(0, limit);
}

let current = -1;
let token = 0;
let clickedWindow = null;
let lastArgs = null;

export function setupFamily(common, update) {
  const list = el('family-list');
  const { models, pairs } = common;
  const fams = [...new Set(pairs.model)];
  for (const m of fams.slice(0, 8000)) {
    const o = document.createElement('option');
    o.value = models.name[m];
    o.label = models.acc[m];
    list.append(o);
  }
  const counts = new Map();
  for (const m of pairs.model) counts.set(m, (counts.get(m) || 0) + 1);
  const withHits = [...counts].filter(([, c]) => c >= 3).map(([m]) => m);
  let top = withHits[0];
  for (const m of withHits) if (counts.get(m) > counts.get(top)) top = m;
  el('family').value = models.name[top];
  el('family-random').onclick = () => { el('family').value = models.name[withHits[Math.floor(Math.random() * withHits.length)]]; update(); };
  el('family').addEventListener('change', update);
  el('family-read').addEventListener('change', update);
  document.getElementById('plot-family').addEventListener('click', () => {});
}

function findModel(models, text) {
  const t = text.trim().toLowerCase();
  let i = models.name.findIndex((n) => n.toLowerCase() === t);
  if (i < 0) i = models.acc.findIndex((n) => n.toLowerCase() === t);
  return i;
}

export async function drawFamily(p, res, d, common) {
  const { models, pairs } = common;
  const m = findModel(models, el('family').value);
  if (m < 0) { el('family-note').textContent = 'No family of that name.'; return; }
  if (m !== current) clickedWindow = null;
  current = m;
  const mine = ++token;
  const lo = await loadProfile(m, models);
  if (mine !== token) return;
  const length = models.length[m];
  const cols = common.meta.shapes[p.shape];
  const byProb = p.method === 'topP' || p.method === 'sample';
  const { out: dist, wins, NB, zero, span } = windowDistributions(lo, length, cols, byProb);
  const info = columnInfo(lo, length);
  const windowInfo = Array.from({ length: wins }, (_, w) => cols.reduce((a, o) => a + info[w + o], 0));
  const order = windowInfo.map((v, w) => [v, w]).sort((a, b) => b[0] - a[0]);
  const keepCount = Math.ceil((p.windowKeep ?? 1) * wins);
  const kept = new Uint8Array(wins);
  for (let i = 0; i < keepCount; i++) kept[order[i][1]] = 1;
  const offset = res.offsets[m];
  const budget = Math.max(1, (p.value || 1) * 2 ** -offset);
  const score = (b) => (b - zero) * STEP;
  // Entries per window, and the score at which that window's list stops.
  const entries = new Float64Array(wins);
  const cut = new Array(wins).fill(null);
  for (let w = 0; w < wins; w++) {
    if (!kept[w]) continue;
    const c = dist[w];
    if (p.method === 'lr') {
      const T = res.T[m];
      for (let b = NB - 1; b >= 0 && score(b) >= T; b--) entries[w] += c[b];
      cut[w] = T;
    } else if (p.method === 'topLR') {
      let left = budget;
      for (let b = NB - 1; b >= 0 && left > 0; b--) { const t = Math.min(left, c[b]); entries[w] += t; left -= t; if (t > 0) cut[w] = score(b); }
    } else if (p.method === 'topP') {
      let left = budget;
      for (let b = 0; b < NB && left > 0; b++) { const t = Math.min(left, c[b]); entries[w] += t; left -= t; if (t > 0) cut[w] = score(b); }
    } else if (p.method === 'sample') {
      for (let b = 0; b < NB; b++) if (c[b]) entries[w] += c[b] * (1 - (1 - Math.min(1 - 1e-12, 2 ** -score(b))) ** budget);
    }
    entries[w] /= p.sparsity || 1;
  }
  // Anchors from real sequences of this family, as the policy keeps them.
  const anchorsAll = new Float64Array(wins);
  const anchorsKept = new Float64Array(wins);
  const select = selector(p, d.truth, pairs.model, res.T, res.offsets);
  let hits = 0, seeded = 0;
  for (let g = 0; g < pairs.n; g++) {
    if (pairs.model[g] !== m) continue;
    if (p.set !== 'all' && pairs.set[g] !== p.set) continue;
    if (p.set === 'all' && pairs.set[g] === 'reads') continue;
    hits++;
    let any = false;
    for (let r = d.truth.offs[g]; r < d.truth.offs[g + 1]; r++) {
      const c = d.truth.col[r];
      if (c >= wins) continue;
      if (d.truth.score[r] / 4 >= 6) anchorsAll[c]++;
      if (select(r, m)) { anchorsKept[c]++; any = true; }
    }
    if (any) seeded++;
  }
  // Clustering: blocks of windows holding anchors, gaps, the share of a read's positions served.
  const has = Array.from(entries, (e) => (p.method === 'words' ? 0 : e >= 1 ? 1 : 0));
  let blocks = 0, gap = 0, longestGap = 0;
  for (let w = 0; w < wins; w++) {
    if (has[w] && !has[w - 1]) blocks++;
    gap = has[w] ? 0 : gap + 1;
    longestGap = Math.max(longestGap, gap);
  }
  const sum = entries.reduce((a, b) => a + b, 0);
  const sorted = Array.from(entries).sort((a, b) => b - a);
  const top10 = sorted.slice(0, Math.max(1, Math.round(wins * 0.1))).reduce((a, b) => a + b, 0) / Math.max(sum, 1e-9);
  const L = +el('family-read').value;
  let one = 0, two = 0, starts = 0;
  for (let s = 0; s + L - span <= wins; s++) {
    starts++;
    let first = -1, pair = false, any = false;
    for (let w = s; w <= s + L - span && w < wins; w++) {
      if (!has[w]) continue;
      any = true;
      if (first < 0) first = w;
      else if (w - first >= span) { pair = true; break; }
    }
    if (any) one++;
    if (pair) two++;
  }
  if (!starts) { starts = 1; one = has.some(Boolean) ? 1 : 0; }
  const item = (title, value, sub) => `<div>${title}<b>${value}</b><small>${sub || ''}</small></div>`;
  el('family-readouts').innerHTML = [
    item(models.name[m], models.acc[m], `${length} columns; GA ${models.ga[m]}; ${models.nseq[m]} seed seqs${models.clan[m] ? `; ${models.clan[m]}` : ''}`),
    item('K-mers kept', compact(sum), `${compact(sum / Math.max(wins, 1))} per window`),
    item('Windows with an anchor', `${((100 * has.reduce((a, b) => a + b, 0)) / Math.max(wins, 1)).toFixed(0)}%`, `${blocks} blocks; longest gap ${longestGap}`),
    item('Top 10% of windows hold', `${(100 * top10).toFixed(0)}%`, 'of this family\'s k-mers'),
    item(`${L} aa reads with an anchor`, `${((100 * one) / starts).toFixed(0)}%`, 'of positions on the model'),
    item(`${L} aa reads with two`, `${((100 * two) / starts).toFixed(0)}%`, 'non-overlapping, ignoring indels'),
    item('Real hits seeded', hits ? `${seeded} of ${hits}` : 'no hits', 'in the chosen sequences'),
  ].join('');

  // Score distribution per window, as a heatmap; the policy's line through it.
  let lowest = NB, highest = 0;
  for (const c of dist) for (let b = 0; b < NB; b++) if (c[b] >= 0.5) { if (b < lowest) lowest = b; if (b > highest) highest = b; }
  const top = byProb ? Math.min(highest, lowest + 80) : highest;
  const bottom = byProb ? lowest : Math.max(lowest, zero); // log-odds below 0 are never chosen
  const ybins = [];
  for (let b = bottom; b <= top; b++) ybins.push(score(b));
  const z = ybins.map((_, i) => dist.map((c) => (c[bottom + i] > 0 ? Math.log10(c[bottom + i]) : null)));
  const xs = Array.from({ length: wins }, (_, w) => w + 1);
  const yLabel = byProb ? 'k-mer −log₂ P(emission)' : 'k-mer log-odds, bits';
  const traces = [
    { type: 'heatmap', x: xs, y: ybins, z, colorscale: 'Blues', reversescale: true, colorbar: { title: 'log₁₀ k-mers', len: 0.35, y: 0.82, thickness: 12 }, yaxis: 'y', hovertemplate: 'window %{x}<br>score %{y}<br>log₁₀ k-mers %{z:.1f}<extra></extra>' },
    { x: xs, y: cut, mode: 'lines', name: byProb ? 'Selection limit' : 'Selection threshold', line: { color: '#c8553d', width: 1.5 }, yaxis: 'y', connectgaps: false },
    { x: xs, y: Array.from(entries, (e) => (e > 0 ? e : null)), type: 'bar', name: 'K-mers kept, per window', marker: { color: '#3366aa' }, yaxis: 'y2' },
    { x: xs, y: anchorsAll, type: 'bar', name: 'Real seeds ≥ 6 bits', marker: { color: '#ccc' }, yaxis: 'y3' },
    { x: xs, y: anchorsKept, type: 'bar', name: 'Real seeds kept (anchors)', marker: { color: '#2a9d5c' }, yaxis: 'y3' },
    { x: Array.from(info, (_, c) => c + 1), y: Array.from(info), type: 'bar', name: 'Information per column', marker: { color: '#999' }, yaxis: 'y4' },
  ];
  Plotly.react('plot-family', traces, {
    margin: { t: 16, r: 16, b: 44, l: 64 }, font: { size: 12 }, barmode: 'overlay', hovermode: 'closest',
    legend: { orientation: 'h', y: -0.08 },
    xaxis: { title: 'Model column (window start)', range: [0, length + 1], anchor: 'y4' },
    yaxis: { domain: [0.62, 1], title: yLabel },
    yaxis2: { domain: [0.42, 0.59], title: 'K-mers kept', type: 'log' },
    yaxis3: { domain: [0.2, 0.39], title: 'Seeds in real hits' },
    yaxis4: { domain: [0, 0.17], title: 'Bits' },
  }, { responsive: true });
  const plot = /** @type {any} */ (document.getElementById('plot-family'));
  if (plot.removeAllListeners) plot.removeAllListeners('plotly_click');
  lastArgs = { lo, cols, wins, res, m, p, entries };
  plot.on('plotly_click', (e) => { clickedWindow = Math.round(e.points[0].x) - 1; showKmers(); });
  if (clickedWindow === null || clickedWindow >= wins) {
    let best = 0;
    for (let w = 1; w < wins; w++) if (windowInfo[w] > windowInfo[best]) best = w;
    clickedWindow = best;
  }
  showKmers();
  el('family-note').textContent = p.method === 'words'
    ? 'Seed-alignment words have no fixed window, so the kept-k-mer track is empty; the anchors track still shows where they land. Click a window for its top k-mers.'
    : `${byProb ? 'Selection here ranks by emission probability.' : 'Selection here ranks by log-odds: k-mers above the line are kept.'} Hash sparsity divides the kept counts, and an anchor track shows the same thinning in real sequences. Click a window for its top k-mers.`;
}

function showKmers() {
  const { lo, cols, wins, res, m, p, entries } = lastArgs;
  if (!wins) { el('family-kmers').innerHTML = ''; return; }
  const w = clickedWindow;
  const rows = topKmers(lo, w, cols, 12);
  const T = p.method === 'lr' ? res.T[m] : null;
  el('family-kmers').innerHTML = `<tr><th colspan="4">Window at column ${w + 1}: top k-mers, of ${compact(entries[w])} kept</th></tr><tr><th>K-mer</th><th>Log-odds</th><th>−log₂ P</th><th>Kept by T</th></tr>${rows.map((r) => `<tr><td><code>${r.kmer}</code></td><td>${r.score.toFixed(2)}</td><td>${r.logp.toFixed(1)}</td><td>${T === null ? '–' : r.score >= T ? 'yes' : 'no'}</td></tr>`).join('')}`;
}
