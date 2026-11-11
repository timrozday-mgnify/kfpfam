// @ts-check
/* global Plotly */
// Whole-sample sketches: one k-mer set per sample, no read structure. Numbers from measure/sketch_score.py.
const el = (id) => /** @type {HTMLInputElement} */ (document.getElementById(id));
const LAYOUT = { margin: { t: 34, r: 16, b: 48, l: 56 }, font: { size: 12 }, legend: { orientation: 'h', y: -0.28 } };
const BLUES = ['#c3d4ec', '#9db9e0', '#7a9fd3', '#5a86c4', '#3f6eb2', '#2c5898', '#1f4479', '#14305a'];
const BOUND = { diagonal: '#c8553d', gapped: '#8a3a2a' };
// Per-read costs from the plan's Progress log, core-ms per read, on the same reads.
const PER_READ = [
  { name: 'hmmsearch, no filter', sub: 70, all: 70, note: 'measured at 50,000 reads; linear in reads' },
  { name: 'Per-read pseudo-alignment (k4 T7 index)', sub: 29, all: 29, note: 'measured at 50,000 reads, 24 to 29 by run; linear in reads' },
  { name: 'nail, whole search', sub: 7.4, all: 1.1, note: 'measured at both depths; fixed cost dominates when small' },
];
const HMMSEARCH_MS = 70;

let data;

const curve = (run, rule) => (rule === 'null' ? run.null : run.best);
const measure = (row, by) => (by === 'hits' ? row.hits[0] : by === 'recall5' ? row.recall5 : row.recall);
const label = (run) => `k${run.k}${run.scaled > 1 ? ` 1 in ${run.scaled}` : ''}`;

function atRecall(rows, by, target) {
  // Highest bar that still reaches the target.
  let found = null;
  for (const r of rows) if (measure(r, by) >= target) found = r;
  return found;
}

function selected(depth) {
  return data.runs.filter((r) => r.depth === depth).sort((a, b) => a.k - b.k || a.scaled - b.scaled);
}

function floorBits(run) {
  return Math.log2(run.residues * data.median_len);
}

function plotTrade(depth, rule, by) {
  const runs = selected(depth);
  const traces = runs.map((run, i) => ({
    x: curve(run, rule).map((r) => r.kept), y: curve(run, rule).map((r) => measure(r, by)),
    text: curve(run, rule).map((r) => `${label(run)}, bar ${r.bar}; ${r.kept} models; ${r.decoy} on reversed reads`),
    hoverinfo: 'text+y', mode: 'lines', name: label(run), line: { width: 2, color: BLUES[Math.min(i + 1, BLUES.length - 1)], dash: run.scaled > 1 ? 'dot' : 'solid' },
  }));
  if (depth === 'sub') {
    for (const [name, rows] of Object.entries(data.bound)) {
      traces.push({ x: rows.map((r) => r.kept), y: rows.map((r) => measure(r, by)), text: rows.map((r) => `read-level ${name}, bar ${r.bar}; ${r.kept} models`), hoverinfo: 'text+y', mode: 'lines', name: `Bound: ${name}`, line: { width: 2, color: BOUND[name] } });
    }
  }
  const ceiling = data.present[depth];
  traces.push({ x: [ceiling, ceiling], y: [0, 1], mode: 'lines', name: 'Models present', hoverinfo: 'name', line: { width: 1, color: '#999', dash: 'dash' } });
  Plotly.react('plot-sketch-trade', traces, { ...LAYOUT, title: { text: 'Recall against models kept', font: { size: 13 } }, xaxis: { type: 'log', title: 'Models kept, of 30,134' }, yaxis: { title: by === 'hits' ? 'Share of GA hits kept' : 'Share of present models kept', range: [0, 1.02] } }, { responsive: true });
}

function plotNoise(depth) {
  const runs = selected(depth).filter((r) => r.scaled === 1);
  const traces = runs.flatMap((run, i) => {
    const color = BLUES[Math.min(i + 1, BLUES.length - 1)];
    return [
      { x: run.best.map((r) => r.bar), y: run.best.map((r) => Math.max(r.kept, 1)), mode: 'lines', name: label(run), line: { width: 2, color } },
      { x: run.best.map((r) => r.bar), y: run.best.map((r) => Math.max(r.decoy, 1)), mode: 'lines', name: `${label(run)} reversed`, line: { width: 2, color, dash: 'dot' }, showlegend: false },
    ];
  });
  const shapes = [];
  const known = runs.find((r) => r.residues);
  if (known) {
    const f = floorBits(known);
    shapes.push({ type: 'line', x0: f, x1: f, yref: 'paper', y0: 0, y1: 1, line: { color: '#999', dash: 'dash', width: 1 } });
  }
  const annotations = known ? [{ x: floorBits(known), yref: 'paper', y: 1, text: 'log₂(residues × columns)', showarrow: false, xanchor: 'left', font: { size: 11, color: '#666' } }] : [];
  Plotly.react('plot-sketch-noise', traces, { ...LAYOUT, shapes, annotations, title: { text: 'Models with a k-mer at this score (dotted: reversed reads)', font: { size: 13 } }, xaxis: { title: 'Best k-mer score in the sample, bits' }, yaxis: { type: 'log', title: 'Models' } }, { responsive: true });
}

function fmtX(n) { return `${(data.models / n).toFixed(1)}×`; }

function drawTable(depth, rule, by) {
  const targets = [0.99, 0.95, 0.9, 0.8];
  const head = `<tr><th>Sketch</th><th>T</th>${targets.map((t) => `<th colspan="3">${t * 100}% recall: kept, fewer, reversed</th>`).join('')}</tr>`;
  const row = (name, rows, t) => `<tr><td>${name}</td><td>${t}</td>${targets.map((target) => {
    const r = atRecall(rows, by, target);
    return r ? `<td>${r.kept.toLocaleString()}</td><td>${fmtX(r.kept)}</td><td>${r.decoy.toLocaleString()}</td>` : '<td colspan="3">not reached</td>';
  }).join('')}</tr>`;
  let body = selected(depth).map((run) => row(label(run), curve(run, rule), run.T)).join('');
  if (depth === 'sub') for (const [name, rows] of Object.entries(data.bound)) body += row(`Bound: read-level ${name}`, rows, '');
  el('sketch-table').innerHTML = head + body;
}

function drawCost(depth, by) {
  const reads = data.reads[depth];
  const key = depth === 'sub' ? 'sub' : 'all';
  const rows = PER_READ.map((c) => `<tr><td class="text">${c.name}</td><td>${c[key].toFixed(1)}</td><td></td><td></td><td class="text">${c.note}</td></tr>`);
  for (const run of selected(depth)) {
    if (run.score_core_s == null) continue;
    const ms = (run.sample_core_s + run.score_core_s) * 1000 / reads;
    const r = atRecall(run.best, by, 0.95);
    const total = r ? ms + HMMSEARCH_MS * r.kept / data.models : null;
    rows.push(`<tr><td class="text">Sketch ${label(run)} T ${run.T}: ${(run.sample_core_s * 1e6 / reads).toFixed(1)} µs per read + ${run.score_core_s.toFixed(0)} core-s per sample</td><td>${ms.toFixed(2)}</td><td>${(29 / ms).toFixed(1)}×</td><td>${total == null ? 'not reached' : total.toFixed(1)}</td><td class="text">measured sketch; hmmsearch scaled by models kept</td></tr>`);
  }
  el('sketch-cost').innerHTML = `<tr><th>Method</th><th>Core-ms per read</th><th>Against per-read pseudo-alignment</th><th>Plus hmmsearch on kept models, at 95% recall</th><th>Basis</th></tr>${rows.join('')}`;
}

function draw() {
  const depth = el('sketch-depth').value;
  const rule = el('sketch-rule').value;
  const by = el('sketch-by').value;
  plotTrade(depth, rule, by);
  plotNoise(depth);
  drawTable(depth, rule, by);
  drawCost(depth, by);
  const present = data.present[depth];
  el('sketch-readouts').innerHTML = [
    ['Models present at GA', present.toLocaleString(), depth === 'sub' ? 'hmmsearch truth' : 'nail truth'],
    ['Most any sample-level filter can save', fmtX(present), 'perfect filter: 30,134 ÷ present'],
    ['Present models with a single GA hit', `${(100 * data.singletons[depth] / present).toFixed(0)}%`, 'one read must clear the noise'],
    ['GA hits', data.ga_hits[depth].toLocaleString(), `${data.reads[depth].toLocaleString()} reads`],
  ].map(([k, v, n]) => `<div><b>${v}</b>${k}<br><small>${n}</small></div>`).join('');
}

export async function setupSketch() {
  const res = await fetch('data/sketch.json');
  if (!res.ok) { el('sketch-readouts').textContent = 'data/sketch.json is missing; run measure/sketch_score.py.'; return; }
  data = await res.json();
  for (const id of ['sketch-depth', 'sketch-rule', 'sketch-by']) el(id).addEventListener('input', draw);
  draw();
}
