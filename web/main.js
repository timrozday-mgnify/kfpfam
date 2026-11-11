// @ts-check
/* global Plotly */
import { barSweep, byStrength, estimateMs, evaluate, FEATURES, N_GRID, parseRows } from './policy.js';
import { loadCommon, loadReal, loadShape, withSet } from './load.js';
import { drawFamily, setupFamily } from './family.js';
import { setupSketch } from './sketch.js';
import { setupPromiscuous } from './promiscuous.js';

const el = (id) => /** @type {HTMLInputElement} */ (document.getElementById(id));
const compact = (x) => (x >= 1e4 || (x > 0 && x < 0.01) ? x.toExponential(1).replace('e+', 'e') : x >= 100 ? x.toFixed(0) : x.toPrecision(3));
const pct = (x) => `${(100 * x).toFixed(1)}%`;
const size = (gb) => (gb >= 1000 ? `${(gb / 1000).toPrecision(3)} TB` : gb >= 1 ? `${gb.toPrecision(3)} GB` : `${(gb * 1000).toPrecision(3)} MB`);
const LAYOUT = { margin: { t: 34, r: 16, b: 48, l: 56 }, font: { size: 12 }, legend: { orientation: 'h', y: -0.2 } };
const COLORS = { main: '#3366aa', second: '#c8553d', muted: '#999', accent: '#2a9d5c' };

let common;
const shapes = {};
const SHAPE_LABELS = { k4: '4-mer', k5: '5-mer', k6: '6-mer', k7: '7-mer', s5: 'spaced 5 of 7 (1101011)', s6: 'spaced 6 of 8 (11101011)', p33: 'pair 3+3, gap 3 (111000111)' };
const VALUE = {
  lr: { name: 'T, bits', min: 4, max: 24, step: 0.25, init: 9, fmt: (v) => v.toFixed(2) },
  topLR: { name: 'N per window', min: 0, max: 4, step: 0.05, init: 2, log: true },
  topP: { name: 'N per window', min: 0, max: 4, step: 0.05, init: 2, log: true },
  sample: { name: 'n draws per window', min: 0, max: 4, step: 0.05, init: 2, log: true },
  words: null,
};
const PRESETS = {
  'Plan default: 5-mers, T 9, bar 22': { shape: 'k5', method: 'lr', value: 9, bar: 22, rule: 'diag', target: 'fixed' },
  'Small index: 4-mers, T 7': { shape: 'k4', method: 'lr', value: 7, bar: 22, rule: 'diag', target: 'fixed' },
  'Same spend per column, 250': { shape: 'k5', method: 'lr', target: 'budget', budget: 250, bar: 22, rule: 'diag' },
  'Same emission coverage': { shape: 'k5', method: 'lr', target: 'coverage', coverage: 0.05, bar: 22, rule: 'diag' },
  'Sampled, 300 per window': { shape: 'k5', method: 'sample', value: 300, bar: 22, rule: 'diag' },
  'Words in seed alignment': { shape: 'k5', method: 'words', bar: 22, rule: 'diag' },
  'Two seeds, one diagonal': { shape: 'k5', method: 'lr', value: 10, rule: 'anchors', anchorN: 2, anchorWindow: 0, anchorSum: 0, target: 'fixed' },
  'Anchor chain, window ±3': { shape: 'k5', method: 'lr', value: 10, rule: 'anchors', anchorN: 2, anchorWindow: 3, anchorSum: 0, target: 'fixed' },
  'Diagonal, or anchors over an indel': { shape: 'k5', method: 'lr', value: 10, rule: 'either', bar: 24, anchorN: 2, anchorWindow: 3, anchorSum: 24, target: 'fixed' },
};
const CONTROL_IDS = ['shape', 'method', 'value', 'target', 'budget', 'coverage', 'sparsity', 'windowKeep', 'rule', 'bar', 'anchorN', 'anchorWindow', 'anchorSum', 'rankBits', 'clanAware', 'truth', 'set', 'c_lookup', 'c_hit', 'c_cand', 'c_bytes'];

function buildSlopes(holder, prefix) {
  const div = el(holder);
  for (const [f, spec] of Object.entries(FEATURES)) {
    const id = `${prefix}-${f}`;
    const l = document.createElement('label');
    l.innerHTML = `${spec.label} <span id="${id}-out"></span><input type="range" id="${id}" min="-4" max="4" step="0.25" value="0">`;
    div.append(l);
    CONTROL_IDS.push(id);
  }
}

function setValueRange(method) {
  const spec = VALUE[method];
  el('value-row').hidden = !spec;
  if (!spec) return;
  el('value-name').textContent = spec.name;
  const v = el('value');
  v.min = String(spec.min); v.max = String(spec.max); v.step = String(spec.step);
  v.value = String(spec.init);
}

function policy() {
  const method = el('method').value;
  const spec = VALUE[method];
  const raw = +el('value').value;
  const slopes = (prefix) => Object.fromEntries(Object.keys(FEATURES).map((f) => [f, +el(`${prefix}-${f}`).value]));
  return {
    shape: el('shape').value,
    method,
    value: spec?.log ? 10 ** raw : raw,
    target: method === 'lr' ? el('target').value : 'fixed',
    budget: 10 ** +el('budget').value,
    coverage: +el('coverage').value,
    sparsity: 2 ** +el('sparsity').value,
    windowKeep: +el('windowKeep').value,
    rule: el('rule').value,
    bar: +el('bar').value,
    anchorN: +el('anchorN').value,
    anchorWindow: +el('anchorWindow').value,
    anchorSum: +el('anchorSum').value,
    rankBits: el('rankBits').value === '' ? null : +el('rankBits').value,
    clanAware: el('clanAware').checked,
    truth: el('truth').value,
    set: el('set').value,
    seedSlopes: slopes('seed'),
    barSlopes: slopes('bar'),
  };
}
function costs() {
  return { lookup: +el('c_lookup').value, hit: +el('c_hit').value, candidate: +el('c_cand').value * 1000, bytes: +el('c_bytes').value };
}

function showInputs(p) {
  const spec = VALUE[p.method];
  if (spec) el('value-out').textContent = spec.log ? compact(p.value) : p.value.toFixed(2);
  el('target-row').hidden = p.method !== 'lr';
  el('budget-row').hidden = p.target !== 'budget';
  el('coverage-row').hidden = p.target !== 'coverage';
  el('value-row').hidden = !spec || p.target !== 'fixed';
  el('budget-out').textContent = compact(p.budget);
  el('coverage-out').textContent = pct(p.coverage);
  el('sparsity-out').textContent = String(p.sparsity);
  el('bar-row').hidden = !['diag', 'either', 'both'].includes(p.rule);
  el('anchor-rows').hidden = !['anchors', 'either', 'both'].includes(p.rule);
  el('anchorWindow-out').textContent = `±${p.anchorWindow}`;
  el('anchorSum-out').textContent = p.anchorSum ? `≥ ${p.anchorSum}` : 'any';
  el('bar-out').textContent = p.bar.toFixed(1);
  for (const id of ['c_lookup', 'c_hit', 'c_cand', 'c_bytes']) el(`${id}-out`).textContent = el(id).value;
  for (const pre of ['seed', 'bar']) {
    for (const f of Object.keys(FEATURES)) el(`${pre}-${f}-out`).textContent = (+el(`${pre}-${f}`).value).toFixed(2);
  }
}

const dataFor = (shape) => withSet(common, shapes[shape], el('set').value);
async function ensure(shape) {
  if (!shapes[shape]) shapes[shape] = await loadShape(shape, common);
  if (el('set').value === 'real' && !shapes[shape].real) shapes[shape].real = await loadReal(shape, common);
}

function summary(p, res) {
  const c = costs();
  const ms = estimateMs(p, res.hits, res.cand.mean, res.wRead, c);
  return {
    recall: res.rec.recall, familyMean: res.rec.familyMean, familyP10: res.rec.familyP10, below50: res.rec.below50,
    zero: res.rec.zero, families: res.rec.families.length, solidFamilies: res.rec.solidFamilies,
    candidates: res.cand.mean, measured: res.cand.measured, decoy: res.decoyCand,
    gb: (res.cost.entries * c.bytes) / 1e9, hits: res.hits, ms, speedup: common.hmmsearchMs / ms,
  };
}

function readouts(s) {
  const item = (title, value, sub, ok) => `<div class="${ok === undefined ? '' : ok ? 'good' : 'bad'}">${title}<b>${value}</b><small>${sub || ''}</small></div>`;
  el('readouts').innerHTML = [
    item('Read recall', pct(s.recall), 'abundance-weighted; gate ≥ 97%', s.recall >= 0.97),
    item('Family recall, p10', pct(s.familyP10), `${s.solidFamilies} families with 3+ hits; gate ≥ 85%`, s.familyP10 >= 0.85),
    item('Family recall, mean', pct(s.familyMean), `${s.families} families`),
    item('Families under 50%', pct(s.below50), `${pct(s.zero)} lose every hit`),
    item('Candidate models per read', compact(s.candidates), s.measured ? `of ${common.models.acc.length}; measured` : 'model, not measured'),
    item('Decoy candidates per read', s.decoy == null ? '–' : compact(s.decoy), s.decoy == null ? '' : `precision about ${pct(Math.max(0, 1 - s.decoy / Math.max(s.candidates, 1e-9)))}`),
    item('Index', size(s.gb), 'gate ≤ 16 GB', s.gb <= 16),
    item('Seed hits per read', compact(s.hits), 'from the background model'),
    item('Estimated time per read', `${compact(s.ms)} ms`, `hmmsearch ${common.hmmsearchMs} ms measured`),
    item('Estimated speed-up', `${compact(s.speedup)}×`, 'gate ≥ 20×; to be measured', s.speedup >= 20),
  ].join('');
}

function plotBars(p, res, d) {
  const bars = [];
  for (let b = 18; b <= 50; b += 1) bars.push(b);
  const sw = barSweep(res, p, d, bars);
  const key = p.truth === 'ga' ? 'recall_ga' : 'recall_e3';
  const traces = [
    { x: sw.map((r) => r.candidates), y: sw.map((r) => r.recall), text: sw.map((r) => `bar ${r.bar}`), mode: 'lines+markers', name: 'This k-mer set, read recall', line: { color: COLORS.main } },
    { x: sw.map((r) => r.candidates), y: sw.map((r) => r.familyP10), text: sw.map((r) => `bar ${r.bar}`), mode: 'lines', name: 'Family p10', line: { color: COLORS.second, dash: 'dot' } },
  ];
  const b = common.baselines.filter((x) => x[key] != null && x.candidates > 0);
  for (const tool of new Set(b.map((x) => x.tool))) {
    const pts = b.filter((x) => x.tool === tool);
    traces.push({ x: pts.map((x) => x.candidates), y: pts.map((x) => x[key]), text: pts.map((x) => `${tool} ${x.settings}, ${compact(x.ms_per_read)} ms/read`), mode: 'markers', name: tool, marker: { size: 9, symbol: tool.startsWith('nail') ? 'diamond' : 'square' } });
  }
  Plotly.react('plot-bars', traces, { ...LAYOUT, title: { text: 'Recall against candidates per read', font: { size: 13 } }, xaxis: { type: 'log', title: 'Candidate models per read' }, yaxis: { title: 'Recall', range: [0, 1.02] } }, { responsive: true });
}

let sweepToken = 0;
function plotSweep(p) {
  const token = ++sweepToken;
  const spec = VALUE[p.method];
  if (!spec || p.target !== 'fixed') {
    Plotly.react('plot-sweep', [], { ...LAYOUT, title: { text: 'Sweep needs a fixed value', font: { size: 13 } } }, { responsive: true });
    return;
  }
  const values = spec.log ? N_GRID : [4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16, 18, 20];
  const d = dataFor(p.shape);
  const rows = [];
  const step = (i) => {
    if (token !== sweepToken) return;
    if (i >= values.length) {
      Plotly.react('plot-sweep', [
        { x: rows.map((r) => r.gb), y: rows.map((r) => r.recall), text: rows.map((r) => `${compact(r.v)}; ${compact(r.candidates)} candidates; ${compact(r.hits)} hits`), mode: 'lines+markers', name: 'Read recall', line: { color: COLORS.main } },
        { x: rows.map((r) => r.gb), y: rows.map((r) => r.familyP10), text: rows.map((r) => compact(r.v)), mode: 'lines+markers', name: 'Family p10', line: { color: COLORS.second } },
      ], { ...LAYOUT, title: { text: 'Recall against index size', font: { size: 13 } }, xaxis: { type: 'log', title: 'Index, GB' }, yaxis: { title: 'Recall', range: [0, 1.02] } }, { responsive: true });
      return;
    }
    const q = { ...p, value: values[i] };
    const r = summary(q, evaluate(q, d));
    rows.push({ v: values[i], ...r });
    setTimeout(() => step(i + 1), 0);
  };
  step(0);
}

function plotStrength(p, res, d) {
  const scheme = el('scheme').value;
  const bins = byStrength(res, p, d, scheme);
  const labels = bins.map((b) => b.label);
  const traces = [{ x: labels, y: bins.map((b) => (b.weight ? b.ok / b.weight : null)), mode: 'lines+markers', name: 'This policy', line: { color: COLORS.main, width: 3 } }];
  for (const t of common.byBin?.tools || []) {
    const series = t.recall?.[p.set === 'all' ? 'fragments' : p.set]?.[scheme];
    // Relaxed nail settings start off the plot, one click away in the legend
    if (series) traces.push({ x: labels, y: series, mode: 'lines+markers', name: t.tool, line: { dash: 'dot' }, visible: /^nail (seeds|report)$|^MMseqs2/.test(t.tool) ? true : 'legendonly' });
  }
  Plotly.react('plot-strength', traces, { ...LAYOUT, title: { text: 'Recall by match strength', font: { size: 13 } }, yaxis: { title: 'Recall', range: [0, 1.02] }, xaxis: { type: 'category', tickangle: -30, automargin: true } }, { responsive: true });
  const rated = d.meta.truthReads;
  const prec = (b) => (b.real > 0 ? Math.max(0, 1 - b.decoy / b.real) : null);
  Plotly.react('plot-precision', [
    { x: labels, y: bins.map((b) => b.real), type: 'bar', name: 'Reads', marker: { color: COLORS.main } },
    { x: labels, y: bins.map((b) => b.decoy), type: 'bar', name: 'Reversed reads', marker: { color: COLORS.second } },
  ], { ...LAYOUT, barmode: 'group', title: { text: 'Candidates per read, by how hmmsearch rates them', font: { size: 13 } }, yaxis: { title: 'Candidates per read', type: 'log' }, xaxis: { type: 'category', tickangle: -30, automargin: true } }, { responsive: true });
  el('strength-table').innerHTML = `<tr><th>Bin</th><th>Hits in sample</th><th>Recall</th><th>Candidates / read</th><th>Decoy candidates / read</th><th>Precision (est.)</th>${rated ? '<th>hmmsearch hits / read</th><th>Decoy hits / read</th><th>hmmsearch precision (est.)</th>' : ''}</tr>${
    bins.map((b) => `<tr><td>${b.label}</td><td>${b.n}</td><td>${b.weight ? pct(b.ok / b.weight) : '–'}</td><td>${compact(b.real)}</td><td>${compact(b.decoy)}</td><td>${prec(b) === null ? '–' : pct(prec(b))}</td>${rated ? `<td>${compact(b.hits)}</td><td>${compact(b.decoyHits)}</td><td>${b.hits ? pct(Math.max(0, 1 - b.decoyHits / b.hits)) : '–'}</td>` : ''}</tr>`).join('')}`;
}

function familyRows(res) {
  const { models } = common;
  return res.rec.families.map((f) => ({ ...f, n_: f.n, entries: res.cost.perModelN[f.m], perColumn: res.cost.perModelN[f.m] / models.length[f.m] }));
}

function plotFamilies(res) {
  const { models } = common;
  const fam = familyRows(res);
  const f = el('feature').value;
  const spec = FEATURES[f];
  const xOf = (r) => models[f][r.m];
  const sorted = [...fam].sort((a, b) => xOf(a) - xOf(b));
  const nb = 10;
  const bins = Array.from({ length: nb }, (_, i) => sorted.slice(Math.floor((i * sorted.length) / nb), Math.floor(((i + 1) * sorted.length) / nb)));
  const labels = bins.map((b) => (b.length ? `${compact(xOf(b[0]))}–${compact(xOf(b[b.length - 1]))}` : ''));
  const mean = (b, fn) => b.reduce((a, r) => a + fn(r), 0) / Math.max(b.length, 1);
  const p10 = (b) => { const v = b.map((r) => r.recall).sort((x, y) => x - y); return v.length ? v[Math.floor(0.1 * v.length)] : null; };
  Plotly.react('plot-bins', [
    { x: labels, y: bins.map((b) => mean(b, (r) => r.recall)), type: 'bar', name: 'Mean family recall', marker: { color: COLORS.main } },
    { x: labels, y: bins.map(p10), mode: 'lines+markers', name: 'p10', line: { color: COLORS.second } },
  ], { ...LAYOUT, title: { text: `Recall by ${spec.label.toLowerCase()}`, font: { size: 13 } }, yaxis: { range: [0, 1.02], title: 'Recall' }, xaxis: { tickangle: -30, automargin: true } }, { responsive: true });
  Plotly.react('plot-hist', [{ x: fam.map((r) => r.recall), type: 'histogram', xbins: { start: 0, end: 1.0001, size: 0.05 }, marker: { color: COLORS.main } }], { ...LAYOUT, title: { text: 'Family recall', font: { size: 13 } }, xaxis: { title: 'Recall in the family', range: [0, 1.02] }, yaxis: { title: 'Families' }, showlegend: false }, { responsive: true });
  Plotly.react('plot-scatter', [{ x: fam.map(xOf), y: fam.map((r) => r.recall + (Math.random() - 0.5) * 0.02), text: fam.map((r) => `${models.name[r.m]} ${models.acc[r.m]}, ${r.n} hits`), mode: 'markers', marker: { size: 4, opacity: 0.4, color: COLORS.main } }], { ...LAYOUT, title: { text: 'Each family', font: { size: 13 } }, xaxis: { title: spec.label, type: spec.log ? 'log' : 'linear' }, yaxis: { title: 'Recall (jittered)', range: [-0.03, 1.05] }, showlegend: false }, { responsive: true });
  Plotly.react('plot-spend', [
    { x: labels, y: bins.map((b) => mean(b, (r) => r.perColumn)), type: 'bar', name: 'Entries per column', marker: { color: COLORS.muted } },
  ], { ...LAYOUT, title: { text: 'Index spent per column, same bins', font: { size: 13 } }, yaxis: { type: 'log', title: 'Entries per column' }, xaxis: { tickangle: -30, automargin: true }, showlegend: false }, { responsive: true });
  const worst = fam.filter((r) => r.n >= 3).sort((a, b) => a.recall - b.recall || b.weight - a.weight).slice(0, 25);
  el('worst').innerHTML = `<tr><th>Family</th><th>Hits</th><th>Recall</th><th>Length</th><th>GA</th><th>Seed seqs</th><th>Mean info</th><th>Best window</th><th>T</th><th>Entries/col</th></tr>${
    worst.map((r) => `<tr data-model="${r.m}" style="cursor:pointer"><td>${models.name[r.m]} <small>${models.acc[r.m]}</small></td><td>${r.n}</td><td>${r.recall.toFixed(2)}</td><td>${models.length[r.m]}</td><td>${models.ga[r.m]}</td><td>${models.nseq[r.m]}</td><td>${models.mean_info[r.m].toFixed(2)}</td><td>${models.top_info[r.m].toFixed(1)}</td><td>${res.T[r.m].toFixed(1)}</td><td>${compact(r.perColumn)}</td></tr>`).join('')}`;
}

const saved = (() => { try { return JSON.parse(localStorage.getItem('kfpfam-saved') || '[]'); } catch { return []; } })();
function persist() { try { localStorage.setItem('kfpfam-saved', JSON.stringify(saved)); } catch { /* optional convenience */ } }
function tableRows(rows, withLoad) {
  return `<tr><th>Policy</th><th>Recall</th><th>Family p10</th><th>Under 50%</th><th>Candidates</th><th>Index</th><th>Hits</th><th>ms/read</th><th>Speed-up</th>${withLoad ? '<th></th>' : ''}</tr>${
    rows.map((r, i) => `<tr><td>${r.name}</td><td>${pct(r.s.recall)}</td><td>${pct(r.s.familyP10)}</td><td>${pct(r.s.below50)}</td><td>${compact(r.s.candidates)}</td><td>${size(r.s.gb)}</td><td>${compact(r.s.hits)}</td><td>${compact(r.s.ms)}</td><td>${compact(r.s.speedup)}×</td>${withLoad ? `<td><button type="button" data-load="${i}">Load</button> <button type="button" data-del="${i}">Remove</button></td>` : ''}</tr>`).join('')}`;
}
function drawSaved() { el('saved').innerHTML = saved.length ? tableRows(saved, true) : ''; }

function apply(state) {
  for (const [k, v] of Object.entries(state)) {
    const c = el(k);
    if (!c) continue;
    if (k === 'method') setValueRange(String(v));
    if (c.type === 'checkbox') c.checked = !!v;
    else c.value = String(v);
  }
  if (state.method && VALUE[state.method]?.log && state.value) el('value').value = String(Math.log10(state.value));
  if (state.budget) el('budget').value = String(Math.log10(state.budget));
}
function formState() {
  return Object.fromEntries(CONTROL_IDS.map((id) => [id, el(id).type === 'checkbox' ? el(id).checked : el(id).value]));
}

let pending = 0;
async function update() {
  const token = ++pending;
  const p = policy();
  showInputs(p);
  el('status').textContent = `Loading ${p.shape}…`;
  await ensure(p.shape);
  if (token !== pending) return;
  const d = dataFor(p.shape);
  const t0 = performance.now();
  const res = evaluate(p, d);
  const s = summary(p, res);
  readouts(s);
  const floor = common.meta.anchorFloor?.[p.shape];
  el('caveat').textContent = p.rule === 'seed' ? 'One-seed candidates come from the background model.'
    : ['anchors', 'either', 'both'].includes(p.rule) && (p.method !== 'lr' || p.target !== 'fixed' || p.value < floor)
      ? `Background anchors on weak diagonals were kept only at ≥ ${floor} bits for this seed, so candidates read low for anchors below that.` : '';
  plotBars(p, res, d);
  plotStrength(p, res, d);
  plotFamilies(res);
  drawFamily(p, res, d, common);
  el('status').textContent = `Evaluated in ${(performance.now() - t0).toFixed(0)} ms; ${d.truth.rows.toLocaleString()} seeds in ${common.pairs.n.toLocaleString()} hits, ${d.bg.rows.rows.toLocaleString()} background seeds over ${d.meta.reads} reads.`;
  clearTimeout(update.t);
  update.t = setTimeout(() => plotSweep(p), 250);
}

async function compareShapes() {
  const p = policy();
  const rows = [];
  for (const shape of Object.keys(common.meta.shapes)) {
    el('status').textContent = `Loading ${shape}…`;
    await ensure(shape);
    const q = { ...p, shape };
    rows.push({ name: SHAPE_LABELS[shape], s: summary(q, evaluate(q, dataFor(shape))) });
  }
  el('shape-table').innerHTML = tableRows(rows, false);
  el('status').textContent = 'Compared at the same policy values; a longer seed at the same T holds a very different index, so read the Index column.';
}

async function main() {
  common = await loadCommon();
  const shapeSel = el('shape');
  for (const s of Object.keys(common.meta.shapes)) shapeSel.append(new Option(SHAPE_LABELS[s], s));
  shapeSel.value = 'k5';
  buildSlopes('seed-slopes', 'seed');
  buildSlopes('bar-slopes', 'bar');
  const feat = el('feature');
  for (const [f, spec] of Object.entries(FEATURES)) feat.append(new Option(spec.label, f));
  setValueRange('lr');
  el('provenance').textContent = common.provenance || '';
  const presets = el('presets');
  for (const [name, state] of Object.entries(PRESETS)) {
    const b = document.createElement('button');
    b.type = 'button'; b.textContent = name;
    b.onclick = () => { apply(state); update(); };
    presets.append(b);
  }
  el('method').addEventListener('change', () => setValueRange(el('method').value));
  el('controls').addEventListener('input', update);
  feat.addEventListener('input', update);
  el('scheme').addEventListener('input', update);
  el('save').onclick = () => {
    const p = policy();
    const name = el('name').value || `${p.shape} ${p.method} ${saved.length + 1}`;
    saved.push({ name, s: summary(p, evaluate(p, dataFor(p.shape))), form: formState(), policy: p });
    persist(); drawSaved();
  };
  el('saved').addEventListener('click', (e) => {
    const t = /** @type {HTMLElement} */ (e.target);
    if (t.dataset.load) { const st = saved[+t.dataset.load].form; setValueRange(st.method); for (const [k, v] of Object.entries(st)) { const c = el(k); if (c.type === 'checkbox') c.checked = v; else c.value = v; } update(); }
    if (t.dataset.del) { saved.splice(+t.dataset.del, 1); persist(); drawSaved(); }
  });
  el('export').onclick = () => {
    const a = document.createElement('a');
    a.href = URL.createObjectURL(new Blob([JSON.stringify(saved.length ? saved : [{ policy: policy() }], null, 1)], { type: 'application/json' }));
    a.download = 'kfpfam-policies.json';
    a.click();
  };
  el('shapes-all').onclick = compareShapes;
  setupFamily(common, update);
  el('worst').addEventListener('click', (e) => { const tr = /** @type {HTMLElement} */ (e.target).closest('tr'); if (tr?.dataset.model) { el('family').value = common.models.name[+tr.dataset.model]; update(); el('plot-family').scrollIntoView(); } });
  drawSaved();
  setupSketch();
  setupPromiscuous();
  await update();
}
main().catch((e) => { el('status').textContent = `Failed: ${e.message}`; throw e; });
void parseRows;
