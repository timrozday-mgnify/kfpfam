// @ts-check
/* global Plotly */
// Rejecting promiscuous k-mers on top of k4 T7. Numbers from measure/promiscuous.py.
const el = (id) => /** @type {HTMLInputElement} */ (document.getElementById(id));
const LAYOUT = { margin: { t: 34, r: 16, b: 48, l: 56 }, font: { size: 12 }, legend: { orientation: 'h', y: -0.28 } };
const RULES = {
  index: ['Longest index lists', '#7a9fd3'],
  expected: ['List length × background frequency', '#2c5898'],
  sample: ['Hits measured on this sample', '#c8553d'],
  repeat: ['Masking tandem repeats in reads', '#2a9d5c'],
};
let data;

const speed = (removed) => `${(1 / (1 - removed)).toFixed(2)}×`;

function draw() {
  const bar = el('prom-bar').value;
  const bin = +el('prom-bin').value;
  const base = data.base[bar];
  const series = (field, title, yTitle, id, range) => {
    const traces = Object.entries(RULES).map(([rule, [name, color]]) => {
      const rows = data.rules.filter((r) => r.rule === rule);
      return {
        x: [0, ...rows.map((r) => r.removed)], y: [field(base), ...rows.map((r) => field(r.bars[bar]))],
        text: ['no rejection', ...rows.map((r) => (r.period ? `repeat within ${r.period}` : `${r.kmers.toLocaleString()} 4-mers rejected`))],
        hoverinfo: 'text+x+y', mode: 'lines+markers', name, line: { width: 2, color }, marker: { size: 8, color },
      };
    });
    Plotly.react(id, traces, { ...LAYOUT, title: { text: title, font: { size: 13 } }, xaxis: { title: 'Share of hits removed', range: [0, 0.75], tickformat: '.0%' }, yaxis: { title: yTitle, ...range } }, { responsive: true });
  };
  series((r) => r.recall[bin], `Recall, ${data.bins[bin]}`, 'Recall', 'plot-prom-recall', { range: [0, 1.02] });
  series((r) => r.bray, 'Family abundance change against no rejection', 'Bray-Curtis', 'plot-prom-bray', {});
  const row = (name, removed, r) => `<tr><td class="text">${name}</td><td>${(removed * 100).toFixed(0)}%</td><td>${speed(removed)}</td>${r.recall.map((v) => `<td>${v.toFixed(3)}</td>`).join('')}<td>${r.candidates.toFixed(1)}</td><td>${r.decoy.toFixed(1)}</td><td>${r.p10.toFixed(2)}</td><td>${r.families_lost}</td><td>${r.bray.toFixed(3)}</td><td class="text">${r.worst.slice(0, 3).map((w) => `${w[0]} ${w[1]}→${w[2]}`).join(', ')}</td></tr>`;
  const head = `<tr><th>Rule</th><th>Hits removed</th><th>Speed, seeding</th>${data.bins.map((b) => `<th>${b}</th>`).join('')}<th>Candidates per read</th><th>On reversed reads</th><th>Family p10</th><th>Families lost</th><th>Bray-Curtis</th><th>Most hits lost</th></tr>`;
  const rows = [row('No rejection', 0, base), ...data.rules.map((r) => row(`${RULES[r.rule][0]}${r.period ? `, within ${r.period}` : `, ${r.kmers.toLocaleString()} 4-mers`}`, r.removed, r.bars[bar]))];
  el('prom-table').innerHTML = head + rows.join('');
}

export async function setupPromiscuous() {
  const res = await fetch('data/promiscuous.json');
  if (!res.ok) { el('prom-table').textContent = 'data/promiscuous.json is missing; run measure/promiscuous.py.'; return; }
  data = await res.json();
  for (const b of data.bars) el('prom-bar').append(new Option(`Diagonal ≥ ${b} bits`, String(b)));
  data.bins.forEach((b, i) => el('prom-bin').append(new Option(b, String(i))));
  for (const id of ['prom-bar', 'prom-bin']) el(id).addEventListener('input', draw);
  draw();
}
