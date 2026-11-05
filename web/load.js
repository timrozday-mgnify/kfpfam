// @ts-check
import { FEATURES, parseRows } from './policy.js';

const isNode = typeof window === 'undefined';
/** Bytes of a file, gunzipped natively if it ends in .gz. */
async function bytes(name) {
  if (isNode) {
    const fs = await import('node:fs');
    const zlib = await import('node:zlib');
    const url = new URL(`data/${name}`, import.meta.url);
    const b = fs.readFileSync(url);
    const out = name.endsWith('.gz') ? zlib.gunzipSync(b) : b;
    return out.buffer.slice(out.byteOffset, out.byteOffset + out.byteLength);
  }
  const r = await fetch(`data/${name}`);
  if (!r.ok) throw new Error(`${name}: ${r.status}`);
  if (!name.endsWith('.gz')) return r.arrayBuffer();
  return new Response(r.body.pipeThrough(new DecompressionStream('gzip'))).arrayBuffer();
}
const json = async (name) => JSON.parse(new TextDecoder().decode(await bytes(name)));

/** A model's log-odds in bits, length x 20, from its piece of the profile. */
const pieces = new Map();
export async function loadProfile(m, models) {
  const n = models.acc.length;
  const CHUNKS = 200;
  let c = Math.floor((m * CHUNKS) / n);
  while (Math.floor(((c + 1) * n) / CHUNKS) <= m) c++;
  while (Math.floor((c * n) / CHUNKS) > m) c--;
  if (!pieces.has(c)) pieces.set(c, new Int8Array(await bytes(`profile/${c}.bin.gz`)));
  const piece = pieces.get(c);
  let at = 0;
  for (let i = Math.floor((c * n) / CHUNKS); i < m; i++) at += models.length[i];
  return piece.subarray(at * 20, (at + models.length[m]) * 20);
}

export async function loadCommon() {
  const [models, rawPairs, meta, baselines] = await Promise.all([json('models.json.gz'), json('pairs.json.gz'), json('meta.json'), json('baselines.json')]);
  models.z = {};
  for (const [f, spec] of Object.entries(FEATURES)) {
    const v = models[f].map((x) => (spec.log ? Math.log2(Math.max(x, 1e-6)) : x));
    const mean = v.reduce((a, b) => a + b, 0) / v.length;
    const sd = Math.sqrt(v.reduce((a, b) => a + (b - mean) ** 2, 0) / v.length) || 1;
    models.z[f] = Float32Array.from(v, (x) => (x - mean) / sd);
  }
  const pairs = { ...rawPairs, n: rawPairs.model.length, model: Uint16Array.from(rawPairs.model), weight: Float32Array.from(rawPairs.weight), best: Float32Array.from(rawPairs.best), bestModel: Uint16Array.from(rawPairs.bestModel) };
  // Real reads have their own background; absent until measure/run_real_policy.sh has run
  const realMeta = await json('real/meta.json').catch(() => null);
  return { models, pairs, meta, realMeta, baselines: baselines.tools, byBin: baselines.byBin, hmmsearchMs: baselines.hmmsearch_ms_per_read, provenance: baselines.provenance };
}

/** Background rows per (read, model), with each one's rating by hmmsearch if it was reported. */
function backgroundSet(rows, g, reads, best, bestModel, hits) {
  const read = new Uint16Array(rows.groups);
  const model = new Uint16Array(rows.groups);
  const bestByGroup = new Float32Array(rows.groups);
  const bestModelByGroup = new Uint16Array(rows.groups);
  const evalue = new Float32Array(rows.groups).fill(1e9);
  const score = new Float32Array(rows.groups).fill(-999);
  const at = new Map();
  for (let i = 0; i < rows.groups; i++) {
    read[i] = g[i * 2];
    model[i] = g[i * 2 + 1];
    bestByGroup[i] = best[read[i]];
    bestModelByGroup[i] = bestModel[read[i]];
    at.set(read[i] * 40000 + model[i], i);
  }
  for (const [r, m, e, s] of hits || []) {
    const i = at.get(r * 40000 + m);
    if (i !== undefined) { evalue[i] = e; score[i] = s; }
  }
  return { rows, read, model, reads, bestByGroup, bestModelByGroup, evalue, score };
}

/** One read set's background and decoys. */
async function loadReads(shape, meta, dir = '') {
  const [bgBuf, groupBuf] = await Promise.all([`${dir}${shape}.bg.gz`, `${dir}${shape}.bggroups.gz`].map(bytes));
  const bg = backgroundSet(parseRows(bgBuf), new Uint16Array(groupBuf), meta.reads, meta.readBest, meta.readBestModel, meta.bgHits);
  if (!meta.decoyHits) return { bg, decoy: null };
  const [dBuf, dGroups] = await Promise.all([`${dir}${shape}.decoy.gz`, `${dir}${shape}.decoygroups.gz`].map(bytes));
  return { bg, decoy: backgroundSet(parseRows(dBuf), new Uint16Array(dGroups), meta.reads, meta.decoyBest, meta.decoyBestModel, meta.decoyHits) };
}

export const loadReal = (shape, common) => loadReads(shape, common.realMeta, 'real/');

/** The data a policy is evaluated on: real reads bring their own background, decoys and rates. */
export const withSet = (common, shape, set) => (set === 'real' && shape.real
  ? { ...common, ...shape, ...shape.real, meta: common.realMeta }
  : { ...common, ...shape, meta: common.meta });

export async function loadShape(shape, common) {
  const [truthBuf, tabBuf, wordsBuf, strataBuf, reads] = await Promise.all([
    ...[`${shape}.rows.gz`, `${shape}.tab.gz`, `${shape}.words.bin`, `${shape}.strata.bin`].map(bytes),
    loadReads(shape, common.meta),
  ]);
  const nm = common.models.acc.length;
  const ng = 21 + 27;
  const logs = new Int16Array(tabBuf, 0, nm * ng * 2);
  const cov = new Uint16Array(tabBuf, nm * ng * 4, nm * ng);
  const data = new Float32Array(nm * ng * 3);
  for (let i = 0; i < nm * ng; i++) {
    data[i * 3] = logs[i * 2] / 8;
    data[i * 3 + 1] = logs[i * 2 + 1] / 8;
    data[i * 3 + 2] = cov[i];
  }
  const words = new Float32Array(nm * 2);
  const wv = new DataView(wordsBuf);
  for (let m = 0; m < nm; m++) {
    words[m * 2] = wv.getUint32(m * 8, true);
    words[m * 2 + 1] = wv.getFloat32(m * 8 + 4, true);
  }
  const strata = new Float32Array(strataBuf);
  const { fracs } = common.meta;
  const tab = {
    ng, data, words,
    strata: (keep, t) => {
      const fi = Math.max(0, fracs.indexOf(keep));
      const ti = Math.max(0, Math.min(20, Math.round(t) - 4));
      return { n: strata[(fi * 21 + ti) * 2], f: strata[(fi * 21 + ti) * 2 + 1] };
    },
  };
  const shapeCols = common.meta.shapes[shape];
  return {
    tab,
    truth: parseRows(truthBuf),
    ...reads,
    span: shapeCols[shapeCols.length - 1] + 1,
  };
}
