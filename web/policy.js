// @ts-check
// Evaluate a k-mer policy from measured seeds, in the browser and in checks/test_policy.mjs.
// A seed is a k-mer of a read at a model window. It carries what any policy could ask of it, so
// changing a policy needs no index. See measure/scan and docs/plan.md.

export const T_GRID = Array.from({ length: 21 }, (_, i) => i + 4);
export const N_GRID = [1, 3, 10, 30, 100, 300, 1000, 3000, 10000];
export const FEATURES = {
  length: { label: 'Model length', log: true },
  ga: { label: 'Gathering threshold (bits)', log: false },
  nseq: { label: 'Sequences in seed alignment', log: true },
  effn: { label: 'Effective sequences', log: true },
  mean_info: { label: 'Mean information per column', log: false },
  top_info: { label: 'Best 5-column window (bits)', log: false },
  conserved: { label: 'Share of columns ≥ 1 bit', log: false },
  mu: { label: 'MSV μ', log: false },
};

/** Offset per model from a linear rule on standardised features: sum of slope × z-score. */
export function modelOffsets(models, slopes) {
  const n = models.acc.length;
  const out = new Float32Array(n);
  for (const [f, slope] of Object.entries(slopes || {})) {
    if (!slope) continue;
    const z = models.z[f];
    for (let i = 0; i < n; i++) out[i] += slope * z[i];
  }
  return out;
}

/** Decode the bytes of scan's write_rows. */
export function parseRows(buf) {
  const v = new DataView(buf);
  const groups = v.getUint32(0, true);
  const rows = v.getUint32(4, true);
  let p = 8;
  const take = (Type, n) => {
    const out = new Type(buf.slice(p, p + n * Type.BYTES_PER_ELEMENT));
    p += n * Type.BYTES_PER_ELEMENT;
    return out;
  };
  const offs = take(Uint32Array, groups + 1);
  return {
    groups, rows, offs,
    col: take(Uint16Array, rows), hash: take(Uint16Array, rows),
    j: take(Uint8Array, rows), score: take(Uint8Array, rows), diag: take(Uint8Array, rows),
    rankLR: take(Uint8Array, rows), rankP: take(Uint8Array, rows), logp: take(Uint8Array, rows),
    info: take(Uint8Array, rows), word: take(Uint8Array, rows),
  };
}

/** One model's value at a grid position, log-linear between grid points. Tables hold log2(x + 1e-12). */
function interp(tab, model, stat, x, grid, base) {
  const ng = tab.ng;
  const at = (i) => tab.data[(model * ng + base + i) * 3 + stat];
  if (x <= grid[0]) return at(0);
  const last = grid.length - 1;
  if (x >= grid[last]) return at(last);
  let i = 0;
  while (grid[i + 1] < x) i++;
  const f = (x - grid[i]) / (grid[i + 1] - grid[i]);
  return at(i) * (1 - f) + at(i + 1) * f;
}

const NB = N_GRID.length;
const NT = T_GRID.length;
const BASE = { lr: 0, topLR: NT, topP: NT + NB, sample: NT + 2 * NB };

/**
 * T per model. Base T plus offset; or solved per model for a budget of entries per column, or a
 * coverage of the model's own emissions.
 */
export function thresholds(policy, models, tab, offsets) {
  const n = models.acc.length;
  const T = new Float32Array(n);
  const target = policy.target || 'fixed';
  for (let m = 0; m < n; m++) {
    if (target === 'fixed') {
      T[m] = policy.value + offsets[m];
    } else {
      const want = target === 'budget'
        ? Math.log2(policy.budget * models.length[m] + 1e-12)
        : policy.coverage;
      let t = T_GRID[NT - 1];
      for (let i = 0; i < NT - 1; i++) {
        const stat = target === 'budget' ? 0 : 2;
        let hi = tab.data[(m * tab.ng + i) * 3 + stat];
        let lo = tab.data[(m * tab.ng + i + 1) * 3 + stat];
        if (target === 'coverage') {
          const w = Math.max(models.length[m], 1);
          hi = hi / w; lo = lo / w;
        }
        if (lo <= want) {
          t = T_GRID[i] + (hi === lo ? 0 : Math.max(0, Math.min(1, (hi - want) / (hi - lo))));
          break;
        }
      }
      T[m] = t + offsets[m];
    }
  }
  return T;
}

/** Index entries and background hits per query k-mer, summed over models. */
export function indexCost(policy, models, tab, T, offsets) {
  const n = models.acc.length;
  const perModelN = new Float32Array(n);
  const perModelF = new Float32Array(n);
  const keep = policy.windowKeep ?? 1;
  let entries = 0, f = 0;
  for (let m = 0; m < n; m++) {
    let e, q;
    if (policy.method === 'lr') {
      e = 2 ** interp(tab, m, 0, T[m], T_GRID, BASE.lr);
      q = 2 ** interp(tab, m, 1, T[m], T_GRID, BASE.lr);
    } else if (policy.method === 'words') {
      e = tab.words[m * 2];
      q = tab.words[m * 2 + 1];
    } else {
      const x = Math.max(1, policy.value * 2 ** -offsets[m]);
      e = 2 ** interp(tab, m, 0, x, N_GRID, BASE[policy.method]);
      q = 2 ** interp(tab, m, 1, x, N_GRID, BASE[policy.method]);
    }
    if (keep < 1 && tab.strata) {
      const t = policy.method === 'lr' ? T[m] : 10;
      const ratio = tab.strata(keep, t);
      e *= ratio.n; q *= ratio.f;
    }
    const s = policy.sparsity || 1;
    e /= s; q /= s;
    perModelN[m] = e; perModelF[m] = q;
    entries += e; f += q;
  }
  return { entries, f, perModelN, perModelF };
}

/** Does this seed belong to the policy's k-mer set? */
export function selector(policy, rows, modelOfGroup, T, offsets) {
  const keepInfo = Math.floor(255 * (1 - (policy.windowKeep ?? 1)));
  const sparsity = policy.sparsity || 1;
  const hashMax = 256 / sparsity;
  const { score, rankLR, rankP, logp, info, word, hash } = rows;
  const pre = (r) => info[r] >= keepInfo && (hash[r] >> 8) < hashMax;
  switch (policy.method) {
    case 'lr':
      return (r, m) => pre(r) && score[r] / 4 >= T[m];
    case 'topLR':
      return (r, m) => pre(r) && rankLR[r] / 8 <= Math.log2(Math.max(1, policy.value * 2 ** -offsets[m]));
    case 'topP':
      return (r, m) => pre(r) && rankP[r] / 8 <= Math.log2(Math.max(1, policy.value * 2 ** -offsets[m]));
    case 'sample':
      return (r, m) => {
        if (!pre(r)) return false;
        const n = Math.max(1, policy.value * 2 ** -offsets[m]);
        const p = 2 ** (-logp[r] / 4);
        return (hash[r] & 255) / 256 < 1 - (1 - p) ** n;
      };
    case 'words':
      return (r) => pre(r) && (word[r] & 1) === 1;
    default:
      throw new Error(policy.method);
  }
}

const usesAnchors = (rule) => rule === 'anchors' || rule === 'either' || rule === 'both';
const usesDiagonal = (rule) => rule === 'diag' || rule === 'either' || rule === 'both';

/**
 * Pseudo-alignment of a read and model from their anchors: the best-scoring chain of at least n
 * selected seeds, in read order, not overlapping, each within `window` diagonals of the last.
 * No hashing of distances, so an indel between anchors costs nothing. Returns the chain's summed
 * seed score, or 0.
 */
function chainScore(idx, rows, span, n, window, dp) {
  const { j, col, score } = rows;
  idx.sort((a, b) => j[a] - j[b]);
  const len = idx.length;
  if (dp.length < len * 3) dp.length = len * 3;
  let best = 0;
  for (let i = 0; i < len; i++) {
    const r = idx[i];
    const s = score[r] / 4;
    const d = col[r] - j[r];
    dp[i * 3] = s;
    dp[i * 3 + 1] = 0;
    dp[i * 3 + 2] = 0;
    for (let p = 0; p < i; p++) {
      const q = idx[p];
      if (j[r] - j[q] < span) continue;
      const delta = d - (col[q] - j[q]);
      if (delta > window || delta < -window) continue;
      if (dp[p * 3] + s > dp[i * 3 + 1]) dp[i * 3 + 1] = dp[p * 3] + s;
      if (dp[p * 3 + 1] > 0 && dp[p * 3 + 1] + s > dp[i * 3 + 2]) dp[i * 3 + 2] = dp[p * 3 + 1] + s;
    }
    const v = n === 2 ? Math.max(dp[i * 3 + 1], dp[i * 3 + 2]) : dp[i * 3 + 2];
    if (v > best) best = v;
  }
  return best;
}

/** Per group: best diagonal among selected seeds (0 if none), and the best anchor chain. */
export function groupEvidence(policy, rows, modelOfGroup, T, offsets, span) {
  const select = selector(policy, rows, modelOfGroup, T, offsets);
  const best = new Float32Array(rows.groups);
  const chain = new Float32Array(rows.groups);
  const { offs, diag } = rows;
  const needChain = usesAnchors(policy.rule);
  const n = policy.anchorN || 2;
  const window = policy.anchorWindow ?? 3;
  const idx = [];
  const dp = [];
  for (let g = 0; g < rows.groups; g++) {
    const m = modelOfGroup[g];
    idx.length = 0;
    for (let r = offs[g]; r < offs[g + 1]; r++) {
      if (!select(r, m)) continue;
      const d = diag[r] / 2;
      if (d > best[g]) best[g] = d;
      if (needChain) idx.push(r);
    }
    if (idx.length >= n) chain[g] = chainScore(idx, rows, span, n, window, dp);
  }
  return { best, chain };
}

/** Whether a group passes the read hit rule. */
export function passes(policy, ev, g, m, bars, bestAny, bestModel, models) {
  const seeded = ev.best[g] > 0;
  if (!seeded) return false;
  const diagonal = ev.best[g] >= bars[m];
  const anchors = ev.chain[g] > 0 && ev.chain[g] >= (policy.anchorSum || 0) + (bars[m] - policy.bar);
  let ok;
  if (policy.rule === 'seed') ok = true;
  else if (policy.rule === 'anchors') ok = anchors;
  else if (policy.rule === 'either') ok = diagonal || anchors;
  else if (policy.rule === 'both') ok = diagonal && anchors;
  else ok = diagonal;
  if (ok && policy.rankBits != null && bestAny) {
    const within = ev.best[g] >= bestAny[g] - policy.rankBits;
    const sameClan = policy.clanAware && models.clan[m] && models.clan[m] === models.clan[bestModel[g]];
    ok = within || sameClan;
  }
  return ok;
}

export function barsFor(policy, offsets) {
  const b = new Float32Array(offsets.length);
  for (let m = 0; m < b.length; m++) b[m] = policy.bar + offsets[m];
  return b;
}

const isTruth = (truth, pairs, g) => {
  if (truth === 'ga') return !!pairs.ga[g];
  if (truth === 'e3') return pairs.evalue[g] <= 1e-3;
  if (truth === 'e10') return pairs.evalue[g] <= 10;
  return true; // every pair: for fragments, their family even where hmmsearch misses it
};

export const E_BINS = [
  { label: '≤ 1e-10', top: 1e-10 }, { label: '1e-10 – 1e-5', top: 1e-5 }, { label: '1e-5 – 1e-3', top: 1e-3 },
  { label: '1e-3 – 1', top: 1 }, { label: '1 – 10', top: 10 }, { label: '10 – 1000', top: 1000 }, { label: 'not reported', top: Infinity },
];
export const SCORE_BINS = [
  { label: '≥ 60', low: 60 }, { label: '40 – 60', low: 40 }, { label: '30 – 40', low: 30 }, { label: '25 – 30', low: 25 },
  { label: '20 – 25', low: 20 }, { label: '15 – 20', low: 15 }, { label: '< 15', low: -900 }, { label: 'not reported', low: -Infinity },
];
export const GA_BINS = [
  { label: 'GA + 30 or more', low: 30 }, { label: 'GA + 10 to 30', low: 10 }, { label: 'GA to + 10', low: 0 },
  { label: 'GA − 5 to GA', low: -5 }, { label: 'GA − 10 to − 5', low: -10 }, { label: 'GA − 20 to − 10', low: -20 },
  { label: 'below GA − 20', low: -900 }, { label: 'not reported', low: -Infinity },
];
export const binsFor = (scheme) => ({ ga: GA_BINS, score: SCORE_BINS, evalue: E_BINS }[scheme]);
export function binOf(scheme, evalue, score, ga) {
  if (scheme === 'ga') return GA_BINS.findIndex((b) => (score <= -900 ? -Infinity : score - ga) >= b.low);
  if (scheme === 'score') return SCORE_BINS.findIndex((b) => score >= b.low);
  return E_BINS.findIndex((b) => evalue <= b.top);
}

/** Recall by match strength, and how candidates of a read set and its decoys spread over the bins. */
export function byStrength(res, policy, d, scheme) {
  const bins = binsFor(scheme);
  const out = bins.map((b) => ({ label: b.label, n: 0, weight: 0, ok: 0, real: 0, decoy: 0, hits: 0, decoyHits: 0 }));
  const want = policy.set || 'all';
  const { pairs } = d;
  // When hmmsearch rated more reads than were scanned, reported hits get their rates from all of them
  const rated = d.meta.truthReads;
  for (let g = 0; g < pairs.n; g++) {
    const isDecoy = rated && pairs.set[g] === `${want}decoy`;
    if (!isDecoy && (want !== 'all' ? pairs.set[g] !== want : !pairs.set[g].startsWith('f'))) continue;
    const b = out[binOf(scheme, pairs.evalue[g], pairs.score[g], d.models.ga[pairs.model[g]])];
    const w = pairs.weight[g];
    const ok = passes(policy, res.truthEv, g, pairs.model[g], res.bars, pairs.best, pairs.bestModel, d.models);
    if (isDecoy) {
      b.decoyHits += w / rated;
      if (ok) b.decoy += w / rated;
      continue;
    }
    b.n++; b.weight += w;
    if (rated) b.hits += w / rated;
    if (ok) { b.ok += w; if (rated) b.real += w / rated; }
  }
  const spread = (set, ev, key) => {
    for (let g = 0; g < set.rows.groups; g++) {
      const m = set.model[g];
      if (rated && set.score[g] > -900) continue;
      if (ev.best[g] && passes(policy, ev, g, m, res.bars, set.bestByGroup, set.bestModelByGroup, d.models)) {
        out[binOf(scheme, set.evalue[g], set.score[g], d.models.ga[m])][key] += 1 / set.reads;
      }
    }
  };
  spread(d.bg, res.bgEv, 'real');
  if (d.decoy) spread(d.decoy, res.decoyEv, 'decoy');
  return out;
}

const quantile = (sorted, q) => (sorted.length ? sorted[Math.min(sorted.length - 1, Math.floor(q * sorted.length))] : NaN);

/** Recall overall and per family, for pairs of the chosen set and truth. */
export function recall(policy, pairs, truthRows, ev, bars, models) {
  const fam = new Map();
  let tot = 0, ok = 0;
  const want = policy.set || 'all';
  for (let g = 0; g < pairs.n; g++) {
    if (want !== 'all' && pairs.set[g] !== want) continue;
    if (want === 'all' && !pairs.set[g].startsWith('f')) continue;
    if (!isTruth(policy.truth, pairs, g)) continue;
    const m = pairs.model[g];
    const pass = passes(policy, ev, g, m, bars, pairs.best, pairs.bestModel, models);
    const w = pairs.weight[g];
    tot += w;
    if (pass) ok += w;
    let f = fam.get(m);
    if (!f) fam.set(m, (f = { n: 0, ok: 0, w: 0, wok: 0 }));
    f.n++; f.w += w;
    if (pass) { f.ok++; f.wok += w; }
  }
  const rec = [];
  for (const [m, f] of fam) rec.push({ m, n: f.n, recall: f.ok / f.n, weight: f.w });
  const sorted = rec.map((r) => r.recall).sort((a, b) => a - b);
  // Tail figures only from families with 3+ hits: one or two hits make a recall of 0 or 1.
  const solid = rec.filter((r) => r.n >= 3).map((r) => r.recall).sort((a, b) => a - b);
  return {
    recall: ok / Math.max(tot, 1e-9),
    families: rec,
    familyMean: sorted.reduce((a, b) => a + b, 0) / Math.max(sorted.length, 1),
    familyP10: quantile(solid, 0.1),
    solidFamilies: solid.length,
    familyP01: quantile(solid, 0.01),
    zero: sorted.filter((x) => x === 0).length / Math.max(sorted.length, 1),
    below50: sorted.filter((x) => x < 0.5).length / Math.max(sorted.length, 1),
  };
}

/** Candidate models per read. */
export function candidates(policy, bg, ev, bars, cost, wRead, models) {
  if (policy.rule === 'seed' && policy.rankBits == null) {
    let c = 0;
    for (let m = 0; m < cost.perModelF.length; m++) c += 1 - Math.exp(-cost.perModelF[m] * wRead);
    return { mean: c, perModel: null, measured: false };
  }
  const perModel = new Float32Array(models.acc.length);
  let c = 0;
  for (let g = 0; g < bg.rows.groups; g++) {
    const m = bg.model[g];
    const read = bg.read[g];
    if (!ev.best[g]) continue;
    const ok = passes(policy, ev, g, m, bars, bg.bestByGroup, bg.bestModelByGroup, models);
    if (ok) { c++; perModel[m]++; void read; }
  }
  return { mean: c / bg.reads, perModel, measured: true };
}

/** Estimated milliseconds per read; each constant is an assumption until the prototype measures it. */
export function estimateMs(policy, hits, cands, wRead, c) {
  const lookups = wRead / (policy.sparsity || 1);
  return (lookups * c.lookup + hits * c.hit + cands * c.candidate) / 1e6;
}

export function evaluate(policy, d) {
  const offsets = modelOffsets(d.models, policy.seedSlopes);
  const barOffsets = modelOffsets(d.models, policy.barSlopes);
  const T = thresholds(policy, d.models, d.tab, offsets);
  const cost = indexCost(policy, d.models, d.tab, T, offsets);
  const bars = barsFor(policy, barOffsets);
  const truthEv = groupEvidence(policy, d.truth, d.pairs.model, T, offsets, d.span);
  const bgEv = groupEvidence(policy, d.bg.rows, d.bg.model, T, offsets, d.span);
  const rec = recall(policy, d.pairs, d.truth, truthEv, bars, d.models);
  const wRead = d.meta.residuesPerRead - (d.span - 1) * d.meta.orfsPerRead;
  const hits = cost.f * wRead;
  const cand = candidates(policy, d.bg, bgEv, bars, cost, wRead, d.models);
  const decoyEv = d.decoy ? groupEvidence(policy, d.decoy.rows, d.decoy.model, T, offsets, d.span) : null;
  const decoyCand = d.decoy && policy.rule !== 'seed' ? candidates(policy, d.decoy, decoyEv, bars, cost, wRead, d.models).mean : null;
  return { policy, T, offsets, barOffsets, bars, cost, rec, hits, cand, truthEv, bgEv, decoyEv, decoyCand, wRead };
}

/** Recall and candidates at other bars, from the same evidence (no re-selection). */
export function barSweep(res, policy, d, barsList) {
  return barsList.map((bar) => {
    const p = { ...policy, rule: usesDiagonal(policy.rule) ? policy.rule : 'diag', bar };
    const bars = barsFor(p, res.barOffsets);
    const rec = recall(p, d.pairs, d.truth, res.truthEv, bars, d.models);
    const cand = candidates(p, d.bg, res.bgEv, bars, res.cost, res.wRead, d.models);
    return { bar, recall: rec.recall, familyP10: rec.familyP10, familyMean: rec.familyMean, candidates: cand.mean };
  });
}
