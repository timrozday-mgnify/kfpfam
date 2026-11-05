//! Whole-Pfam seed measurements for the dashboard. No index: a k-mer is "in the index" for a
//! model exactly when a policy would select it at some window, so every seed in a read can be
//! recorded with what a policy needs (score, ranks in its window, emission probability, hash,
//! window information, diagonal score) and policies are evaluated later, in the browser.
//!
//! scan tables  <dir> <out>             per-model entries, background and coverage of each method
//! scan truth   <dir> <pairs.tsv> <out> seeds of each (sequence, model) pair
//! scan bg      <dir> <orfs.tsv> <bar> <out>  seeded diagonals >= bar over all models, per read
//! scan best    <dir> <seqs.tsv> <out>  best ungapped diagonal score of each sequence, over all models
//! <dir> holds profile.f32 (states x 21 log-odds bits), offs.u32 and seed.gz (Pfam-A.seed).
use flate2::{read::MultiGzDecoder, write::GzEncoder, Compression};
use rayon::prelude::*;
use std::{
    collections::{HashMap, HashSet},
    fs,
    io::{BufRead, BufReader, Write},
};

const A: &[u8] = b"ACDEFGHIKLMNPQRSTVWY";
const Q: [f64; 20] = [
    0.0787945, 0.0151600, 0.0535222, 0.0668298, 0.0397062, 0.0695071, 0.0229198, 0.0590092,
    0.0594422, 0.0963728, 0.0237718, 0.0414386, 0.0482904, 0.0395639, 0.0540978, 0.0683364,
    0.0540687, 0.0673417, 0.0114135, 0.0304133,
];
/// (name, columns, lowest score recorded). Spaced shapes and fixed-distance pairs are one thing.
const SHAPES: &[(&str, &[usize], f32)] = &[
    ("k4", &[0, 1, 2, 3], 6.0),
    ("k5", &[0, 1, 2, 3, 4], 6.0),
    ("k6", &[0, 1, 2, 3, 4, 5], 8.0),
    ("k7", &[0, 1, 2, 3, 4, 5, 6], 8.0),
    ("s5", &[0, 1, 3, 5, 6], 6.0),
    ("s6", &[0, 1, 2, 4, 6, 7], 8.0),
    ("p33", &[0, 1, 2, 6, 7, 8], 8.0),
];
pub const T_GRID: std::ops::RangeInclusive<i32> = 4..=24;
const N_GRID: [f64; 9] = [1., 3., 10., 30., 100., 300., 1000., 3000., 10000.];

struct Pfam {
    lo: Vec<f32>,
    offs: Vec<u32>,
    nbg: [f32; 20],
}
impl Pfam {
    fn load(dir: &str) -> Pfam {
        let b = fs::read(format!("{dir}/profile.f32")).unwrap();
        let lo = b.chunks_exact(4).map(|c| f32::from_le_bytes(c.try_into().unwrap())).collect();
        let b = fs::read(format!("{dir}/offs.u32")).unwrap();
        let offs = b.chunks_exact(4).map(|c| u32::from_le_bytes(c.try_into().unwrap())).collect();
        let mut qs = Q;
        let t: f64 = qs.iter().sum();
        qs.iter_mut().for_each(|q| *q /= t);
        Pfam { lo, offs, nbg: qs.map(|q| q as f32) }
    }
    fn models(&self) -> usize {
        self.offs.len() - 1
    }
    #[inline]
    fn range(&self, m: usize) -> (usize, usize) {
        (self.offs[m] as usize, self.offs[m + 1] as usize)
    }
}
fn encode(s: &[u8]) -> Vec<u8> {
    let mut t = [20u8; 256];
    for (i, c) in A.iter().enumerate() {
        t[*c as usize] = i as u8;
    }
    s.iter().map(|c| t[*c as usize]).collect()
}

const LR_BINS: usize = 640; // 0.25 bit, -40..120
const LR_ZERO: i32 = 160;
const P_BINS: usize = 400; // -log2 p in 0.25 bit
/// Binned distribution of a window's k-mers by log-odds and by emission probability.
struct Tab {
    lr: [Vec<f64>; 3],
    p: [Vec<f64>; 3],
}
#[inline]
fn lr_bin(s: f64) -> usize {
    ((s * 4.0).round() as i32 + LR_ZERO).clamp(0, LR_BINS as i32 - 1) as usize
}
#[inline]
fn p_bin(lp: f64) -> usize {
    ((lp * 4.0).round() as usize).min(P_BINS - 1)
}
fn window_tab(pf: &Pfam, st: usize, shape: &[usize], want_full: bool) -> Tab {
    let mut lr = [vec![0.0; LR_BINS], vec![0.0; LR_BINS], vec![0.0; LR_BINS]];
    let mut p = [vec![0.0; P_BINS], vec![0.0; P_BINS], vec![0.0; P_BINS]];
    lr[0][LR_ZERO as usize] = 1.0;
    lr[1][LR_ZERO as usize] = 1.0;
    lr[2][LR_ZERO as usize] = 1.0;
    p[0][0] = 1.0;
    p[1][0] = 1.0;
    p[2][0] = 1.0;
    let mut nl = [vec![0.0; LR_BINS], vec![0.0; LR_BINS], vec![0.0; LR_BINS]];
    let mut np = [vec![0.0; P_BINS], vec![0.0; P_BINS], vec![0.0; P_BINS]];
    for &o in shape {
        let row = &pf.lo[(st + o) * 21..(st + o) * 21 + 20];
        nl.iter_mut().for_each(|a| a.fill(0.0));
        np.iter_mut().for_each(|a| a.fill(0.0));
        for a in 0..20 {
            let s = row[a] as f64;
            let q = pf.nbg[a] as f64;
            let pp = q * 2f64.powf(s);
            let sh = (s * 4.0).round() as i32;
            let lp = (-(pp.log2()) * 4.0).round().max(0.0) as usize;
            for b in 0..LR_BINS {
                if lr[0][b] == 0.0 {
                    continue;
                }
                let t = (b as i32 + sh).clamp(0, LR_BINS as i32 - 1) as usize;
                nl[0][t] += lr[0][b];
                nl[1][t] += lr[1][b] * q;
                nl[2][t] += lr[2][b] * pp;
            }
            if want_full {
                for b in 0..P_BINS {
                    if p[0][b] == 0.0 {
                        continue;
                    }
                    let t = (b + lp).min(P_BINS - 1);
                    np[0][t] += p[0][b];
                    np[1][t] += p[1][b] * q;
                    np[2][t] += p[2][b] * pp;
                }
            } else {
                for b in 0..P_BINS {
                    if p[0][b] == 0.0 {
                        continue;
                    }
                    np[0][(b + lp).min(P_BINS - 1)] += p[0][b];
                }
            }
        }
        std::mem::swap(&mut lr, &mut nl);
        std::mem::swap(&mut p, &mut np);
    }
    Tab { lr, p }
}
fn write_gz(path: &str, b: &[u8]) {
    let mut e = GzEncoder::new(fs::File::create(path).unwrap(), Compression::new(6));
    e.write_all(b).unwrap();
    e.finish().unwrap();
}

/// Per model × grid: entries, summed background hit rate, mean coverage of the emissions.
/// Grid: T 4..24; top-N by log-odds; top-N by probability; n samples from emissions.
fn tables(pf: &Pfam, out: &str) {
    let ts: Vec<f64> = T_GRID.map(|t| t as f64).collect();
    let ng = ts.len() + 3 * N_GRID.len();
    fs::create_dir_all(out).unwrap();
    for (name, shape, _) in SHAPES {
        let span = shape.last().unwrap() + 1;
        let t0 = std::time::Instant::now();
        // keep windows by information: everything / top half, quarter, tenth, within each model
        let fracs = [1.0, 0.5, 0.25, 0.1];
        let res: Vec<(Vec<f64>, Vec<f64>)> = (0..pf.models())
            .into_par_iter()
            .map(|m| {
                let (a, b) = pf.range(m);
                let mut row = vec![0.0f64; ng * 3];
                let mut strata = vec![0.0f64; fracs.len() * ts.len() * 3];
                if b - a < span {
                    return (row, strata);
                }
                let wins = b - a - span + 1;
                let mut order: Vec<(f64, usize)> = (0..wins)
                    .map(|w| {
                        let mut info = 0.0;
                        for &o in *shape {
                            let r = &pf.lo[(a + w + o) * 21..(a + w + o) * 21 + 20];
                            for c in 0..20 {
                                let lo = r[c] as f64;
                                info += pf.nbg[c] as f64 * 2f64.powf(lo) * lo;
                            }
                        }
                        (-info, w)
                    })
                    .collect();
                order.sort_by(|x, y| x.partial_cmp(y).unwrap());
                for (rank, (_, w)) in order.iter().enumerate() {
                    let tab = window_tab(pf, a + w, shape, true);
                    // cumulative from the top of the log-odds range
                    let mut cum = [vec![0.0; LR_BINS + 1], vec![0.0; LR_BINS + 1], vec![0.0; LR_BINS + 1]];
                    for i in 0..3 {
                        for b2 in (0..LR_BINS).rev() {
                            cum[i][b2] = cum[i][b2 + 1] + tab.lr[i][b2];
                        }
                    }
                    for (ti, t) in ts.iter().enumerate() {
                        let bi = lr_bin(*t);
                        for i in 0..3 {
                            row[ti * 3 + i] += cum[i][bi];
                            for (fi, fr) in fracs.iter().enumerate() {
                                if (rank as f64) < fr * wins as f64 {
                                    strata[(fi * ts.len() + ti) * 3 + i] += cum[i][bi];
                                }
                            }
                        }
                    }
                    let base = ts.len();
                    // top N by log-odds
                    for (ni, n) in N_GRID.iter().enumerate() {
                        let mut left = *n;
                        let (mut c, mut qq, mut pp) = (0.0, 0.0, 0.0);
                        for b2 in (0..LR_BINS).rev() {
                            let cnt = tab.lr[0][b2];
                            if cnt <= 0.0 {
                                continue;
                            }
                            let f = (left / cnt).min(1.0);
                            c += cnt * f;
                            qq += tab.lr[1][b2] * f;
                            pp += tab.lr[2][b2] * f;
                            left -= cnt * f;
                            if left <= 0.0 {
                                break;
                            }
                        }
                        let k = (base + ni) * 3;
                        row[k] += c;
                        row[k + 1] += qq;
                        row[k + 2] += pp;
                        // top N by probability
                        let mut left = *n;
                        let (mut c, mut qq, mut pp) = (0.0, 0.0, 0.0);
                        for b2 in 0..P_BINS {
                            let cnt = tab.p[0][b2];
                            if cnt <= 0.0 {
                                continue;
                            }
                            let f = (left / cnt).min(1.0);
                            c += cnt * f;
                            qq += tab.p[1][b2] * f;
                            pp += tab.p[2][b2] * f;
                            left -= cnt * f;
                            if left <= 0.0 {
                                break;
                            }
                        }
                        let k = (base + N_GRID.len() + ni) * 3;
                        row[k] += c;
                        row[k + 1] += qq;
                        row[k + 2] += pp;
                        // n samples: P(in the sampled set) is 1-(1-p)^n
                        let (mut c, mut qq, mut pp) = (0.0, 0.0, 0.0);
                        for b2 in 0..P_BINS {
                            let cnt = tab.p[0][b2];
                            if cnt <= 0.0 {
                                continue;
                            }
                            let pr = 2f64.powf(-(b2 as f64) / 4.0);
                            let inc = 1.0 - (1.0 - pr.min(1.0 - 1e-12)).powf(*n);
                            c += cnt * inc;
                            qq += tab.p[1][b2] * inc;
                            pp += tab.p[2][b2] * inc;
                        }
                        let k = (base + 2 * N_GRID.len() + ni) * 3;
                        row[k] += c;
                        row[k + 1] += qq;
                        row[k + 2] += pp;
                    }
                }
                (row, strata)
            })
            .collect();
        let mut b = Vec::with_capacity(res.len() * ng * 12);
        for (row, _) in &res {
            for v in row {
                b.extend_from_slice(&(*v as f32).to_le_bytes());
            }
        }
        write_gz(&format!("{out}/{name}.tab.f32.gz"), &b);
        let mut tot = vec![0.0; fracs.len() * ts.len() * 3];
        for (_, s) in &res {
            for (t, v) in tot.iter_mut().zip(s) {
                *t += v;
            }
        }
        let mut b = Vec::new();
        for v in &tot {
            b.extend_from_slice(&(*v as f64).to_le_bytes());
        }
        fs::write(format!("{out}/{name}.strata.f64"), b).unwrap();
        eprintln!("{name}: {:.0}s", t0.elapsed().as_secs_f32());
    }
}

/// Background seeds on weak diagonals are kept only as anchors: at least this score, with a
/// non-overlapping partner within MAX_WINDOW diagonals, so window rules can be evaluated later.
const ANCHOR_FLOOR: [f32; 7] = [9.0, 9.5, 10.5, 11.5, 9.5, 10.5, 10.5];
const MAX_WINDOW: i64 = 6;

/// One seed in a sequence (see the module docs); packed for the browser later.
#[derive(Clone)]
struct Row {
    group: u32,
    col: u16,
    j: u8,
    score: f32,
    diag: f32,
    kmer: u64,
    rank_lr: f32,
    rank_p: f32,
    logp: f32,
    info: f32,
    word: bool,
    anchor: bool,
}
#[inline]
fn kadane(v: &[f32]) -> f32 {
    let (mut c, mut b) = (0.0f32, 0.0f32);
    for &x in v {
        c = (c + x).max(0.0);
        if c > b {
            b = c;
        }
    }
    b
}
fn kmer_of(f: &[u8], at: usize, shape: &[usize]) -> Option<u64> {
    let mut h = 0u64;
    for &o in shape {
        let c = f[at + o];
        if c >= 20 {
            return None;
        }
        h = h * 20 + c as u64;
    }
    Some(h)
}
/// Stable 16-bit hash of a k-mer (a policy's sparsity and sampling draws).
#[inline]
fn hash16(k: u64) -> u16 {
    let mut x = k.wrapping_mul(0x9E37_79B9_7F4A_7C15);
    x ^= x >> 29;
    x = x.wrapping_mul(0xBF58_476D_1CE4_E5B9);
    (x >> 48) as u16
}
/// Percentile of each window's information within its model (0..255), per shape.
fn info_percentiles(pf: &Pfam, m: usize, shape: &[usize]) -> Vec<u8> {
    let (a, b) = pf.range(m);
    let span = shape.last().unwrap() + 1;
    if b - a < span {
        return vec![];
    }
    let wins = b - a - span + 1;
    let mut v: Vec<(f64, usize)> = (0..wins)
        .map(|w| {
            let mut info = 0.0;
            for &o in shape {
                let r = &pf.lo[(a + w + o) * 21..(a + w + o) * 21 + 20];
                for c in 0..20 {
                    info += pf.nbg[c] as f64 * 2f64.powf(r[c] as f64) * r[c] as f64;
                }
            }
            (info, w)
        })
        .collect();
    v.sort_by(|x, y| x.partial_cmp(y).unwrap());
    let mut out = vec![0u8; wins];
    for (rank, (_, w)) in v.iter().enumerate() {
        out[*w] = (255 * rank / wins.max(1)) as u8;
    }
    out
}
fn fill_ranks(pf: &Pfam, a: usize, r: &mut Row, shape: &[usize], f: &[u8], memo: &mut HashMap<u16, Tab>) {
    let t = memo.entry(r.col).or_insert_with(|| window_tab(pf, a + r.col as usize, shape, false));
    let bi = lr_bin(r.score as f64);
    r.rank_lr = t.lr[0][bi..].iter().sum::<f64>().max(1.0) as f32;
    let mut lp = 0.0f64;
    for &o in shape {
        let c = f[r.j as usize + o] as usize;
        let s = pf.lo[(a + r.col as usize + o) * 21 + c] as f64;
        lp += -(pf.nbg[c] as f64 * 2f64.powf(s)).log2();
    }
    r.logp = lp as f32;
    r.rank_p = t.p[0][..=p_bin(lp)].iter().sum::<f64>().max(1.0) as f32;
}
/// All diagonals of one model against a sequence. Seeds are recorded where the diagonal scores
/// at least `bar`; returns the best diagonal score.
#[allow(clippy::too_many_arguments)]
fn scan_model(pf: &Pfam, m: usize, f: &[u8], si: usize, bar: f32, group: u32, out: &mut Vec<Row>, info: Option<&[u8]>, buf: &mut Vec<f32>, anchors: bool) -> f32 {
    let (a, b) = pf.range(m);
    let big = b - a;
    let l = f.len();
    let (_, shape, tmin) = SHAPES[si];
    let span = shape.last().unwrap() + 1;
    let mut best = 0.0f32;
    let mut memo = HashMap::new();
    let floor = ANCHOR_FLOOR[si];
    // (diagonal, read position, score, diagonal score, emitted already)
    let mut anchor_list: Vec<(i64, usize, f32, f32, bool)> = vec![];
    for d in (1 - l as i64)..(big as i64) {
        let x0 = (-d).max(0) as usize;
        let x1 = l.min((big as i64 - d) as usize);
        if x1 <= x0 {
            continue;
        }
        buf.clear();
        for x in x0..x1 {
            let col = (d + x as i64) as usize;
            buf.push(pf.lo[(a + col) * 21 + f[x] as usize]);
        }
        let ds = kadane(buf);
        if ds > best {
            best = ds;
        }
        if buf.len() < span {
            continue;
        }
        if ds < bar {
            if anchors {
                for p in 0..=(buf.len() - span) {
                    let mut s = 0.0;
                    for &o in shape {
                        s += buf[p + o];
                    }
                    if s >= floor {
                        anchor_list.push((d, x0 + p, s, ds, false));
                    }
                }
            }
            continue;
        }
        for p in 0..=(buf.len() - span) {
            let mut s = 0.0;
            for &o in shape {
                s += buf[p + o];
            }
            if s >= tmin {
                let j = x0 + p;
                let col = (d + j as i64) as usize;
                let Some(kmer) = kmer_of(f, j, shape) else { continue };
                let mut r = Row { group, col: col as u16, j: j as u8, score: s, diag: ds, kmer, rank_lr: 0.0, rank_p: 0.0, logp: 0.0, info: info.map_or(0.0, |v| v[col] as f32), word: false, anchor: false };
                fill_ranks(pf, a, &mut r, shape, f, &mut memo);
                out.push(r);
                if anchors && s >= floor {
                    anchor_list.push((d, j, s, ds, true));
                }
            }
        }
    }
    // Weak-diagonal anchors with a partner in the window join the record.
    for i in 0..anchor_list.len() {
        let (d, j, s, ds, emitted) = anchor_list[i];
        if emitted {
            continue;
        }
        let partner = anchor_list.iter().any(|&(d2, j2, ..)| (d2 - d).abs() <= MAX_WINDOW && j2.abs_diff(j) >= span);
        if !partner {
            continue;
        }
        let col = (d + j as i64) as usize;
        let Some(kmer) = kmer_of(f, j, shape) else { continue };
        let mut r = Row { group, col: col as u16, j: j as u8, score: s, diag: ds, kmer, rank_lr: 0.0, rank_p: 0.0, logp: 0.0, info: info.map_or(0.0, |v| v[col] as f32), word: false, anchor: true };
        fill_ranks(pf, a, &mut r, shape, f, &mut memo);
        out.push(r);
    }
    best
}
/// Pfam-A.seed k-mers per family, by accession (without version).
fn seed_words(dir: &str) -> HashMap<String, Vec<Vec<u8>>> {
    let r = BufReader::new(MultiGzDecoder::new(fs::File::open(format!("{dir}/seed.gz")).unwrap()));
    let mut out = HashMap::new();
    let mut acc = String::new();
    let mut seqs: Vec<Vec<u8>> = vec![];
    for line in r.split(b'\n') {
        let line = line.unwrap();
        if line.starts_with(b"#=GF AC") {
            let s = String::from_utf8_lossy(&line[7..]);
            acc = s.trim().split('.').next().unwrap().to_string();
        } else if line.starts_with(b"//") {
            out.insert(std::mem::take(&mut acc), std::mem::take(&mut seqs));
        } else if !line.is_empty() && line[0] != b'#' {
            if let Some(p) = line.iter().rposition(|c| *c == b' ') {
                let s: Vec<u8> = line[p + 1..].iter().filter(|c| c.is_ascii_alphabetic()).map(|c| c.to_ascii_uppercase()).collect();
                seqs.push(encode(&s));
            }
        }
    }
    out
}
fn accs(dir: &str) -> Vec<String> {
    fs::read_to_string(format!("{dir}/accs.txt")).unwrap().lines().map(|s| s.to_string()).collect()
}
fn family_words(seqs: &[Vec<u8>], shape: &[usize]) -> HashSet<u64> {
    let span = shape.last().unwrap() + 1;
    let mut set = HashSet::new();
    for s in seqs {
        if s.len() < span {
            continue;
        }
        for at in 0..=(s.len() - span) {
            if let Some(k) = kmer_of(s, at, shape) {
                set.insert(k);
            }
        }
    }
    set
}
/// A k-mer of the seed alignment at the window where it scores best, if the diagonal clears bar.
fn word_row(pf: &Pfam, m: usize, f: &[u8], j: usize, shape: &[usize], kmer: u64, group: u32, bar: f32, info: &[u8], buf: &mut Vec<f32>, memo: &mut HashMap<u16, Tab>) -> Option<Row> {
    let (a, b) = pf.range(m);
    let span = shape.last().unwrap() + 1;
    if b - a < span {
        return None;
    }
    let (mut bs, mut bc) = (f32::MIN, 0);
    for c in 0..=(b - a - span) {
        let mut s = 0.0;
        for &o in shape {
            s += pf.lo[(a + c + o) * 21 + f[j + o] as usize];
        }
        if s > bs {
            bs = s;
            bc = c;
        }
    }
    let d = bc as i64 - j as i64;
    buf.clear();
    for x in 0..f.len() {
        let col = d + x as i64;
        if col >= 0 && (col as usize) < b - a {
            buf.push(pf.lo[(a + col as usize) * 21 + f[x] as usize]);
        }
    }
    let ds = kadane(buf);
    if ds < bar {
        return None;
    }
    let mut r = Row { group, col: bc as u16, j: j as u8, score: bs, diag: ds, kmer, rank_lr: 0., rank_p: 0., logp: 0., info: info[bc] as f32, word: true, anchor: false };
    fill_ranks(pf, a, &mut r, shape, f, memo);
    Some(r)
}
/// groups u32, rows u32, group offsets u32[g+1]; columns col u16, hash u16, j, score(0.25), diag(0.5),
/// rankLR, rankP (log2*8), logp (*4), info, flags (1: in the seed alignment; 2: kept as an anchor on a
/// weak diagonal): one byte each.
fn write_rows(path: &str, rows: &mut [Row], groups: usize) {
    rows.sort_by_key(|r| (r.group, r.col, r.j));
    let mut offs = vec![0u32; groups + 1];
    for r in rows.iter() {
        offs[r.group as usize + 1] += 1;
    }
    for i in 0..groups {
        offs[i + 1] += offs[i];
    }
    let mut b = Vec::with_capacity(rows.len() * 12 + groups * 4);
    b.extend_from_slice(&(groups as u32).to_le_bytes());
    b.extend_from_slice(&(rows.len() as u32).to_le_bytes());
    offs.iter().for_each(|o| b.extend_from_slice(&o.to_le_bytes()));
    rows.iter().for_each(|r| b.extend_from_slice(&r.col.to_le_bytes()));
    rows.iter().for_each(|r| b.extend_from_slice(&hash16(r.kmer).to_le_bytes()));
    rows.iter().for_each(|r| b.push(r.j));
    rows.iter().for_each(|r| b.push((r.score * 4.0).round().clamp(0.0, 255.0) as u8));
    rows.iter().for_each(|r| b.push((r.diag * 2.0).round().clamp(0.0, 255.0) as u8));
    rows.iter().for_each(|r| b.push(((r.rank_lr as f64).log2() * 8.0).round().clamp(0.0, 255.0) as u8));
    rows.iter().for_each(|r| b.push(((r.rank_p as f64).log2() * 8.0).round().clamp(0.0, 255.0) as u8));
    rows.iter().for_each(|r| b.push((r.logp * 4.0).round().clamp(0.0, 255.0) as u8));
    rows.iter().for_each(|r| b.push(r.info as u8));
    rows.iter().for_each(|r| b.push(r.word as u8 | (r.anchor as u8) << 1));
    write_gz(path, &b);
}
fn read_tsv(p: &str) -> Vec<Vec<String>> {
    fs::read_to_string(p).unwrap().lines().map(|l| l.split('\t').map(|s| s.to_string()).collect()).collect()
}
/// Input: seq, model index, per line; a group per line.
fn truth(pf: &Pfam, dir: &str, input: &str, out: &str) {
    let lines = read_tsv(input);
    let pairs: Vec<(Vec<u8>, usize)> = lines.iter().map(|l| (encode(l[0].as_bytes()), l[1].parse().unwrap())).collect();
    let words = seed_words(dir);
    let names = accs(dir);
    let mut by_model: HashMap<usize, Vec<usize>> = HashMap::new();
    for (i, (_, m)) in pairs.iter().enumerate() {
        by_model.entry(*m).or_default().push(i);
    }
    let by_model: Vec<(usize, Vec<usize>)> = by_model.into_iter().collect();
    for (si, (name, shape, _)) in SHAPES.iter().enumerate() {
        let t0 = std::time::Instant::now();
        let parts: Vec<Vec<Row>> = by_model
            .par_iter()
            .map(|(m, idx)| {
                let info = info_percentiles(pf, *m, shape);
                let set = words.get(&names[*m]).map(|s| family_words(s, shape)).unwrap_or_default();
                let mut out = vec![];
                let mut buf = vec![];
                let mut memo = HashMap::new();
                for &i in idx {
                    let f = &pairs[i].0;
                    let before = out.len();
                    scan_model(pf, *m, f, si, -1.0, i as u32, &mut out, if info.is_empty() { None } else { Some(&info) }, &mut buf, false);
                    // Words below the recording bar are added; recorded ones are flagged.
                    let mut seen = HashSet::new();
                    for r in out[before..].iter_mut() {
                        if set.contains(&r.kmer) {
                            r.word = true;
                            seen.insert(r.j);
                        }
                    }
                    let span = shape.last().unwrap() + 1;
                    if f.len() >= span && !info.is_empty() {
                        for j in 0..=(f.len() - span) {
                            if seen.contains(&(j as u8)) {
                                continue;
                            }
                            if let Some(k) = kmer_of(f, j, shape) {
                                if set.contains(&k) {
                                    if let Some(r) = word_row(pf, *m, f, j, shape, k, i as u32, -1.0, &info, &mut buf, &mut memo) {
                                        out.push(r);
                                    }
                                }
                            }
                        }
                    }
                }
                out
            })
            .collect();
        let mut rows: Vec<Row> = parts.into_iter().flatten().collect();
        write_rows(&format!("{out}/{name}.rows.gz"), &mut rows, pairs.len());
        eprintln!("truth {name}: {} rows, {:.0}s", rows.len(), t0.elapsed().as_secs_f32());
    }
}
/// Input: read, orf per line. A group for each (read, model) with a seeded diagonal >= bar.
fn background(pf: &Pfam, dir: &str, input: &str, bar: f32, out: &str) {
    let lines = read_tsv(input);
    let orfs: Vec<(u32, Vec<u8>)> = lines.iter().map(|l| (l[0].parse().unwrap(), encode(l[1].as_bytes()))).collect();
    let reads = orfs.iter().map(|o| o.0).max().unwrap() as usize + 1;
    let words = seed_words(dir);
    let names = accs(dir);
    for (si, (name, shape, _)) in SHAPES.iter().enumerate() {
        let t0 = std::time::Instant::now();
        let span = shape.last().unwrap() + 1;
        // per-model: rows and distinct seed words (the words method's index size)
        let parts: Vec<(Vec<Row>, u32, f64)> = (0..pf.models())
            .into_par_iter()
            .map(|m| {
                let mut rows = vec![];
                let mut buf = vec![];
                let info = info_percentiles(pf, m, shape);
                let set = words.get(&names[m]).map(|s| family_words(s, shape)).unwrap_or_default();
                let mut memo = HashMap::new();
                for (read, f) in &orfs {
                    let g = (*read as u64 * 40000 + m as u64) as u32;
                    let before = rows.len();
                    scan_model(pf, m, f, si, bar, g, &mut rows, if info.is_empty() { None } else { Some(&info) }, &mut buf, true);
                    let mut seen = HashSet::new();
                    for r in rows[before..].iter_mut() {
                        if set.contains(&r.kmer) {
                            r.word = true;
                            seen.insert(r.j);
                        }
                    }
                    if f.len() >= span && !info.is_empty() && !set.is_empty() {
                        for j in 0..=(f.len() - span) {
                            if seen.contains(&(j as u8)) {
                                continue;
                            }
                            if let Some(k) = kmer_of(f, j, shape) {
                                if set.contains(&k) {
                                    if let Some(r) = word_row(pf, m, f, j, shape, k, g, bar, &info, &mut buf, &mut memo) {
                                        rows.push(r);
                                    }
                                }
                            }
                        }
                    }
                }
                let mut qsum = 0.0;
                for k in &set {
                    let mut x = *k;
                    let mut p = 1.0f64;
                    for _ in 0..shape.len() {
                        p *= pf.nbg[(x % 20) as usize] as f64;
                        x /= 20;
                    }
                    qsum += p;
                }
                (rows, set.len() as u32, qsum)
            })
            .collect();
        let mut wb = Vec::new();
        for (_, n, q) in &parts {
            wb.extend_from_slice(&n.to_le_bytes());
            wb.extend_from_slice(&(*q as f32).to_le_bytes());
        }
        fs::write(format!("{out}/{name}.words.bin"), wb).unwrap();
        let mut rows: Vec<Row> = parts.into_iter().flat_map(|p| p.0).collect();
        // renumber groups; the table of (read, model) is written beside the rows
        let mut keys: Vec<u32> = rows.iter().map(|r| r.group).collect();
        keys.sort_unstable();
        keys.dedup();
        let map: HashMap<u32, u32> = keys.iter().enumerate().map(|(i, k)| (*k, i as u32)).collect();
        rows.iter_mut().for_each(|r| r.group = map[&r.group]);
        write_rows(&format!("{out}/{name}.bg.gz"), &mut rows, keys.len());
        let mut b = Vec::new();
        for k in &keys {
            b.extend_from_slice(&((k / 40000) as u16).to_le_bytes());
            b.extend_from_slice(&((k % 40000) as u16).to_le_bytes());
        }
        write_gz(&format!("{out}/{name}.bggroups.gz"), &b);
        eprintln!("bg {name}: {} rows, {} (read, model), {reads} reads, {:.0}s", rows.len(), keys.len(), t0.elapsed().as_secs_f32());
    }
}
/// Input: id, seq per line (ids may repeat: one read's ORFs). Output: id, best score, its model, its second-best other model.
fn best(pf: &Pfam, input: &str, out: &str) {
    let lines = read_tsv(input);
    let seqs: Vec<(String, Vec<u8>)> = lines.iter().map(|l| (l[0].clone(), encode(l[1].as_bytes()))).collect();
    let res: Vec<(f32, usize)> = seqs
        .par_iter()
        .map(|(_, f)| {
            let mut buf = vec![];
            let (mut bs, mut bm) = (0.0f32, 0usize);
            let mut unused = vec![];
            for m in 0..pf.models() {
                let s = scan_model(pf, m, f, 0, f32::MAX, 0, &mut unused, None, &mut buf, false);
                if s > bs {
                    bs = s;
                    bm = m;
                }
            }
            (bs, bm)
        })
        .collect();
    let mut o = fs::File::create(out).unwrap();
    for ((id, _), (s, m)) in seqs.iter().zip(res) {
        writeln!(o, "{id}\t{s:.2}\t{m}").unwrap();
    }
}
fn main() {
    let a: Vec<String> = std::env::args().collect();
    let pf = Pfam::load(&a[2]);
    match a[1].as_str() {
        "tables" => tables(&pf, &a[3]),
        "truth" => truth(&pf, &a[2], &a[3], &a[4]),
        "bg" => background(&pf, &a[2], &a[3], a[4].parse().unwrap(), &a[5]),
        "best" => best(&pf, &a[3], &a[4]),
        _ => panic!("tables | truth | bg | best"),
    }
}
