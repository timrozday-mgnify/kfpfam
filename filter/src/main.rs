//! The first real index: every k-mer scoring at least T against a column, then an ungapped diagonal
//! score from each hit, then a gapped score on what passes.
//!
//! filter build  <prep> <k> <T> <out.idx>
//! filter search <prep> <out.idx> <orfs.faa> <bar> <out.tsv> [seeds.tsv]   query, model, diagonal, gapped; seeds: query, position, k-mer, model, diagonal, for each hit at the bar
//! filter window <prep> <out.idx> <orfs.faa> <bar> <out prefix>   speed plan step 1: window widths and quarter-bit scores
//! filter tile <prep> <out.idx> <out.tidx>   the index in tiles of models, by k-mer within a tile
//! filter join <prep> <out.tidx> <bar> <gate bits> <block ORFs> <orfs.faa> <out.tsv>   speed plan step 2: postings per block, joined from the model's side
//! filter sketch <prep> <orfs.faa> <k> <T> <scaled> <fwd|rev> <out.tsv>   per model, what a whole sample's k-mer set holds
//! <prep> holds profile.f32 (columns x 21 log-odds bits) and offs.u32, from measure/prep.py.
use rayon::prelude::*;
use std::{
    fs,
    io::{Read, Write},
    sync::atomic::{AtomicU64, Ordering::Relaxed},
    time::Instant,
};

const A: &[u8] = b"ACDEFGHIKLMNPQRSTVWY";
const GAP_OPEN: f32 = 6.0;
const GAP_EXTEND: f32 = 1.0;

struct Pfam {
    lo: Vec<f32>,
    offs: Vec<u32>,
}
impl Pfam {
    fn load(dir: &str) -> Pfam {
        // read straight into the floats, so the profile is never held twice (native little-endian)
        let path = format!("{dir}/profile.f32");
        let mut lo = vec![0f32; fs::metadata(&path).unwrap().len() as usize / 4];
        fs::File::open(&path).unwrap().read_exact(bytes_of_mut(&mut lo)).unwrap();
        let b = fs::read(format!("{dir}/offs.u32")).unwrap();
        let offs = b.chunks_exact(4).map(|c| u32::from_le_bytes(c.try_into().unwrap())).collect();
        Pfam { lo, offs }
    }
    fn models(&self) -> usize {
        self.offs.len() - 1
    }
}

/// Calls `visit` with each k-mer scoring >= t in the k columns from `at`, and its score.
/// `keep(i, prefix)` can drop a prefix of length i, and everything under it.
fn enumerate(lo: &[f32], at: usize, k: usize, t: f32, keep: &impl Fn(usize, usize) -> bool, visit: &mut impl FnMut(usize, f32)) {
    let mut rest = [0.0f32; 16];
    for i in (0..k).rev() {
        let top = (0..20).map(|r| lo[(at + i) * 21 + r]).fold(f32::MIN, f32::max);
        rest[i] = rest[i + 1] + top;
    }
    #[allow(clippy::too_many_arguments)]
    fn go(lo: &[f32], at: usize, i: usize, k: usize, code: usize, sum: f32, t: f32, rest: &[f32], keep: &impl Fn(usize, usize) -> bool, visit: &mut impl FnMut(usize, f32)) {
        if i == k {
            visit(code, sum);
            return;
        }
        for r in 0..20 {
            let s = sum + lo[(at + i) * 21 + r];
            let next = code * 20 + r;
            if s + rest[i + 1] >= t && keep(i + 1, next) {
                go(lo, at, i + 1, k, next, s, t, rest, keep, visit);
            }
        }
    }
    go(lo, at, 0, k, 0, 0.0, t, &rest, keep, visit);
}

fn build(pf: &Pfam, k: usize, t_all: f32, out: &str) {
    // FILTER_T_FILE: one T per model, one per line, in place of the global T
    let per_model: Vec<f32> = match std::env::var("FILTER_T_FILE") {
        Ok(p) => fs::read_to_string(p).unwrap().lines().map(|l| l.trim().parse().unwrap()).collect(),
        Err(_) => vec![t_all; pf.models()],
    };
    assert_eq!(per_model.len(), pf.models());
    let t0 = Instant::now();
    let space = 20usize.pow(k as u32);
    let mut starts = vec![0u64; space + 1];
    let windows = |m: usize| {
        let (a, b) = (pf.offs[m] as usize, pf.offs[m + 1] as usize);
        (a, if b - a >= k { b - a - k + 1 } else { 0 })
    };
    for m in 0..pf.models() {
        let (a, n) = windows(m);
        for c in 0..n {
            enumerate(&pf.lo, a + c, k, per_model[m], &|_, _| true, &mut |code, _| starts[code + 1] += 1);
        }
    }
    for i in 0..space {
        starts[i + 1] += starts[i];
    }
    let total = starts[space] as usize;
    eprintln!("{total} entries, {:.1} GB, counted in {:.0}s", total as f64 * 4.0 / 1e9, t0.elapsed().as_secs_f32());
    let mut entries = vec![0u32; total];
    let mut cursor = starts.clone();
    for m in 0..pf.models() {
        let (a, n) = windows(m);
        for c in 0..n {
            assert!(m < 65536 && c < 65536);
            enumerate(&pf.lo, a + c, k, per_model[m], &|_, _| true, &mut |code, _| {
                entries[cursor[code] as usize] = (m as u32) << 16 | c as u32;
                cursor[code] += 1;
            });
        }
    }
    let mut f = std::io::BufWriter::new(fs::File::create(out).unwrap());
    f.write_all(&(k as u64).to_le_bytes()).unwrap();
    f.write_all(&(total as u64).to_le_bytes()).unwrap();
    f.write_all(bytes_of(&starts)).unwrap();
    f.write_all(bytes_of(&entries)).unwrap();
    eprintln!("built in {:.0}s", t0.elapsed().as_secs_f32());
}

fn bytes_of<T>(v: &[T]) -> &[u8] {
    // ponytail: native-endian dump of integers; an index is rebuilt rather than moved between machines
    unsafe { std::slice::from_raw_parts(v.as_ptr() as *const u8, std::mem::size_of_val(v)) }
}
fn bytes_of_mut<T>(v: &mut [T]) -> &mut [u8] {
    unsafe { std::slice::from_raw_parts_mut(v.as_mut_ptr() as *mut u8, std::mem::size_of_val(v)) }
}

struct Index {
    k: usize,
    starts: Vec<u64>,
    entries: Vec<u32>,
}
fn load_index(path: &str) -> Index {
    let mut f = fs::File::open(path).unwrap();
    let mut h = [0u64; 2];
    f.read_exact(bytes_of_mut(&mut h)).unwrap();
    let k = h[0] as usize;
    let mut starts = vec![0u64; 20usize.pow(k as u32) + 1];
    f.read_exact(bytes_of_mut(&mut starts)).unwrap();
    let mut entries = vec![0u32; h[1] as usize];
    f.read_exact(bytes_of_mut(&mut entries)).unwrap();
    Index { k, starts, entries }
}

fn read_orfs(path: &str) -> Vec<(String, Vec<u8>)> {
    let mut t = [20u8; 256];
    for (i, c) in A.iter().enumerate() {
        t[*c as usize] = i as u8;
    }
    let text = fs::read_to_string(path).unwrap();
    let mut out = vec![];
    let mut lines = text.lines();
    while let (Some(name), Some(seq)) = (lines.next(), lines.next()) {
        out.push((name[1..].split_whitespace().next().unwrap().to_string(), seq.bytes().map(|c| t[c as usize]).collect()));
    }
    out
}

/// ORFs packed into two buffers, as millions of small allocations cost more than the sequences.
struct Packed {
    names: String,
    name_at: Vec<usize>,
    seqs: Vec<u8>,
    seq_at: Vec<usize>,
}
impl Packed {
    /// The next n ORFs (or fewer at the end) of a one-line-per-sequence FASTA, so a sample is held one
    /// block at a time rather than whole (about 0.46 GB per million reads otherwise).
    fn next(lines: &mut impl Iterator<Item = std::io::Result<String>>, n: usize) -> Packed {
        let mut t = [20u8; 256];
        for (i, c) in A.iter().enumerate() {
            t[*c as usize] = i as u8;
        }
        let mut p = Packed { names: String::new(), name_at: vec![0], seqs: vec![], seq_at: vec![0] };
        while p.len() < n {
            let Some(line) = lines.next() else { break };
            let line = line.unwrap();
            if let Some(name) = line.strip_prefix('>') {
                p.names.push_str(name.split_whitespace().next().unwrap());
                p.name_at.push(p.names.len());
            } else {
                p.seqs.extend(line.bytes().map(|c| t[c as usize]));
                p.seq_at.push(p.seqs.len());
            }
        }
        assert_eq!(p.name_at.len(), p.seq_at.len());
        p
    }
    fn len(&self) -> usize {
        self.seq_at.len() - 1
    }
    fn name(&self, i: usize) -> &str {
        &self.names[self.name_at[i]..self.name_at[i + 1]]
    }
    fn seq(&self, i: usize) -> &[u8] {
        &self.seqs[self.seq_at[i]..self.seq_at[i + 1]]
    }
}

/// max without f32::max's NaN rule, which costs baseline x86 three instructions where maxss is one.
#[inline(always)]
fn fmax(a: f32, b: f32) -> f32 {
    if a > b { a } else { b }
}

/// Best local alignment with affine gaps, in bits. Whole matrix and scalar: about 16 core-us per pair.
fn gapped(lo: &[f32], a: usize, cols: usize, f: &[u8], h: &mut Vec<f32>, del: &mut Vec<f32>) -> f32 {
    h.clear();
    h.resize(cols + 1, 0.0);
    del.clear();
    del.resize(cols + 1, f32::NEG_INFINITY);
    let mut best = 0.0f32;
    for &r in f {
        let (mut diag, mut left, mut ins) = (0.0f32, 0.0f32, f32::NEG_INFINITY);
        for c in 1..=cols {
            let up = h[c];
            del[c] = (up - GAP_OPEN).max(del[c] - GAP_EXTEND);
            ins = (left - GAP_OPEN).max(ins - GAP_EXTEND);
            let v = (diag + lo[(a + c - 1) * 21 + r as usize]).max(del[c]).max(ins).max(0.0);
            diag = up;
            h[c] = v;
            left = v;
            best = best.max(v);
        }
    }
    best
}

/// gapped() restricted to a band of +-band columns around diagonal d (column - read position).
/// Cells outside the band count as unaligned. Same scores as gapped() when the band holds the alignment.
fn gapped_band(lo: &[f32], a: usize, cols: usize, f: &[u8], d: i64, band: i64, h: &mut Vec<f32>, del: &mut Vec<f32>) -> f32 {
    h.clear();
    h.resize(cols + 2, 0.0);
    del.clear();
    del.resize(cols + 2, f32::NEG_INFINITY);
    let mut best = 0.0f32;
    for (i, &r) in f.iter().enumerate() {
        // 1-based columns c hold model column c - 1; the band is model columns i + d +- band
        let c_lo = (i as i64 + d - band + 1).max(1);
        let c_hi = (i as i64 + d + band + 1).min(cols as i64);
        if c_lo > c_hi {
            continue;
        }
        let (c_lo, c_hi) = (c_lo as usize, c_hi as usize);
        // the cell entering the band on the right was not computed last row
        h[c_hi] = if c_hi as i64 == (i as i64 - 1) + d + band + 1 { h[c_hi] } else { 0.0 };
        if (c_hi as i64) > (i as i64 - 1) + d + band + 1 || i == 0 {
            del[c_hi] = f32::NEG_INFINITY;
        }
        let mut diag = h[c_lo - 1];
        if c_lo as i64 - 1 < (i as i64 - 1) + d - band + 1 {
            diag = 0.0; // left of last row's band
        }
        // a gap opens from m, the cell before `ins`: opening from a gap never beats extending it, and the
        // row's dependency chain is then ins alone
        let (mut m, mut ins) = (0.0f32, f32::NEG_INFINITY);
        let n = c_hi - c_lo + 1;
        let (h, del, lo) = (&mut h[c_lo..][..n], &mut del[c_lo..][..n], &lo[(a + c_lo - 1) * 21..][..n * 21]);
        for c in 0..n {
            let up = h[c];
            let dl = fmax(up - GAP_OPEN, del[c] - GAP_EXTEND);
            del[c] = dl;
            ins = fmax(m - GAP_OPEN, ins - GAP_EXTEND);
            m = fmax(fmax(diag + lo[c * 21 + r as usize], dl), 0.0);
            diag = up;
            let v = fmax(m, ins);
            h[c] = v;
            best = fmax(best, v);
        }
    }
    best
}

/// What the gapped bound needs of a SIMD register: 16 signed bytes.
#[cfg(any(target_arch = "aarch64", target_arch = "x86_64"))]
trait Bytes16: Copy {
    unsafe fn load(p: *const i8) -> Self;
    unsafe fn splat(v: i8) -> Self;
    unsafe fn adds(self, o: Self) -> Self;
    unsafe fn subs(self, o: Self) -> Self;
    unsafe fn max(self, o: Self) -> Self;
    unsafe fn min(self, o: Self) -> Self;
    /// Bytes N..N + 16 of self followed by o.
    unsafe fn ext<const N: i32>(self, o: Self) -> Self;
    unsafe fn hmax(self) -> i8;
}
#[cfg(target_arch = "aarch64")]
impl Bytes16 for std::arch::aarch64::int8x16_t {
    #[inline(always)]
    unsafe fn load(p: *const i8) -> Self { std::arch::aarch64::vld1q_s8(p) }
    #[inline(always)]
    unsafe fn splat(v: i8) -> Self { std::arch::aarch64::vdupq_n_s8(v) }
    #[inline(always)]
    unsafe fn adds(self, o: Self) -> Self { std::arch::aarch64::vqaddq_s8(self, o) }
    #[inline(always)]
    unsafe fn subs(self, o: Self) -> Self { std::arch::aarch64::vqsubq_s8(self, o) }
    #[inline(always)]
    unsafe fn max(self, o: Self) -> Self { std::arch::aarch64::vmaxq_s8(self, o) }
    #[inline(always)]
    unsafe fn min(self, o: Self) -> Self { std::arch::aarch64::vminq_s8(self, o) }
    #[inline(always)]
    unsafe fn ext<const N: i32>(self, o: Self) -> Self { std::arch::aarch64::vextq_s8::<N>(self, o) }
    #[inline(always)]
    unsafe fn hmax(self) -> i8 { std::arch::aarch64::vmaxvq_s8(self) }
}
#[cfg(target_arch = "x86_64")]
impl Bytes16 for std::arch::x86_64::__m128i {
    #[inline(always)]
    unsafe fn load(p: *const i8) -> Self { std::arch::x86_64::_mm_loadu_si128(p as *const _) }
    #[inline(always)]
    unsafe fn splat(v: i8) -> Self { std::arch::x86_64::_mm_set1_epi8(v) }
    #[inline(always)]
    unsafe fn adds(self, o: Self) -> Self { std::arch::x86_64::_mm_adds_epi8(self, o) }
    #[inline(always)]
    unsafe fn subs(self, o: Self) -> Self { std::arch::x86_64::_mm_subs_epi8(self, o) }
    #[inline(always)]
    unsafe fn max(self, o: Self) -> Self { std::arch::x86_64::_mm_max_epi8(self, o) }
    #[inline(always)]
    unsafe fn min(self, o: Self) -> Self { std::arch::x86_64::_mm_min_epi8(self, o) }
    #[inline(always)]
    unsafe fn ext<const N: i32>(self, o: Self) -> Self { std::arch::x86_64::_mm_alignr_epi8::<N>(o, self) }
    #[inline(always)]
    unsafe fn hmax(self) -> i8 {
        let mut b = [0i8; 16];
        std::arch::x86_64::_mm_storeu_si128(b.as_mut_ptr() as *mut _, self);
        b.into_iter().max().unwrap()
    }
}

/// An upper bound on gapped_band's score, in bytes, one read residue per step with the band in one
/// register, so a pair that cannot reach `min` is dropped without the f32 matrix. Scores are rounded up,
/// gap costs down and the band's 17th column is merged into the 16th, so it never drops a pair
/// gapped_band would keep.
struct Bound {
    prof: Vec<i8>, // per model, 21 residues x (PAD + columns + PAD), residue-major
    at: Vec<usize>,
    scale: f32,
    pass: i8, // the lowest bound that may reach min, offset by -128
}
const PAD: usize = 24; // columns either side of a model: the band may hang off both ends
const BAND: i64 = 8;
impl Bound {
    /// None when there is nothing to bound (no min), or no SIMD for it here.
    fn new(pf: &Pfam, min: f32) -> Option<Bound> {
        if !(16.0..=250.0).contains(&min) || std::env::var("FILTER_SIMD").as_deref() == Ok("scalar") {
            return None;
        }
        #[cfg(target_arch = "x86_64")]
        if !is_x86_feature_detected!("avx2") {
            return None;
        }
        #[cfg(not(any(target_arch = "x86_64", target_arch = "aarch64")))]
        return None;
        let scale = (254.0 / min).floor(); // at most 15 here, so 8 gap extensions fit in a byte
        let mut at = vec![0usize];
        for m in 0..pf.models() {
            at.push(at[m] + 21 * ((pf.offs[m + 1] - pf.offs[m]) as usize + 2 * PAD));
        }
        let mut prof = vec![i8::MIN; *at.last().unwrap()];
        let mut parts = Vec::with_capacity(pf.models());
        let mut rest: &mut [i8] = &mut prof;
        for m in 0..pf.models() {
            let (head, tail) = std::mem::take(&mut rest).split_at_mut(at[m + 1] - at[m]);
            parts.push(head);
            rest = tail;
        }
        parts.into_par_iter().enumerate().for_each(|(m, out)| {
            let (a, stride) = (pf.offs[m] as usize, out.len() / 21);
            for c in 0..stride - 2 * PAD {
                for r in 0..21 {
                    out[r * stride + PAD + c] = (pf.lo[(a + c) * 21 + r] * scale).ceil().clamp(-128.0, 127.0) as i8;
                }
            }
        });
        Some(Bound { prof, at, scale, pass: ((min * scale - 0.01).floor() - 128.0) as i8 })
    }
    /// Whether each pair (model, read, diagonal as column - read position) may reach min.
    fn may_reach(&self, pairs: [(usize, &[u8], i64); PAIRS]) -> [bool; PAIRS] {
        let jobs = pairs.map(|(m, f, d)| (&self.prof[self.at[m]..self.at[m + 1]], f, d));
        #[cfg(target_arch = "aarch64")]
        return unsafe { bound::<std::arch::aarch64::int8x16_t>(jobs, self.scale) }.map(|b| b >= self.pass);
        #[cfg(target_arch = "x86_64")]
        {
            // AVX2 for its three-operand encoding; the registers are 16 bytes either way
            #[target_feature(enable = "avx2")]
            unsafe fn run(jobs: [Job; PAIRS], scale: f32) -> [i8; PAIRS] {
                bound::<std::arch::x86_64::__m128i>(jobs, scale)
            }
            return unsafe { run(jobs, self.scale) }.map(|b| b >= self.pass);
        }
        #[allow(unreachable_code)]
        [true; PAIRS]
    }
}
const PAIRS: usize = 4; // scored side by side, as each pair's rows are one long dependency chain
type Job<'a> = (&'a [i8], &'a [u8], i64); // a model's rows of Bound::prof, read, diagonal

/// The band as lanes: lane l is column i + d - 8 + l at read position i, so last row's cell on the same
/// diagonal is the same lane and the cell above is the next lane. Zero is -128, so saturation is the local
/// score's floor. Gaps along the model are a prefix maximum, taken in doubling steps.
#[cfg(any(target_arch = "aarch64", target_arch = "x86_64"))]
#[inline(always)]
unsafe fn bound<V: Bytes16>(jobs: [Job; PAIRS], scale: f32) -> [i8; PAIRS] {
    let neg = V::splat(i8::MIN);
    let ext = (GAP_EXTEND * scale) as i8;
    let open = V::splat((GAP_OPEN * scale) as i8);
    let (e1, e2, e4, e8) = (V::splat(ext), V::splat(2 * ext), V::splat(4 * ext), V::splat(8 * ext));
    let mut last = [i8::MIN; 16];
    last[15] = i8::MAX;
    let last = V::load(last.as_ptr());
    let (mut h, mut del, mut best) = ([neg; PAIRS], [neg; PAIRS], [neg; PAIRS]);
    // each pair runs from the first read position whose band touches its model to the last
    let rows = jobs.map(|(prof, f, d)| {
        let cols = (prof.len() / 21 - 2 * PAD) as i64;
        ((-d - BAND).max(0), (f.len() as i64).min(cols - d + BAND))
    });
    let steps = rows.iter().map(|&(from, to)| to - from).max().unwrap().max(0);
    for step in 0..steps {
        for p in 0..PAIRS {
            let (prof, f, d) = jobs[p];
            let i = rows[p].0 + step;
            // a pair that is done reads padding, the first bytes of its profile
            let at = if i < rows[p].1 {
                *f.get_unchecked(i as usize) as usize * (prof.len() / 21) + (PAD as i64 + i + d - BAND) as usize
            } else {
                0
            };
            debug_assert!(at + 17 <= prof.len());
            let s = V::load(prof.as_ptr().add(at)).max(V::load(prof.as_ptr().add(at + 1)).min(last));
            del[p] = h[p].ext::<1>(neg).subs(open).max(del[p].ext::<1>(neg).subs(e1));
            let m = h[p].adds(s).max(del[p]);
            let mut x = neg.ext::<15>(m.subs(open));
            x = x.max(neg.ext::<15>(x).subs(e1));
            x = x.max(neg.ext::<14>(x).subs(e2));
            x = x.max(neg.ext::<12>(x).subs(e4));
            x = x.max(neg.ext::<8>(x).subs(e8));
            h[p] = m.max(x);
            best[p] = best[p].max(h[p]);
        }
    }
    best.map(|b| b.hmax())
}

fn search(pf: &Pfam, idx: &Index, path: &str, bar: f32, out: &str, seeds: Option<&String>) {
    let orfs = read_orfs(path);
    let (hits, passed, align_ns) = (AtomicU64::new(0), AtomicU64::new(0), AtomicU64::new(0));
    let t0 = Instant::now();
    let chunks: Vec<(String, String)> = orfs
        .par_chunks(64)
        .map(|chunk| {
            let (mut text, mut seed_text, mut h, mut del, mut found) = (String::new(), String::new(), vec![], vec![], Vec::<(u32, f32)>::new());
            let (mut n_hits, mut n_passed, mut ns) = (0u64, 0u64, 0u64);
            for (name, f) in chunk {
                found.clear();
                if f.len() < idx.k {
                    continue;
                }
                for j in 0..=(f.len() - idx.k) {
                    let mut code = 0usize;
                    let mut ok = true;
                    for &c in &f[j..j + idx.k] {
                        ok &= c < 20;
                        code = code * 20 + c as usize;
                    }
                    if !ok {
                        continue;
                    }
                    let list = &idx.entries[idx.starts[code] as usize..idx.starts[code + 1] as usize];
                    n_hits += list.len() as u64;
                    for &e in list {
                        let (m, col) = ((e >> 16) as usize, (e & 0xffff) as usize);
                        let a = pf.offs[m] as usize;
                        let cols = pf.offs[m + 1] as usize - a;
                        let back = j.min(col);
                        let len = (f.len() - j + back).min(cols - col + back);
                        let mut at = (a + col - back) * 21;
                        let (mut cur, mut best) = (0.0f32, 0.0f32);
                        for &r in &f[j - back..j - back + len] {
                            cur = (cur + pf.lo[at + r as usize]).max(0.0);
                            best = best.max(cur);
                            at += 21;
                        }
                        if best >= bar {
                            if seeds.is_some() {
                                seed_text += &format!("{name}\t{j}\t{code}\t{m}\t{best:.2}\n");
                            }
                            match found.iter_mut().find(|x| x.0 == m as u32) {
                                Some(x) => x.1 = x.1.max(best),
                                None => found.push((m as u32, best)),
                            }
                        }
                    }
                }
                let t = Instant::now();
                for &(m, d) in &found {
                    let a = pf.offs[m as usize] as usize;
                    let cols = pf.offs[m as usize + 1] as usize - a;
                    let g = gapped(&pf.lo, a, cols, f, &mut h, &mut del);
                    text += &format!("{name}\t{m}\t{d:.2}\t{g:.2}\n");
                }
                ns += t.elapsed().as_nanos() as u64;
                n_passed += found.len() as u64;
            }
            hits.fetch_add(n_hits, Relaxed);
            passed.fetch_add(n_passed, Relaxed);
            align_ns.fetch_add(ns, Relaxed);
            (text, seed_text)
        })
        .collect();
    let wall = t0.elapsed().as_secs_f64();
    let threads = rayon::current_num_threads() as f64;
    let align = align_ns.load(Relaxed) as f64 / 1e9;
    let (pairs, seed_rows): (Vec<String>, Vec<String>) = chunks.into_iter().unzip();
    fs::write(out, pairs.concat()).unwrap();
    if let Some(path) = seeds {
        fs::write(path, seed_rows.concat()).unwrap();
    }
    eprintln!(
        "{} ORFs, {:.1}s wall on {threads} threads, {:.0} core-s; alignment {:.1} core-s; {} hits, {} passed",
        orfs.len(), wall, wall * threads, align, hits.load(Relaxed), passed.load(Relaxed)
    );
}

const BINS: usize = 48; // half-bit bins above T
const SET_BITS: u32 = 34; // 2 GB; about 1 false bit in 40 at 2 x 10^8 residues and 2 levels
const MIN_PREFIX: usize = 6; // shorter prefixes are saturated in any sample, so not worth a lookup

/// Where a k-mer of length `len` falls in the sample's set.
fn slot(code: usize, len: usize) -> usize {
    let h = (code as u64 ^ (len as u64) << 58).wrapping_mul(0x9E37_79B9_7F4A_7C15);
    (h >> (64 - SET_BITS)) as usize
}

struct Bits(Vec<AtomicU64>);
impl Bits {
    fn new(n: usize) -> Bits {
        Bits((0..n / 64 + 1).map(|_| AtomicU64::new(0)).collect())
    }
    fn set(&self, i: usize) {
        let (w, b) = (&self.0[i / 64], 1u64 << (i % 64));
        if w.load(Relaxed) & b == 0 {
            w.fetch_or(b, Relaxed);
        }
    }
    fn get(&self, i: usize) -> bool {
        self.0[i / 64].load(Relaxed) >> (i % 64) & 1 == 1
    }
    fn count(&self) -> u64 {
        self.0.iter().map(|w| w.load(Relaxed).count_ones() as u64).sum()
    }
}

/// FracMinHash's rule: keep a k-mer when its hash falls in the lowest 1/scaled of the range.
fn in_sketch(code: usize, scaled: u64) -> bool {
    scaled == 1 || (code as u64).wrapping_mul(0xC2B2_AE3D_27D4_EB4F) >> 32 < (1u64 << 32) / scaled
}

/// One sample's set of k-mers, with no read structure, and what each model finds in it.
/// reverse: use reversed sequences, the null.
#[allow(clippy::too_many_arguments)]
fn sketch(pf: &Pfam, path: &str, k: usize, t: f32, scaled: u64, reverse: bool, out: &str) {
    let t0 = Instant::now();
    let orfs = read_orfs(path);
    let parse = t0.elapsed().as_secs_f64();
    let threads = rayon::current_num_threads() as f64;
    let shrink: Vec<usize> = (0..=k.saturating_sub(MIN_PREFIX)).scan(1usize, |p, _| { let v = *p; *p *= 20; Some(v) }).collect();
    let space = 20usize.pow(k as u32);
    let set = Bits::new(1 << SET_BITS); // k-mers (hash-filtered at scaled > 1) and their prefixes, hashed
    let add = |s: &[u8]| {
        let (mut code, mut run) = (0usize, 0usize);
        for &c in s {
            if c >= 20 {
                run = 0;
                continue;
            }
            code = (code * 20 + c as usize) % space;
            run += 1;
            if run >= k {
                if in_sketch(code, scaled) {
                    set.set(slot(code, k));
                }
                for (d, p) in shrink.iter().enumerate().skip(1) {
                    set.set(slot(code / p, k - d));
                }
            }
        }
    };
    let t1 = Instant::now();
    orfs.par_iter().for_each(|(_, f)| {
        if reverse {
            add(&f.iter().rev().copied().collect::<Vec<_>>());
        } else {
            add(f);
        }
    });
    let sample = t1.elapsed().as_secs_f64();
    let residues: u64 = orfs.iter().map(|o| o.1.len() as u64).sum();
    let filled = set.count();

    let t3 = Instant::now();
    let keep = |i: usize, code: usize| i == k || i < MIN_PREFIX || set.get(slot(code, i));
    let visited = AtomicU64::new(0);
    let rows: Vec<String> = (0..pf.models())
        .into_par_iter()
        .map(|m| {
            let (a, b) = (pf.offs[m] as usize, pf.offs[m + 1] as usize);
            let n = if b - a >= k { b - a - k + 1 } else { 0 };
            let (mut windows, mut hits, mut best, mut leaves) = ([0u32; BINS], [0u32; BINS], f32::MIN, 0u64);
            let mut chain = vec![0.0f32; n + 1]; // best sum over non-overlapping windows
            for c in 0..n {
                let mut top = f32::MIN;
                enumerate(&pf.lo, a + c, k, t, &keep, &mut |code, s| {
                    leaves += 1;
                    if set.get(slot(code, k)) {
                        top = top.max(s);
                        hits[(((s - t) * 2.0) as usize).min(BINS - 1)] += 1;
                    }
                });
                chain[c + 1] = chain[c];
                if top > f32::MIN {
                    windows[(((top - t) * 2.0) as usize).min(BINS - 1)] += 1;
                    best = best.max(top);
                    chain[c + 1] = chain[c].max(chain[(c + 1).saturating_sub(k)] + top);
                }
            }
            visited.fetch_add(leaves, Relaxed);
            let list = |v: &[u32]| v.iter().map(u32::to_string).collect::<Vec<_>>().join(",");
            format!("{m}\t{best:.2}\t{:.2}\t{}\t{}\n", chain[n], list(&windows), list(&hits))
        })
        .collect();
    let score = t3.elapsed().as_secs_f64();
    fs::write(out, rows.concat()).unwrap();
    eprintln!(
        "k={k} T={t} scaled={scaled} rev={reverse} orfs={} residues={residues} set_bits={filled} filled={:.4} leaves={} \
         parse_s={parse:.1} sample_core_s={:.1} score_core_s={:.1} ns_per_residue={:.1} threads={threads} out={out}",
        orfs.len(),
        filled as f64 / (1u64 << SET_BITS) as f64,
        visited.load(Relaxed),
        sample * threads,
        score * threads,
        sample * threads * 1e9 / residues as f64,
    );
}

const WIDTHS: [usize; 5] = [12, 16, 20, 24, 32];
const CENTRED: [usize; 5] = [24, 32, 40, 48, 64]; // FILTER_CENTRED=1: fixed windows centred on the k-mer, not clamped into the ORF
const WBINS: usize = 256; // quarter bits from 0

/// Speed plan step 1: for every hit, the best local score in a window of each width around its k-mer, and
/// the whole diagonal in float and in quarter-bit integers. Writes a histogram of window scores over all
/// hits (width, bin, count) and, per (query, model) with a diagonal near the bar, the best window score of
/// each width among hits whose own diagonal reaches the bar.
fn window(pf: &Pfam, idx: &Index, path: &str, bar: f32, out: &str) {
    let orfs = read_orfs(path);
    let k = idx.k;
    let q: Vec<i32> = pf.lo.iter().map(|&v| (v * 4.0).round().max(-400.0) as i32).collect();
    let centred = std::env::var("FILTER_CENTRED").is_ok();
    let widths = if centred { CENTRED } else { WIDTHS };
    let t0 = Instant::now();
    let parts: Vec<(Vec<u64>, String)> = orfs
        .par_chunks(64)
        .map(|chunk| {
            let mut hist = vec![0u64; WIDTHS.len() * WBINS];
            let mut text = String::new();
            // per model in this ORF: best diagonal, best quantized diagonal, best window per width at the bar
            let mut found: Vec<(u32, f32, i32, [f32; 5])> = vec![];
            for (name, f) in chunk {
                found.clear();
                if f.len() < k {
                    continue;
                }
                for j in 0..=(f.len() - k) {
                    let mut code = 0usize;
                    let mut ok = true;
                    for &c in &f[j..j + k] {
                        ok &= c < 20;
                        code = code * 20 + c as usize;
                    }
                    if !ok {
                        continue;
                    }
                    for &e in &idx.entries[idx.starts[code] as usize..idx.starts[code + 1] as usize] {
                        let (m, col) = ((e >> 16) as usize, (e & 0xffff) as usize);
                        let a = pf.offs[m] as usize;
                        let cols = pf.offs[m + 1] as usize - a;
                        let back = j.min(col);
                        let len = (f.len() - j + back).min(cols - col + back);
                        let (lo_i, hi_i) = (j - back, j - back + len); // read positions on the model
                        let base = (a + col) as isize - j as isize; // column of read position i is base + i
                        let (mut cur, mut best, mut cq, mut bq) = (0.0f32, 0.0f32, 0i32, 0i32);
                        for i in lo_i..hi_i {
                            let at = (base + i as isize) as usize * 21 + f[i] as usize;
                            cur = (cur + pf.lo[at]).max(0.0);
                            best = best.max(cur);
                            cq = (cq + q[at]).max(0);
                            bq = bq.max(cq);
                        }
                        let mut wbest = [0.0f32; 5];
                        for (w, &width) in widths.iter().enumerate() {
                            let (s0, s1) = if centred {
                                let half = (width - k) / 2;
                                (j.saturating_sub(half).max(lo_i), (j + k + half).min(hi_i))
                            } else {
                                let start = if f.len() >= width { (j as isize - ((width - k) / 2) as isize).clamp(0, (f.len() - width) as isize) as usize } else { 0 };
                                (start.max(lo_i), (start + width).min(f.len()).min(hi_i))
                            };
                            let (mut c, mut b) = (0.0f32, 0.0f32);
                            for i in s0..s1 {
                                c = (c + pf.lo[(base + i as isize) as usize * 21 + f[i] as usize]).max(0.0);
                                b = b.max(c);
                            }
                            wbest[w] = b;
                            hist[w * WBINS + ((b * 4.0) as usize).min(WBINS - 1)] += 1;
                        }
                        if best >= bar - 2.0 || bq as f32 >= (bar - 2.0) * 4.0 {
                            let slot = match found.iter().position(|x| x.0 == m as u32) {
                                Some(p) => p,
                                None => {
                                    found.push((m as u32, 0.0, 0, [0.0; 5]));
                                    found.len() - 1
                                }
                            };
                            let x = &mut found[slot];
                            x.1 = x.1.max(best);
                            x.2 = x.2.max(bq);
                            if best >= bar {
                                for w in 0..5 {
                                    x.3[w] = x.3[w].max(wbest[w]);
                                }
                            }
                        }
                    }
                }
                for (m, d, dq, wb) in &found {
                    text += &format!("{name}\t{m}\t{d:.2}\t{:.2}\t{}\n", *dq as f32 / 4.0,
                        wb.iter().map(|v| format!("{v:.2}")).collect::<Vec<_>>().join("\t"));
                }
            }
            (hist, text)
        })
        .collect();
    let mut hist = vec![0u64; WIDTHS.len() * WBINS];
    let mut text = String::new();
    for (h, t) in parts {
        for (a, b) in hist.iter_mut().zip(h) {
            *a += b;
        }
        text += &t;
    }
    fs::write(format!("{out}.pairs.tsv"), text).unwrap();
    let rows: String = widths.iter().enumerate()
        .flat_map(|(w, width)| (0..WBINS).map(move |b| (w, *width, b)))
        .map(|(w, width, b)| format!("{width}\t{b}\t{}\n", hist[w * WBINS + b]))
        .collect();
    fs::write(format!("{out}.hist.tsv"), rows).unwrap();
    eprintln!("{} ORFs in {:.1}s", orfs.len(), t0.elapsed().as_secs_f64());
}

// ---- Speed plan step 2: join from the model's side ----

const WIN: usize = 48; // residues per posting, the k-mer at HALF..HALF + k
const HALF: usize = 22;
const STRIDE: usize = 32; // bytes per profile column: 21 residues, then padding

/// The index in tiles: each tile is a run of whole models whose padded int8 profile has at most
/// TILE_ROWS rows (2 MB), so it fits in L2 and an entry fits in a u16. Within a tile, entries are grouped by k-mer, so each postings list is read once
/// per tile and scored against every entry of the tile that holds its k-mer. An entry is the row of the
/// tile's padded profile: model base + column.
const TILE_ROWS: usize = u16::MAX as usize;
const SURV_COMPACT: usize = 1 << 22; // survivors held before sorting out repeats: 32 MB per thread
const SLACK: f32 = 3.0; // bits: quarter-bit rounding moved a whole diagonal by at most 2.1 bits in step 1
struct Tile {
    m0: usize,
    m1: usize,
    starts: Vec<u32>, // per k-mer, into rows
    rows: Vec<u16>,
}
struct TileIndex {
    k: usize,
    tiles: Vec<Tile>,
}

/// Padded rows of each model in a tile: model m's column c is row base[m - m0] + c.
fn tile_bases(pf: &Pfam, m0: usize, m1: usize) -> Vec<usize> {
    let mut base = vec![0usize; m1 - m0 + 1];
    for m in m0..m1 {
        base[m - m0 + 1] = base[m - m0] + (pf.offs[m + 1] - pf.offs[m]) as usize + WIN;
    }
    base
}

fn tile(pf: &Pfam, idx: &Index, out: &str) {
    let mut bounds = vec![0usize];
    let mut rows_in = 0usize;
    for m in 0..pf.models() {
        let rows_m = (pf.offs[m + 1] - pf.offs[m]) as usize + WIN;
        assert!(rows_m <= TILE_ROWS);
        if rows_in + rows_m > TILE_ROWS {
            bounds.push(m);
            rows_in = 0;
        }
        rows_in += rows_m;
    }
    if *bounds.last().unwrap() != pf.models() {
        bounds.push(pf.models());
    }
    let space = idx.starts.len() - 1;
    let tiles: Vec<Tile> = bounds.windows(2).collect::<Vec<_>>().into_par_iter().map(|b| {
        let (m0, m1) = (b[0], b[1]);
        let base = tile_bases(pf, m0, m1);
        let (mut starts, mut rows) = (vec![0u32; space + 1], vec![]);
        for code in 0..space {
            // a k-mer's entries are in model order, so the tile's are one run
            let list = &idx.entries[idx.starts[code] as usize..idx.starts[code + 1] as usize];
            let lo = list.partition_point(|&e| ((e >> 16) as usize) < m0);
            let hi = list.partition_point(|&e| ((e >> 16) as usize) < m1);
            for &e in &list[lo..hi] {
                rows.push((base[(e >> 16) as usize - m0] + (e & 0xffff) as usize) as u16);
            }
            starts[code + 1] = rows.len() as u32;
        }
        Tile { m0, m1, starts, rows }
    }).collect();
    let mut f = std::io::BufWriter::new(fs::File::create(out).unwrap());
    f.write_all(bytes_of(&[idx.k as u64, tiles.len() as u64])).unwrap();
    for t in &tiles {
        f.write_all(bytes_of(&[t.m0 as u64, t.m1 as u64, t.rows.len() as u64])).unwrap();
        f.write_all(bytes_of(&t.starts)).unwrap();
        f.write_all(bytes_of(&t.rows)).unwrap();
    }
    eprintln!("{} tiles, {} entries", tiles.len(), tiles.iter().map(|t| t.rows.len()).sum::<usize>());
}

fn load_tiles(path: &str) -> TileIndex {
    let mut f = std::io::BufReader::new(fs::File::open(path).unwrap());
    let mut h = [0u64; 2];
    f.read_exact(bytes_of_mut(&mut h)).unwrap();
    let space = 20usize.pow(h[0] as u32);
    let tiles = (0..h[1]).map(|_| {
        let mut t = [0u64; 3];
        f.read_exact(bytes_of_mut(&mut t)).unwrap();
        let mut starts = vec![0u32; space + 1];
        f.read_exact(bytes_of_mut(&mut starts)).unwrap();
        let mut rows = vec![0u16; t[2] as usize];
        f.read_exact(bytes_of_mut(&mut rows)).unwrap();
        Tile { m0: t[0] as usize, m1: t[1] as usize, starts, rows }
    }).collect();
    TileIndex { k: h[0] as usize, tiles }
}

/// One block's postings, LANES at a time: per k-mer, groups of LANES windows stored offset-major
/// (offset o of lane l at o * LANES + l), WIN residues each, 20 outside the ORF, all plus BIAS. Unused
/// lanes hold 20 everywhere, so they score nothing, and orf u32::MAX.
///
/// LANES is the widest register the architecture's kernel fills with one job: 32 bytes on AVX2, where the
/// 16-byte halves of a profile column are broadcast, so nothing is assembled from two jobs. Wider groups
/// are fuller only in larger blocks (hit-weighted, 0.69 full in 300,000 ORFs and 0.93 in 1.6 M).
#[cfg(target_arch = "x86_64")]
const LANES: usize = 32;
#[cfg(not(target_arch = "x86_64"))]
const LANES: usize = 16;
/// Added to every residue in a window, so that `pshufb`, which gives 0 for an index with its top bit set,
/// looks up residues 0..16 as they are (0x70..0x80) and 16..21 with the top bit flipped (0x80..0x85 to
/// 0..5): the two halves of a column are then joined with an OR.
const BIAS: u8 = if cfg!(target_arch = "x86_64") { 0x70 } else { 0 };
/// Bits per lane in a mask, as each architecture's SIMD compare gives them.
const LANE_BITS: usize = if cfg!(target_arch = "x86_64") { 1 } else { 4 };
const LANE: u64 = (1 << LANE_BITS) - 1;
struct Postings {
    gstarts: Vec<u32>, // per k-mer, into groups
    count: Vec<u32>,   // per k-mer, real postings
    windows: Vec<u8>,
    orf: Vec<u32>,
    pos: Vec<u16>,
    covered: Vec<u64>, // per group, a mask of the lanes whose window holds the whole ORF
    span: Vec<(u8, u8)>, // per group, the window offsets that hold a residue in some lane
}

fn postings(orfs: &[&[u8]], k: usize) -> Postings {
    let space = 20usize.pow(k as u32);
    let code_at = |f: &[u8], j: usize| -> Option<usize> {
        let mut code = 0usize;
        for &c in &f[j..j + k] {
            if c >= 20 {
                return None;
            }
            code = code * 20 + c as usize;
        }
        Some(code)
    };
    // Each k-mer's postings sorted by where the ORF sits in the window, so a group's lanes share most of
    // their residue offsets and the kernel can skip offsets that are padding in every lane
    // filled in place at each ORF's offset, so the list is never held twice
    let per_orf: Vec<usize> = orfs.par_iter().map(|f| (0..(f.len() + 1).saturating_sub(k)).filter(|&j| code_at(f, j).is_some()).count()).collect();
    let mut at = Vec::with_capacity(orfs.len() + 1);
    at.push(0usize);
    for n in &per_orf {
        at.push(at.last().unwrap() + n);
    }
    let mut all = vec![(0u32, 0u8, 0u8, 0u32, 0u16); *at.last().unwrap()];
    let mut rest: &mut [(u32, u8, u8, u32, u16)] = &mut all;
    let mut parts = Vec::with_capacity(orfs.len());
    for &n in &per_orf {
        let (head, tail) = std::mem::take(&mut rest).split_at_mut(n);
        parts.push(head);
        rest = tail;
    }
    parts.into_par_iter().enumerate().for_each(|(o, out)| {
        let f = orfs[o];
        let mut i = 0;
        for j in 0..(f.len() + 1).saturating_sub(k) {
            if let Some(c) = code_at(f, j) {
                out[i] = (c as u32, HALF.saturating_sub(j) as u8, (HALF + f.len() - j).min(WIN) as u8, o as u32, j as u16);
                i += 1;
            }
        }
    });
    all.par_sort_unstable();
    let mut count = vec![0u32; space];
    for x in &all {
        count[x.0 as usize] += 1;
    }
    let (mut gstarts, mut first) = (vec![0u32; space + 1], vec![0usize; space + 1]);
    for c in 0..space {
        gstarts[c + 1] = gstarts[c] + count[c].div_ceil(LANES as u32);
        first[c + 1] = first[c] + count[c] as usize;
    }
    let groups = gstarts[space] as usize;
    let mut group_code = vec![0u32; groups];
    for c in 0..space {
        group_code[gstarts[c] as usize..gstarts[c + 1] as usize].fill(c as u32);
    }
    let (mut windows, mut orf, mut pos) = (vec![20 + BIAS; groups * WIN * LANES], vec![u32::MAX; groups * LANES], vec![0u16; groups * LANES]);
    let (mut covered, mut span) = (vec![0u64; groups], vec![(WIN as u8, 0u8); groups]);
    windows
        .par_chunks_mut(WIN * LANES)
        .zip(orf.par_chunks_mut(LANES))
        .zip(pos.par_chunks_mut(LANES))
        .zip(covered.par_iter_mut())
        .zip(span.par_iter_mut())
        .enumerate()
        .for_each(|(g, ((((w, orf), pos), covered), span))| {
            let c = group_code[g] as usize;
            let r0 = (g - gstarts[c] as usize) * LANES;
            let lanes = (count[c] as usize - r0).min(LANES);
            for (lane, &(_, start, end, o, j)) in all[first[c] + r0..first[c] + r0 + lanes].iter().enumerate() {
                orf[lane] = o;
                pos[lane] = j;
                *span = (span.0.min(start), span.1.max(end));
                let f = orfs[o as usize];
                let j = j as usize;
                // the window holds the whole ORF, so its score is the whole diagonal's
                *covered |= (((j <= HALF && f.len() <= j + WIN - HALF) as u64) * LANE) << (LANE_BITS * lane);
                for x in start as usize..end as usize {
                    w[x * LANES + lane] = f[j + x - HALF].min(20) + BIAS;
                }
            }
        });
    Postings { gstarts, count, windows, orf, pos, covered, span }
}

/// Scores four (profile row, group) jobs at once: best local score of each lane's window, in quarter bits
/// saturating at 127, against the profile rows from `row`. Gives two masks per job, the lanes whose best
/// is at least gates.0 and gates.1 (LANE_BITS per lane). Independent chains hide the add-max latency.
///
/// In SIMD a score is held less 128, so the saturating add is also the local score's floor. That lets it
/// run past 127, where the scalar one sticks, but only once it has passed every gate.
type Masks = [(u64, u64); 4];
type Scorer4 = fn(&[i8], [usize; 4], [&[u8]; 4], (usize, usize), (i8, i8), &mut Masks);

/// The oracle: every SIMD version must return the same.
fn score4_scalar(prof: &[i8], rows: [usize; 4], w: [&[u8]; 4], span: (usize, usize), gates: (i8, i8), masks: &mut Masks) {
    for i in 0..4 {
        masks[i] = (0, 0);
        for lane in 0..LANES {
            let (mut cur, mut best) = (0i8, 0i8);
            for o in span.0..span.1 {
                cur = cur.saturating_add(prof[(rows[i] + o) * STRIDE + (w[i][o * LANES + lane] - BIAS) as usize]).max(0);
                best = best.max(cur);
            }
            masks[i].0 |= ((best >= gates.0) as u64 * LANE) << (LANE_BITS * lane);
            masks[i].1 |= ((best >= gates.1) as u64 * LANE) << (LANE_BITS * lane);
        }
    }
}

#[cfg(target_arch = "aarch64")]
fn score4_neon(prof: &[i8], rows: [usize; 4], w: [&[u8]; 4], span: (usize, usize), gates: (i8, i8), masks: &mut Masks) {
    use std::arch::aarch64::*;
    for i in 0..4 {
        assert!((rows[i] + WIN) * STRIDE <= prof.len() && w[i].len() >= WIN * LANES);
    }
    unsafe {
        let (mut cur, mut best) = ([vdupq_n_s8(i8::MIN); 4], [vdupq_n_s8(i8::MIN); 4]);
        let p = rows.map(|r| prof.as_ptr().add(r * STRIDE));
        for o in span.0..span.1 {
            for i in 0..4 {
                let table = vld1q_s8_x2(p[i].add(o * STRIDE));
                let idx = vld1q_u8(w[i].as_ptr().add(o * LANES));
                cur[i] = vqaddq_s8(cur[i], vqtbl2q_s8(table, idx));
                best[i] = vmaxq_s8(best[i], cur[i]);
            }
        }
        let (g0, g1) = (vdupq_n_s8(gates.0 ^ i8::MIN), vdupq_n_s8(gates.1 ^ i8::MIN));
        // compare, then narrow each byte to 4 bits: a 64-bit mask with a nibble per lane
        let mask = |best, gate| vget_lane_u64(vreinterpret_u64_u8(vshrn_n_u16(vreinterpretq_u16_u8(vcgeq_s8(best, gate)), 4)), 0);
        for i in 0..4 {
            masks[i] = (mask(best[i], g0), mask(best[i], g1));
        }
    }
}

/// `pshufb` indexes 16 bytes, so a 32-entry lookup is two shuffles and an OR (see BIAS) where NEON's
/// `vqtbl2q_s8` is one instruction. Here 32 lanes are two passes of 16.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "sse4.1")]
unsafe fn score4_sse41_inner(prof: &[i8], rows: [usize; 4], w: [&[u8]; 4], span: (usize, usize), gates: (i8, i8), masks: &mut Masks) {
    use std::arch::x86_64::*;
    for i in 0..4 {
        assert!((rows[i] + WIN) * STRIDE <= prof.len() && w[i].len() >= WIN * LANES);
    }
    let p = rows.map(|r| prof.as_ptr().add(r * STRIDE));
    let flip = _mm_set1_epi8(i8::MIN);
    let (g0, g1) = (_mm_set1_epi8(gates.0 ^ i8::MIN), _mm_set1_epi8(gates.1 ^ i8::MIN));
    *masks = [(0, 0); 4];
    for half in 0..2 {
        let (mut cur, mut best) = ([flip; 4], [flip; 4]);
        for o in span.0..span.1 {
            for i in 0..4 {
                let col = p[i].add(o * STRIDE) as *const __m128i;
                let idx = _mm_loadu_si128(w[i].as_ptr().add(o * LANES + half * 16) as *const __m128i);
                let lo = _mm_shuffle_epi8(_mm_loadu_si128(col), idx);
                let hi = _mm_shuffle_epi8(_mm_loadu_si128(col.add(1)), _mm_xor_si128(idx, flip));
                cur[i] = _mm_adds_epi8(cur[i], _mm_or_si128(lo, hi));
                best[i] = _mm_max_epi8(best[i], cur[i]);
            }
        }
        // best >= gate, as max(best, gate) == best
        let mask = |best, gate| (_mm_movemask_epi8(_mm_cmpeq_epi8(_mm_max_epi8(best, gate), best)) as u32 as u64) << (16 * half);
        for i in 0..4 {
            masks[i].0 |= mask(best[i], g0);
            masks[i].1 |= mask(best[i], g1);
        }
    }
}

#[cfg(target_arch = "x86_64")]
fn score4_sse41(prof: &[i8], rows: [usize; 4], w: [&[u8]; 4], span: (usize, usize), gates: (i8, i8), masks: &mut Masks) {
    unsafe { score4_sse41_inner(prof, rows, w, span, gates, masks) }
}

/// AVX2: one job's 32 lanes per register. `vpshufb` shuffles each 128-bit half with its own table, so both
/// halves get the same one, broadcast from memory, which costs a load and no shuffle.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn score4_avx2_inner(prof: &[i8], rows: [usize; 4], w: [&[u8]; 4], span: (usize, usize), gates: (i8, i8), masks: &mut Masks) {
    use std::arch::x86_64::*;
    for i in 0..4 {
        assert!((rows[i] + WIN) * STRIDE <= prof.len() && w[i].len() >= WIN * LANES);
    }
    let p = rows.map(|r| prof.as_ptr().add(r * STRIDE));
    let flip = _mm256_set1_epi8(i8::MIN);
    let (mut cur, mut best) = ([flip; 4], [flip; 4]);
    for o in span.0..span.1 {
        for i in 0..4 {
            let col = p[i].add(o * STRIDE) as *const __m128i;
            let idx = _mm256_loadu_si256(w[i].as_ptr().add(o * LANES) as *const __m256i);
            let lo = _mm256_shuffle_epi8(_mm256_broadcastsi128_si256(_mm_loadu_si128(col)), idx);
            let hi = _mm256_shuffle_epi8(_mm256_broadcastsi128_si256(_mm_loadu_si128(col.add(1))), _mm256_xor_si256(idx, flip));
            cur[i] = _mm256_adds_epi8(cur[i], _mm256_or_si256(lo, hi));
            best[i] = _mm256_max_epi8(best[i], cur[i]);
        }
    }
    let (g0, g1) = (_mm256_set1_epi8(gates.0 ^ i8::MIN), _mm256_set1_epi8(gates.1 ^ i8::MIN));
    let mask = |best, gate| _mm256_movemask_epi8(_mm256_cmpeq_epi8(_mm256_max_epi8(best, gate), best)) as u32 as u64;
    for i in 0..4 {
        masks[i] = (mask(best[i], g0), mask(best[i], g1));
    }
}

#[cfg(target_arch = "x86_64")]
fn score4_avx2(prof: &[i8], rows: [usize; 4], w: [&[u8]; 4], span: (usize, usize), gates: (i8, i8), masks: &mut Masks) {
    unsafe { score4_avx2_inner(prof, rows, w, span, gates, masks) }
}

/// Every kernel this CPU has, fastest first. FILTER_SIMD=name forces one.
fn scorers() -> Vec<(&'static str, Scorer4)> {
    let mut all: Vec<(&'static str, Scorer4)> = vec![];
    #[cfg(target_arch = "aarch64")]
    all.push(("neon", score4_neon));
    #[cfg(target_arch = "x86_64")]
    {
        if is_x86_feature_detected!("avx2") {
            all.push(("avx2", score4_avx2));
        }
        if is_x86_feature_detected!("sse4.1") {
            all.push(("sse4.1", score4_sse41));
        }
    }
    all.push(("scalar", score4_scalar));
    all
}

fn scorer4() -> (&'static str, Scorer4) {
    let all = scorers();
    match std::env::var("FILTER_SIMD") {
        Ok(want) if !want.is_empty() => *all.iter().find(|(name, _)| *name == want).unwrap_or_else(|| panic!("no {want} kernel on this CPU")),
        _ => all[0],
    }
}

/// A model's profile in quarter bits, int8, STRIDE bytes per column, with HALF columns of padding before
/// and WIN - HALF after; padding and residue 20 score i8::MIN, which resets a local score.
fn int8_profile(pf: &Pfam, m: usize, out: &mut Vec<i8>) {
    let (a, b) = (pf.offs[m] as usize, pf.offs[m + 1] as usize);
    let at = out.len();
    out.resize(at + (b - a + WIN) * STRIDE, i8::MIN);
    let out = &mut out[at..];
    for c in 0..b - a {
        for r in 0..20 {
            out[(c + HALF) * STRIDE + r] = (pf.lo[(a + c) * 21 + r] * 4.0).round().clamp(-127.0, 127.0) as i8;
        }
    }
}

/// Whole-diagonal ungapped score in bits, as search computes it, for read position j at model column col.
fn diagonal(pf: &Pfam, m: usize, col: usize, f: &[u8], j: usize) -> f32 {
    let a = pf.offs[m] as usize;
    let cols = pf.offs[m + 1] as usize - a;
    let back = j.min(col);
    let len = (f.len() - j + back).min(cols - col + back);
    let mut at = (a + col - back) * 21;
    let (mut cur, mut best) = (0.0f32, 0.0f32);
    for &r in &f[j - back..j - back + len] {
        cur = fmax(cur + pf.lo[at + r as usize], 0.0);
        best = fmax(best, cur);
        at += 21;
    }
    best
}

fn join(pf: &Pfam, ti: &TileIndex, path: &str, bar: f32, gate_bits: f32, block: usize, out: &str) {
    use std::io::BufRead;
    let mut lines = std::io::BufReader::with_capacity(1 << 20, fs::File::open(path).unwrap()).lines();
    let mut n_orfs = 0usize;
    let gate = (gate_bits * 4.0).round() as i8;
    // a survivor is model << 48 | ORF << 17 | diagonal key, so it sorts as one integer
    assert!(pf.models() < 1 << 16 && block < 1 << 31);
    // gapped score in a band of +-8 around the best diagonal (step 4); FILTER_BAND=b or =whole to change it
    let band: Option<i64> = match std::env::var("FILTER_BAND").as_deref() { Ok("whole") => None, Ok(v) => Some(v.parse().unwrap()), Err(_) => Some(8) };
    // FILTER_DROP=n: skip the n k-mers with the longest lists in the index (promiscuous k-mers, from the
    // index alone)
    let space = ti.tiles[0].starts.len() - 1;
    let mut dropped = vec![false; space];
    if let Ok(n) = std::env::var("FILTER_DROP") {
        let mut len = vec![0u64; space];
        for t in &ti.tiles {
            for c in 0..space {
                len[c] += (t.starts[c + 1] - t.starts[c]) as u64;
            }
        }
        let mut order: Vec<usize> = (0..space).collect();
        order.sort_unstable_by_key(|&c| std::cmp::Reverse(len[c]));
        for &c in &order[..n.parse::<usize>().unwrap()] {
            dropped[c] = true;
        }
    }
    let covered_gate = (((bar - SLACK) * 4.0).floor().min(127.0) as i8).max(gate);
    // FILTER_MIN_GAPPED=g: write only pairs whose gapped score is at least g (HMMER's bar), so pairs scored
    // at the diagonal bar but below the hand-off are not formatted and written (324 M rows at diagonal 17).
    // Half a printed digit lower, so a score that prints as g is kept and readers' cut on the text decides.
    let min_gapped: f32 = std::env::var("FILTER_MIN_GAPPED").map_or(f32::NEG_INFINITY, |v| v.parse::<f32>().unwrap() - 0.005);
    let bound = if band == Some(BAND) { Bound::new(pf, min_gapped) } else { None };
    let bounded = AtomicU64::new(0); // pairs the bound lets through to the f32 score
    let (engine, score4) = scorer4();
    let (hits, survivors, pairs_n) = (AtomicU64::new(0), AtomicU64::new(0), AtomicU64::new(0));
    let (surv_ns, unique, gap_ns) = (AtomicU64::new(0), AtomicU64::new(0), AtomicU64::new(0)); // surv: sort, dedupe and exact diagonals
    let (mut t_post, mut t_join) = (0.0f64, 0.0f64);
    let mut sink = std::io::BufWriter::new(fs::File::create(out).unwrap());
    let threads = rayon::current_num_threads() as f64;
    let t0 = Instant::now();
    loop {
        let orfs = Packed::next(&mut lines, block);
        if orfs.len() == 0 {
            break;
        }
        n_orfs += orfs.len();
        let chunk: Vec<&[u8]> = (0..orfs.len()).map(|i| orfs.seq(i)).collect();
        let t = Instant::now();
        let p = postings(&chunk, ti.k);
        t_post += t.elapsed().as_secs_f64();
        let real: u64 = p.count.iter().map(|&c| c as u64).sum();
        eprintln!("block of {} ORFs: {} postings in {} groups, lanes filled {:.3}, {:.0} MB", chunk.len(), real,
            p.gstarts.last().unwrap(), real as f64 / (*p.gstarts.last().unwrap() as f64 * LANES as f64), p.windows.len() as f64 / 1e6);
        let t = Instant::now();
        // each tile's pairs are scored and formatted before the next tile's are found, so a block holds its
        // text, not its pairs
        let rows: Vec<String> = ti.tiles
            .par_iter()
            .map_init(|| (Vec::<i8>::new(), Vec::<u64>::new(), vec![], vec![]), |(prof, surv, h, del), tile| {
                let base = tile_bases(pf, tile.m0, tile.m1);
                prof.clear();
                for m in tile.m0..tile.m1 {
                    int8_profile(pf, m, prof);
                }
                // model of each padded row
                let mut row_model = vec![0u32; *base.last().unwrap()];
                for m in tile.m0..tile.m1 {
                    row_model[base[m - tile.m0]..base[m - tile.m0 + 1]].fill(m as u32);
                }
                surv.clear();
                let mut n_hits = 0u64;
                let mut jobs = [(0usize, 0usize); 4];
                let mut queued = 0usize;
                let mut masks: Masks = [(0, 0); 4];
                let mut span = (WIN, 0usize);
                let drain = |jobs: &[(usize, usize)], masks: &Masks, surv: &mut Vec<u64>| {
                    for (i, &(row, g)) in jobs.iter().enumerate() {
                        // lanes that pass: the gate, or for a window holding its whole ORF the bar less
                        // SLACK, as that window already scores the whole diagonal to within rounding
                        let (pass, high) = masks[i];
                        let mut mask = pass & (!p.covered[g] | high);
                        if mask == 0 {
                            continue;
                        }
                        let m = row_model[row];
                        let c = row - base[m as usize - tile.m0];
                        while mask != 0 {
                            let lane = mask.trailing_zeros() as usize / LANE_BITS;
                            mask &= !(LANE << (LANE_BITS * lane));
                            let q = g * LANES + lane;
                            // diagonal key: column minus read position, offset to stay positive
                            let key = (c as u32 + 65536) - p.pos[q] as u32;
                            surv.push((m as u64) << 48 | (p.orf[q] as u64) << 17 | key as u64);
                        }
                    }
                    // a diagonal survives once per k-mer on it: compact before the list grows large
                    if surv.len() > SURV_COMPACT {
                        surv.sort_unstable();
                        surv.dedup();
                    }
                };
                for code in 0..tile.starts.len() - 1 {
                    let (r0, r1) = (tile.starts[code] as usize, tile.starts[code + 1] as usize);
                    let (g0, g1) = (p.gstarts[code] as usize, p.gstarts[code + 1] as usize);
                    if r0 == r1 || g0 == g1 || dropped[code] {
                        continue;
                    }
                    n_hits += (r1 - r0) as u64 * p.count[code] as u64;
                    // the list is read from memory once here, then stays in cache for every row;
                    // one group against four rows at a time, so a batch shares the group's span
                    for g in g0..g1 {
                        // offset 0 of the window sits at column c - HALF, padded row c
                        for &row in &tile.rows[r0..r1] {
                            jobs[queued] = (row as usize, g);
                            queued += 1;
                            span = (span.0.min(p.span[g].0 as usize), span.1.max(p.span[g].1 as usize));
                            if queued == 4 {
                                let w = jobs.map(|(_, g)| &p.windows[g * WIN * LANES..(g + 1) * WIN * LANES]);
                                score4(prof, jobs.map(|(r, _)| r), w, span, (gate, covered_gate), &mut masks);
                                drain(&jobs, &masks, surv);
                                queued = 0;
                                span = (WIN, 0);
                            }
                        }
                    }
                }
                if queued > 0 {
                    // the last few jobs: the first pads them to four
                    for j in queued..4 {
                        jobs[j] = jobs[0];
                    }
                    let w = jobs.map(|(_, g)| &p.windows[g * WIN * LANES..(g + 1) * WIN * LANES]);
                    score4(prof, jobs.map(|(r, _)| r), w, span, (gate, covered_gate), &mut masks);
                }
                drain(&jobs[..queued], &masks, surv);
                hits.fetch_add(n_hits, Relaxed);
                survivors.fetch_add(surv.len() as u64, Relaxed);
                let ts = Instant::now();
                surv.sort_unstable();
                surv.dedup();
                unique.fetch_add(surv.len() as u64, Relaxed);
                let mut res: Vec<(u32, u32, f32, i64)> = vec![];
                for &x in surv.iter() {
                    let (m, o, key) = ((x >> 48) as u32, (x >> 17) as u32 & 0x7fff_ffff, x as u32 & 0x1ffff);
                    let f = chunk[o as usize];
                    let cols = (pf.offs[m as usize + 1] - pf.offs[m as usize]) as usize;
                    let d = key as i64 - 65536; // column - position
                    // a point on the diagonal inside both the read and the model
                    let j = if d >= 0 { 0usize } else { (-d) as usize };
                    let col = (j as i64 + d) as usize;
                    if j >= f.len() || col >= cols {
                        continue;
                    }
                    let s = diagonal(pf, m as usize, col, f, j);
                    if s >= bar {
                        match res.last_mut() {
                            Some(x) if x.0 == o && x.1 == m => {
                                if s > x.2 {
                                    (x.2, x.3) = (s, d);
                                }
                            }
                            _ => res.push((o, m, s, d)),
                        }
                    }
                }
                surv_ns.fetch_add(ts.elapsed().as_nanos() as u64, Relaxed);
                let ts = Instant::now();
                use std::fmt::Write;
                let (mut text, mut n) = (String::new(), 0u64);
                for four in res.chunks(PAIRS) {
                    let reach = bound.as_ref().map_or([true; PAIRS], |b| {
                        // a short last chunk repeats its first pair
                        b.may_reach(std::array::from_fn(|i| {
                            let &(o, m, _, diag) = four.get(i).unwrap_or(&four[0]);
                            (m as usize, chunk[o as usize], diag)
                        }))
                    });
                    for (&(o, m, d, diag), _) in four.iter().zip(reach).filter(|(_, r)| *r) {
                        let (name, f) = (orfs.name(o as usize), chunk[o as usize]);
                        n += 1;
                        let a = pf.offs[m as usize] as usize;
                        let cols = pf.offs[m as usize + 1] as usize - a;
                        let g = match band {
                            Some(b) => gapped_band(&pf.lo, a, cols, f, diag, b, h, del),
                            None => gapped(&pf.lo, a, cols, f, h, del),
                        };
                        if g >= min_gapped {
                            writeln!(text, "{name}\t{m}\t{d:.2}\t{g:.2}").unwrap();
                        }
                    }
                }
                bounded.fetch_add(n, Relaxed);
                pairs_n.fetch_add(res.len() as u64, Relaxed);
                gap_ns.fetch_add(ts.elapsed().as_nanos() as u64, Relaxed);
                text
            })
            .collect();
        t_join += t.elapsed().as_secs_f64();
        for r in rows {
            sink.write_all(r.as_bytes()).unwrap();
        }
    }
    sink.flush().unwrap();
    let wall = t0.elapsed().as_secs_f64();
    let h = hits.load(Relaxed) as f64;
    let (surv_s, gap_s) = (surv_ns.load(Relaxed) as f64 / 1e9, gap_ns.load(Relaxed) as f64 / 1e9);
    let join_s = t_join * threads - gap_s; // the survivor stage is part of it, and so is waiting for the last tile
    eprintln!(
        "{engine}: {} ORFs, {:.1}s wall on {threads} threads, {:.0} core-s; postings {:.1} s wall, join {:.1} core-s, gapped {:.1} core-s; \
         {} hits, {:.1} ns per hit (join core time), {} survivors ({:.4}), {} unique diagonals, survivor stage {:.1} core-s, {} pairs, {} scored gapped",
        n_orfs, wall, wall * threads, t_post, join_s, gap_s,
        h as u64, join_s * 1e9 / h.max(1.0), survivors.load(Relaxed), survivors.load(Relaxed) as f64 / h.max(1.0), unique.load(Relaxed), surv_s, pairs_n.load(Relaxed), bounded.load(Relaxed)
    );
}

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let pf = Pfam::load(&a[2]);
    match a[1].as_str() {
        "build" => build(&pf, a[3].parse().unwrap(), a[4].parse().unwrap(), &a[5]),
        "search" => search(&pf, &load_index(&a[3]), &a[4], a[5].parse().unwrap(), &a[6], a.get(7)),
        "window" => window(&pf, &load_index(&a[3]), &a[4], a[5].parse().unwrap(), &a[6]),
        "tile" => tile(&pf, &load_index(&a[3]), &a[4]),
        "join" => join(&pf, &load_tiles(&a[3]), &a[7], a[4].parse().unwrap(), a[5].parse().unwrap(), a[6].parse().unwrap(), &a[8]),
        "sketch" => sketch(&pf, &a[3], a[4].parse().unwrap(), a[5].parse().unwrap(), a[6].parse().unwrap(), a[7] == "rev", &a[8]),
        _ => panic!("build | search | sketch"),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn enumerate_prunes_prefixes_and_reports_scores() {
        // Two columns: residue 0 scores 5 in each, all else -9.
        let mut lo = vec![-9.0f32; 42];
        lo[0] = 5.0;
        lo[21] = 5.0;
        let mut seen = vec![];
        enumerate(&lo, 0, 2, 9.0, &|_, _| true, &mut |c, s| seen.push((c, s)));
        assert_eq!(seen, vec![(0, 10.0)]);
        seen.clear();
        enumerate(&lo, 0, 2, 9.0, &|i, _| i < 2, &mut |c, s| seen.push((c, s)));
        assert!(seen.is_empty());
    }

    #[test]
    fn every_kernel_equals_scalar() {
        let mut x = 99u64;
        let mut next = || { x ^= x << 13; x ^= x >> 7; x ^= x << 17; x };
        // scores up to 60 a column, so some lanes saturate
        let prof: Vec<i8> = (0..3 * WIN * STRIDE).map(|i| if i % STRIDE >= 20 { i8::MIN } else { (next() % 90) as i8 - 30 }).collect();
        for t in 0..400 {
            // residues and padding, or only well-scoring ones
            let ws: Vec<Vec<u8>> = (0..4).map(|_| (0..WIN * LANES).map(|_| (next() % 21) as u8 + BIAS).collect()).collect();
            let w = [&ws[0][..], &ws[1][..], &ws[2][..], &ws[3][..]];
            let rows = [0, 1, 2, 3].map(|_| (next() % (2 * WIN as u64)) as usize);
            let (span, gates) = (((next() % 10) as usize, WIN - (next() % 10) as usize), ([66, 10, 127, 0][t % 4], [68, 127, 100, 40][t % 4]));
            let mut want: Masks = [(0, 0); 4];
            score4_scalar(&prof, rows, w, span, gates, &mut want);
            for (name, kernel) in scorers() {
                let mut got: Masks = [(0, 0); 4];
                kernel(&prof, rows, w, span, gates, &mut got);
                assert_eq!(got, want, "{name}");
            }
        }
    }

    #[test]
    fn bound_never_drops_what_the_band_keeps() {
        let mut x = 7u64;
        let mut next = || { x ^= x << 13; x ^= x >> 7; x ^= x << 17; x };
        // models with a strong residue per column; reads that copy a stretch with gaps and changes
        let lens = [9usize, 40, 130];
        let mut offs = vec![0u32];
        let mut lo = vec![];
        let mut top: Vec<Vec<u8>> = vec![];
        for n in lens {
            let t: Vec<u8> = (0..n).map(|_| (next() % 20) as u8).collect();
            for &best in &t {
                lo.extend((0..21u8).map(|r| if r == best { 2.0 + (next() % 200) as f32 / 100.0 } else { -(next() as u32 as f32 / u32::MAX as f32) * 3.0 }));
            }
            offs.push(offs.last().unwrap() + n as u32);
            top.push(t);
        }
        let pf = Pfam { lo, offs };
        let (mut kept, mut dropped) = (0, 0);
        for min in [17.0f32, 21.25, 26.25] {
            let Some(b) = Bound::new(&pf, min - 0.005) else { return }; // no SIMD here
            let (mut h, mut del) = (vec![], vec![]);
            for _ in 0..3000 {
                let mut pairs: Vec<(usize, Vec<u8>, i64)> = vec![];
                for _ in 0..PAIRS {
                    let m = (next() % 3) as usize;
                    let (mut f, mut c) = (vec![], (next() % lens[m] as u64) as usize);
                    let d = c as i64 - (next() % 5) as i64 + (next() % 21) as i64 - 10;
                    while f.len() < 1 + (next() % 60) as usize {
                        match next() % 12 {
                            0 => c += (next() % 4) as usize,
                            1..=3 => f.push((next() % 21) as u8),
                            _ => f.push(top[m].get(c).copied().unwrap_or(20)),
                        }
                        c += 1;
                    }
                    pairs.push((m, f, d));
                }
                let reach = b.may_reach(std::array::from_fn(|i| (pairs[i].0, &pairs[i].1[..], pairs[i].2)));
                for ((m, f, d), reach) in pairs.iter().zip(reach) {
                    let g = gapped_band(&pf.lo, pf.offs[*m] as usize, lens[*m], f, *d, BAND, &mut h, &mut del);
                    assert!(reach || g < min - 0.005, "dropped {g} at {min}");
                    (kept, dropped) = (kept + (g >= min) as usize, dropped + !reach as usize);
                }
            }
        }
        assert!(kept > 100 && dropped > 100, "{kept} kept, {dropped} dropped: the test has stopped testing");
    }

    #[test]
    fn band_wide_enough_equals_whole_matrix() {
        let mut x = 7u64;
        let mut next = || { x ^= x << 13; x ^= x >> 7; x ^= x << 17; x };
        let (mut h, mut del) = (vec![], vec![]);
        for _ in 0..300 {
            let cols = 5 + (next() % 60) as usize;
            let lo: Vec<f32> = (0..cols * 21).map(|_| (next() % 100) as f32 / 10.0 - 6.0).collect();
            let f: Vec<u8> = (0..3 + next() % 40).map(|_| (next() % 20) as u8).collect();
            let whole = gapped(&lo, 0, cols, &f, &mut h, &mut del);
            let wide = gapped_band(&lo, 0, cols, &f, 0, 200, &mut h, &mut del);
            assert!((whole - wide).abs() < 1e-4, "{whole} {wide}");
            for d in [-10i64, 0, 7] {
                for b in [0i64, 3, 8] {
                    // reference: the whole matrix with cells outside the band unaligned
                    let (mut hh, mut dd) = (vec![0.0f32; cols + 1], vec![f32::NEG_INFINITY; cols + 1]);
                    let mut want = 0.0f32;
                    for (i, &r) in f.iter().enumerate() {
                        let (mut diag, mut left, mut ins) = (0.0f32, 0.0f32, f32::NEG_INFINITY);
                        for c in 1..=cols {
                            let up = hh[c];
                            if ((c as i64 - 1) - i as i64 - d).abs() > b {
                                (hh[c], dd[c], diag, left, ins) = (0.0, f32::NEG_INFINITY, up, 0.0, f32::NEG_INFINITY);
                                continue;
                            }
                            dd[c] = (up - GAP_OPEN).max(dd[c] - GAP_EXTEND);
                            ins = (left - GAP_OPEN).max(ins - GAP_EXTEND);
                            let v = (diag + lo[(c - 1) * 21 + r as usize]).max(dd[c]).max(ins).max(0.0);
                            diag = up;
                            hh[c] = v;
                            left = v;
                            want = want.max(v);
                        }
                    }
                    let got = gapped_band(&lo, 0, cols, &f, d, b, &mut h, &mut del);
                    assert!((got - want).abs() < 1e-4, "d {d} b {b}: {got} {want}");
                }
            }
        }
    }

    #[test]
    fn sketch_keeps_the_stated_fraction() {
        let kept = (0..1_000_000).filter(|&c| in_sketch(c, 8)).count();
        assert!((kept as f64 / 125_000.0 - 1.0).abs() < 0.02, "{kept}");
        assert!(in_sketch(12345, 1));
    }
}
