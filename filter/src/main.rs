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
        lo.par_chunks_mut(21).for_each(score_unknown);
        Pfam { lo, offs }
    }
    fn models(&self) -> usize {
        self.offs.len() - 1
    }
}

/// Scores an unknown residue, X, in one column as HMMER does: the column's mean score over the background. prep.py
/// writes -100, which ended every alignment at an X that HMMER carries on through.
fn score_unknown(column: &mut [f32]) {
    let mean = (0..20).map(|r| Q[r] * column[r] as f64).sum::<f64>() / Q.iter().sum::<f64>();
    column[20] = mean as f32;
}

/// Which residues of a stretch make up the word an index holds: bit i for residue i. Contiguous for a k-mer;
/// a spaced word leaves residues out between its first and last, which are both in it.
#[derive(Clone, Copy, PartialEq, Debug)]
struct Shape(u64);
impl Shape {
    /// "4" for a 4-mer, "1011001" for a spaced word.
    fn parse(s: &str) -> Shape {
        if s.len() == 1 {
            let k: u32 = s.parse().expect("k, or a shape of 0s and 1s");
            return Shape((1 << k) - 1);
        }
        assert!(s.len() <= 16 && s.starts_with('1') && s.ends_with('1') && s.bytes().all(|b| b == b'0' || b == b'1'), "{s}: a shape is 0s and 1s, from a 1 to a 1");
        Shape(s.bytes().enumerate().map(|(i, b)| ((b == b'1') as u64) << i).sum())
    }
    /// Residues in the word.
    fn k(self) -> usize { self.0.count_ones() as usize }
    /// Residues from its first to its last.
    fn span(self) -> usize { 64 - self.0.leading_zeros() as usize }
    fn offsets(self) -> Vec<usize> { (0..self.span()).filter(|i| self.0 >> i & 1 == 1).collect() }
    fn contiguous(self) -> bool { self.k() == self.span() }
    /// In an index's header: k for a k-mer, as indexes have always had, or the mask above any k.
    fn header(self) -> u64 { if self.contiguous() { self.k() as u64 } else { self.0 << 8 } }
    fn from_header(h: u64) -> Shape { if h < 256 { Shape((1 << h) - 1) } else { Shape(h >> 8) } }
    /// k, for the steps that read k-mers only.
    fn kmer(self) -> usize {
        assert!(self.contiguous(), "a spaced index is for tile and join only");
        self.k()
    }
}

/// Calls `visit` with each word scoring >= t in the columns `cols` of the profile, and its score.
/// `keep(i, prefix)` can drop a prefix of length i, and everything under it.
fn enumerate(lo: &[f32], cols: &[usize], t: f32, keep: &impl Fn(usize, usize) -> bool, visit: &mut impl FnMut(usize, f32)) {
    let mut rest = [0.0f32; 16];
    for i in (0..cols.len()).rev() {
        let top = (0..20).map(|r| lo[cols[i] * 21 + r]).fold(f32::MIN, f32::max);
        rest[i] = rest[i + 1] + top;
    }
    #[allow(clippy::too_many_arguments)]
    fn go(lo: &[f32], cols: &[usize], i: usize, code: usize, sum: f32, t: f32, rest: &[f32], keep: &impl Fn(usize, usize) -> bool, visit: &mut impl FnMut(usize, f32)) {
        if i == cols.len() {
            visit(code, sum);
            return;
        }
        for r in 0..20 {
            let s = sum + lo[cols[i] * 21 + r];
            let next = code * 20 + r;
            if s + rest[i + 1] >= t && keep(i + 1, next) {
                go(lo, cols, i + 1, next, s, t, rest, keep, visit);
            }
        }
    }
    go(lo, cols, 0, 0, 0.0, t, &rest, keep, visit);
}

fn build(pf: &Pfam, shape: Shape, t_all: f32, out: &str) {
    let (k, span, offsets) = (shape.k(), shape.span(), shape.offsets());
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
        (a, if b - a >= span { b - a - span + 1 } else { 0 })
    };
    let at = |c: usize| offsets.iter().map(|o| c + o).collect::<Vec<_>>();
    for m in 0..pf.models() {
        let (a, n) = windows(m);
        for c in 0..n {
            enumerate(&pf.lo, &at(a + c), per_model[m], &|_, _| true, &mut |code, _| starts[code + 1] += 1);
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
            enumerate(&pf.lo, &at(a + c), per_model[m], &|_, _| true, &mut |code, _| {
                entries[cursor[code] as usize] = (m as u32) << 16 | c as u32;
                cursor[code] += 1;
            });
        }
    }
    let mut f = std::io::BufWriter::new(fs::File::create(out).unwrap());
    f.write_all(&shape.header().to_le_bytes()).unwrap();
    f.write_all(&(total as u64).to_le_bytes()).unwrap();
    f.write_all(bytes_of(&starts)).unwrap();
    f.write_all(bytes_of(&entries)).unwrap();
    // dropping a BufWriter ignores the last write's error: a full disk left a short index, read as garbage
    f.flush().unwrap();
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
    shape: Shape,
    starts: Vec<u64>,
    entries: Vec<u32>,
}
fn load_index(path: &str) -> Index {
    let mut f = fs::File::open(path).unwrap();
    let mut h = [0u64; 2];
    f.read_exact(bytes_of_mut(&mut h)).unwrap();
    let shape = Shape::from_header(h[0]);
    let mut starts = vec![0u64; 20usize.pow(shape.k() as u32) + 1];
    f.read_exact(bytes_of_mut(&mut starts)).unwrap();
    let mut entries = vec![0u32; h[1] as usize];
    f.read_exact(bytes_of_mut(&mut entries)).unwrap();
    Index { shape, starts, entries }
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

/// What the gapped bound needs of a SIMD register: LANES signed bytes.
#[cfg(any(target_arch = "aarch64", target_arch = "x86_64"))]
trait Bytes: Copy {
    const LANES: usize;
    unsafe fn load(p: *const i8) -> Self;
    unsafe fn splat(v: i8) -> Self;
    unsafe fn adds(self, o: Self) -> Self;
    unsafe fn subs(self, o: Self) -> Self;
    unsafe fn max(self, o: Self) -> Self;
    unsafe fn min(self, o: Self) -> Self;
    /// Every byte N lanes up, N a power of two of at most half the lanes, with `below`'s last N bytes under them.
    unsafe fn raise<const N: usize>(self, below: Self) -> Self;
    /// Every byte one lane down, with `above`'s first byte over them.
    unsafe fn lower(self, above: Self) -> Self;
    /// The last byte in the first lane.
    unsafe fn last_first(self) -> Self;
    unsafe fn hmax(self) -> i8;
}
#[cfg(target_arch = "aarch64")]
impl Bytes for std::arch::aarch64::int8x16_t {
    const LANES: usize = 16;
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
    unsafe fn raise<const N: usize>(self, below: Self) -> Self {
        use std::arch::aarch64::vextq_s8;
        match N {
            1 => vextq_s8::<15>(below, self),
            2 => vextq_s8::<14>(below, self),
            4 => vextq_s8::<12>(below, self),
            8 => vextq_s8::<8>(below, self),
            _ => unreachable!(),
        }
    }
    #[inline(always)]
    unsafe fn lower(self, above: Self) -> Self { std::arch::aarch64::vextq_s8::<1>(self, above) }
    #[inline(always)]
    unsafe fn last_first(self) -> Self { std::arch::aarch64::vextq_s8::<15>(self, self) }
    #[inline(always)]
    unsafe fn hmax(self) -> i8 { std::arch::aarch64::vmaxvq_s8(self) }
}
/// Two registers as one, for a CPU without a wider one.
#[cfg(target_arch = "aarch64")]
#[derive(Clone, Copy)]
struct Two<V>(V, V);
#[cfg(target_arch = "aarch64")]
impl<V: Bytes> Bytes for Two<V> {
    const LANES: usize = 2 * V::LANES;
    #[inline(always)]
    unsafe fn load(p: *const i8) -> Self { Two(V::load(p), V::load(p.add(V::LANES))) }
    #[inline(always)]
    unsafe fn splat(v: i8) -> Self { Two(V::splat(v), V::splat(v)) }
    #[inline(always)]
    unsafe fn adds(self, o: Self) -> Self { Two(self.0.adds(o.0), self.1.adds(o.1)) }
    #[inline(always)]
    unsafe fn subs(self, o: Self) -> Self { Two(self.0.subs(o.0), self.1.subs(o.1)) }
    #[inline(always)]
    unsafe fn max(self, o: Self) -> Self { Two(self.0.max(o.0), self.1.max(o.1)) }
    #[inline(always)]
    unsafe fn min(self, o: Self) -> Self { Two(self.0.min(o.0), self.1.min(o.1)) }
    #[inline(always)]
    unsafe fn raise<const N: usize>(self, below: Self) -> Self {
        if N == V::LANES { Two(below.1, self.0) } else { Two(self.0.raise::<N>(below.1), self.1.raise::<N>(self.0)) }
    }
    #[inline(always)]
    unsafe fn lower(self, above: Self) -> Self { Two(self.0.lower(self.1), self.1.lower(above.0)) }
    #[inline(always)]
    unsafe fn last_first(self) -> Self { Two(self.1.last_first(), self.1) }
    #[inline(always)]
    unsafe fn hmax(self) -> i8 { self.0.max(self.1).hmax() }
}
#[cfg(target_arch = "x86_64")]
impl Bytes for std::arch::x86_64::__m128i {
    const LANES: usize = 16;
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
    unsafe fn raise<const N: usize>(self, below: Self) -> Self {
        use std::arch::x86_64::_mm_alignr_epi8;
        match N {
            1 => _mm_alignr_epi8::<15>(self, below),
            2 => _mm_alignr_epi8::<14>(self, below),
            4 => _mm_alignr_epi8::<12>(self, below),
            8 => _mm_alignr_epi8::<8>(self, below),
            _ => unreachable!(),
        }
    }
    #[inline(always)]
    unsafe fn lower(self, above: Self) -> Self { std::arch::x86_64::_mm_alignr_epi8::<1>(above, self) }
    #[inline(always)]
    unsafe fn last_first(self) -> Self { std::arch::x86_64::_mm_srli_si128::<15>(self) }
    #[inline(always)]
    unsafe fn hmax(self) -> i8 {
        let mut b = [0i8; 16];
        std::arch::x86_64::_mm_storeu_si128(b.as_mut_ptr() as *mut _, self);
        b.into_iter().max().unwrap()
    }
}
#[cfg(target_arch = "x86_64")]
impl Bytes for std::arch::x86_64::__m256i {
    const LANES: usize = 32;
    #[inline(always)]
    unsafe fn load(p: *const i8) -> Self { std::arch::x86_64::_mm256_loadu_si256(p as *const _) }
    #[inline(always)]
    unsafe fn splat(v: i8) -> Self { std::arch::x86_64::_mm256_set1_epi8(v) }
    #[inline(always)]
    unsafe fn adds(self, o: Self) -> Self { std::arch::x86_64::_mm256_adds_epi8(self, o) }
    #[inline(always)]
    unsafe fn subs(self, o: Self) -> Self { std::arch::x86_64::_mm256_subs_epi8(self, o) }
    #[inline(always)]
    unsafe fn max(self, o: Self) -> Self { std::arch::x86_64::_mm256_max_epi8(self, o) }
    #[inline(always)]
    unsafe fn min(self, o: Self) -> Self { std::arch::x86_64::_mm256_min_epi8(self, o) }
    #[inline(always)]
    unsafe fn raise<const N: usize>(self, below: Self) -> Self {
        use std::arch::x86_64::{_mm256_alignr_epi8, _mm256_permute2x128_si256};
        // alignr works within each half, so it needs what is under each: below's upper half, self's lower one
        let under = _mm256_permute2x128_si256::<0x21>(below, self);
        match N {
            1 => _mm256_alignr_epi8::<15>(self, under),
            2 => _mm256_alignr_epi8::<14>(self, under),
            4 => _mm256_alignr_epi8::<12>(self, under),
            8 => _mm256_alignr_epi8::<8>(self, under),
            16 => under,
            _ => unreachable!(),
        }
    }
    #[inline(always)]
    unsafe fn lower(self, above: Self) -> Self {
        use std::arch::x86_64::{_mm256_alignr_epi8, _mm256_permute2x128_si256};
        _mm256_alignr_epi8::<1>(_mm256_permute2x128_si256::<0x21>(self, above), self)
    }
    #[inline(always)]
    unsafe fn last_first(self) -> Self {
        use std::arch::x86_64::{_mm256_permute2x128_si256, _mm256_srli_si256};
        _mm256_permute2x128_si256::<0x11>(_mm256_srli_si256::<15>(self), self)
    }
    #[inline(always)]
    unsafe fn hmax(self) -> i8 {
        use std::arch::x86_64::{_mm256_castsi256_si128, _mm256_extracti128_si256};
        Bytes::max(_mm256_castsi256_si128(self), _mm256_extracti128_si256::<1>(self)).hmax()
    }
}

/// An upper bound on gapped_band's score, in bytes, one read residue per step with the band in one
/// register, so a pair that cannot reach `min` is dropped without the f32 matrix. Scores are rounded up,
/// gap costs down, and a band of 8 or 16, a column wider than the register, has its last two columns merged, so it
/// never drops a pair gapped_band would keep. Merged columns score each the other's residues, so that on a whole
/// protein the last lane climbs by chance: a band of 7 or 15 has none, and its bound drops three to seven times as
/// many pairs there. The band is 8 or 16 columns either way, the narrower that holds the one
/// asked for: a wider band only raises a score.
struct Bound {
    prof: Vec<i8>, // per model, 21 residues x (PAD + columns + PAD), residue-major
    at: Vec<usize>,
    scale: f32,
    pass: Vec<i8>, // per model, the lowest bound that may reach its min, offset by -128
    wide: bool, // a band of 16, not 8
    merge: bool, // the band has a column more than there are lanes, and the last lane stands for two
}
const PAD: usize = 40; // columns either side of a model: the band may hang off both ends
impl Bound {
    /// None when there is nothing to bound (no min), or no SIMD for it or for a band this wide here. `mins` has one
    /// per model, and the highest sets the scale, so every one of them fits in a byte.
    fn new(pf: &Pfam, mins: &[f32], band: i64) -> Option<Bound> {
        let min = mins.iter().fold(f32::NEG_INFINITY, |a, &b| a.max(b));
        if !(16.0..=250.0).contains(&min) || !(0..=16).contains(&band) || std::env::var("FILTER_SIMD").as_deref() == Ok("scalar") {
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
        let pass = mins.iter().map(|min| ((min * scale - 0.01).floor() - 128.0) as i8).collect();
        Some(Bound { prof, at, scale, pass, wide: band > 8, merge: band == 8 || band == 16 })
    }
    /// Whether each pair (model, read, diagonal as column - read position) may reach min.
    fn may_reach(&self, pairs: [(usize, &[u8], i64); PAIRS]) -> [bool; PAIRS] {
        let jobs = pairs.map(|(m, f, d)| (&self.prof[self.at[m]..self.at[m + 1]], f, d));
        #[cfg(target_arch = "aarch64")]
        {
            use std::arch::aarch64::int8x16_t;
            // SAFETY: NEON is part of aarch64, and a job's rows have PAD columns either side
            let best = unsafe { if self.wide { bound::<Two<int8x16_t>>(jobs, self.scale, self.merge) } else { bound::<int8x16_t>(jobs, self.scale, self.merge) } };
            return std::array::from_fn(|i| best[i] >= self.pass[pairs[i].0]);
        }
        #[cfg(target_arch = "x86_64")]
        {
            // AVX2 for the wide band's registers, and for its three-operand encoding in the narrow one
            #[target_feature(enable = "avx2")]
            unsafe fn run(jobs: [Job; PAIRS], scale: f32, wide: bool, merge: bool) -> [i8; PAIRS] {
                use std::arch::x86_64::{__m128i, __m256i};
                if wide { bound::<__m256i>(jobs, scale, merge) } else { bound::<__m128i>(jobs, scale, merge) }
            }
            // SAFETY: new() checked for AVX2, and a job's rows have PAD columns either side
            let best = unsafe { run(jobs, self.scale, self.wide, self.merge) };
            return std::array::from_fn(|i| best[i] >= self.pass[pairs[i].0]);
        }
        #[allow(unreachable_code)]
        [true; PAIRS]
    }
}
const PAIRS: usize = 4; // scored side by side, as each pair's rows are one long dependency chain
type Job<'a> = (&'a [i8], &'a [u8], i64); // a model's rows of Bound::prof, read, diagonal

/// The band as lanes, half of them either side: lane l is column i + d - half + l at read position i, so last
/// row's cell on the same diagonal is the same lane and the cell above is the next lane. Zero is -128, so
/// saturation is the local score's floor. Gaps along the model are a prefix maximum, taken in doubling steps.
#[cfg(any(target_arch = "aarch64", target_arch = "x86_64"))]
#[inline(always)]
unsafe fn bound<V: Bytes>(jobs: [Job; PAIRS], scale: f32, merge: bool) -> [i8; PAIRS] {
    let half = V::LANES as i64 / 2;
    let neg = V::splat(i8::MIN);
    let ext = (GAP_EXTEND * scale) as i8;
    let open = V::splat((GAP_OPEN * scale) as i8);
    let (e1, e2, e4, e8) = (V::splat(ext), V::splat(2 * ext), V::splat(4 * ext), V::splat(8 * ext));
    let mut last = [i8::MIN; 32];
    if merge {
        last[V::LANES - 1] = i8::MAX;
    }
    let last = V::load(last.as_ptr());
    let (mut h, mut del, mut best) = ([neg; PAIRS], [neg; PAIRS], [neg; PAIRS]);
    // each pair runs from the first read position whose band touches its model to the last
    let rows = jobs.map(|(prof, f, d)| {
        let cols = (prof.len() / 21 - 2 * PAD) as i64;
        ((-d - half).max(0), (f.len() as i64).min(cols - d + half))
    });
    let steps = rows.iter().map(|&(from, to)| to - from).max().unwrap().max(0);
    for step in 0..steps {
        for p in 0..PAIRS {
            let (prof, f, d) = jobs[p];
            let i = rows[p].0 + step;
            // a pair that is done reads padding, the first bytes of its profile
            let at = if i < rows[p].1 {
                *f.get_unchecked(i as usize) as usize * (prof.len() / 21) + (PAD as i64 + i + d - half) as usize
            } else {
                0
            };
            debug_assert!(at + V::LANES + 1 <= prof.len());
            let s = V::load(prof.as_ptr().add(at)).max(V::load(prof.as_ptr().add(at + 1)).min(last));
            // up one lane, except that the last keeps its own value too: it stands for two columns, and a
            // gap from the second to the first stays in it
            let up = h[p].subs(open).max(del[p].subs(e1));
            del[p] = up.lower(if merge { up.last_first() } else { neg });
            let m = h[p].adds(s).max(del[p]);
            let mut x = m.subs(open).raise::<1>(neg);
            x = x.max(x.raise::<1>(neg).subs(e1));
            x = x.max(x.raise::<2>(neg).subs(e2));
            x = x.max(x.raise::<4>(neg).subs(e4));
            x = x.max(x.raise::<8>(neg).subs(e8));
            if V::LANES > 16 {
                // in two, as 16 gap extensions don't fit in a byte
                x = x.max(x.raise::<16>(neg).subs(e8).subs(e8));
            }
            h[p] = m.max(x);
            best[p] = best[p].max(h[p]);
        }
    }
    best.map(|b| b.hmax())
}

fn search(pf: &Pfam, idx: &Index, path: &str, bar: f32, out: &str, seeds: Option<&String>) {
    let (orfs, k) = (read_orfs(path), idx.shape.kmer());
    let (hits, passed, align_ns) = (AtomicU64::new(0), AtomicU64::new(0), AtomicU64::new(0));
    let t0 = Instant::now();
    let chunks: Vec<(String, String)> = orfs
        .par_chunks(64)
        .map(|chunk| {
            let (mut text, mut seed_text, mut h, mut del, mut found) = (String::new(), String::new(), vec![], vec![], Vec::<(u32, f32)>::new());
            let (mut n_hits, mut n_passed, mut ns) = (0u64, 0u64, 0u64);
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
                enumerate(&pf.lo, &(a + c..a + c + k).collect::<Vec<_>>(), t, &keep, &mut |code, s| {
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
    let k = idx.shape.kmer();
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
const STRIDE: usize = 32; // bytes per profile column: 21 residues, X the last, then padding
const OUTSIDE: u8 = 21; // a window's residue outside its ORF, in the column's padding

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
    shape: Shape,
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
    f.write_all(bytes_of(&[idx.shape.header(), tiles.len() as u64])).unwrap();
    for t in &tiles {
        f.write_all(bytes_of(&[t.m0 as u64, t.m1 as u64, t.rows.len() as u64])).unwrap();
        f.write_all(bytes_of(&t.starts)).unwrap();
        f.write_all(bytes_of(&t.rows)).unwrap();
    }
    f.flush().unwrap();
    eprintln!("{} tiles, {} entries", tiles.len(), tiles.iter().map(|t| t.rows.len()).sum::<usize>());
}

fn load_tiles(path: &str) -> TileIndex {
    let mut f = std::io::BufReader::new(fs::File::open(path).unwrap());
    let mut h = [0u64; 2];
    f.read_exact(bytes_of_mut(&mut h)).unwrap();
    let shape = Shape::from_header(h[0]);
    let space = 20usize.pow(shape.k() as u32);
    let tiles = (0..h[1]).map(|_| {
        let mut t = [0u64; 3];
        f.read_exact(bytes_of_mut(&mut t)).unwrap();
        let mut starts = vec![0u32; space + 1];
        f.read_exact(bytes_of_mut(&mut starts)).unwrap();
        let mut rows = vec![0u16; t[2] as usize];
        f.read_exact(bytes_of_mut(&mut rows)).unwrap();
        Tile { m0: t[0] as usize, m1: t[1] as usize, starts, rows }
    }).collect();
    TileIndex { shape, tiles }
}

/// One block's postings, LANES at a time: per k-mer, groups of LANES windows stored offset-major
/// (offset o of lane l at o * LANES + l), WIN residues each, OUTSIDE beyond the ORF, all plus BIAS. Unused
/// lanes hold OUTSIDE everywhere, so they score nothing, and orf u32::MAX.
///
/// LANES is the widest register the architecture's kernel fills with one job: 32 bytes on AVX2, where the
/// 16-byte halves of a profile column are broadcast, so nothing is assembled from two jobs. Wider groups
/// are fuller only in larger blocks (hit-weighted, 0.69 full in 300,000 ORFs and 0.93 in 1.6 M).
#[cfg(target_arch = "x86_64")]
const LANES: usize = 32;
#[cfg(not(target_arch = "x86_64"))]
const LANES: usize = 16;
/// Added to every residue in a window, so that `pshufb`, which gives 0 for an index with its top bit set,
/// looks up residues 0..16 as they are (0x70..0x80) and 16..22 with the top bit flipped (0x80..0x86 to
/// 0..6): the two halves of a column are then joined with an OR.
const BIAS: u8 = if cfg!(target_arch = "x86_64") { 0x70 } else { 0 };
/// Bits per lane in a mask, as each architecture's SIMD compare gives them.
const LANE_BITS: usize = if cfg!(target_arch = "x86_64") { 1 } else { 4 };
const LANE: u64 = (1 << LANE_BITS) - 1;
struct Postings {
    gstarts: Vec<u32>, // per k-mer, into groups
    count: Vec<u32>,   // per k-mer, real postings
    windows: Vec<u8>,
    orf: Vec<u32>,
    pos: Vec<u32>,
    covered: Vec<u64>, // per group, a mask of the lanes whose window holds the whole ORF
    span: Vec<(u8, u8)>, // per group, the window offsets that hold a residue in some lane
}

/// The word at j is the residues `shape` takes from j on, and sits at the same place in the window as a k-mer would.
fn postings(orfs: &[&[u8]], shape: Shape) -> Postings {
    let (k, span, offsets) = (shape.k(), shape.span(), shape.offsets());
    let space = 20usize.pow(k as u32);
    let code_at = |f: &[u8], j: usize| -> Option<usize> {
        let mut code = 0usize;
        for &o in &offsets {
            let c = f[j + o];
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
    let per_orf: Vec<usize> = orfs.par_iter().map(|f| (0..(f.len() + 1).saturating_sub(span)).filter(|&j| code_at(f, j).is_some()).count()).collect();
    let mut at = Vec::with_capacity(orfs.len() + 1);
    at.push(0usize);
    for n in &per_orf {
        at.push(at.last().unwrap() + n);
    }
    let mut all = vec![(0u32, 0u8, 0u8, 0u32, 0u32); *at.last().unwrap()];
    let mut rest: &mut [(u32, u8, u8, u32, u32)] = &mut all;
    let mut parts = Vec::with_capacity(orfs.len());
    for &n in &per_orf {
        let (head, tail) = std::mem::take(&mut rest).split_at_mut(n);
        parts.push(head);
        rest = tail;
    }
    parts.into_par_iter().enumerate().for_each(|(o, out)| {
        let f = orfs[o];
        let mut i = 0;
        for j in 0..(f.len() + 1).saturating_sub(span) {
            if let Some(c) = code_at(f, j) {
                out[i] = (c as u32, HALF.saturating_sub(j) as u8, (HALF + f.len() - j).min(WIN) as u8, o as u32, j as u32);
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
    let (mut windows, mut orf, mut pos) = (vec![OUTSIDE + BIAS; groups * WIN * LANES], vec![u32::MAX; groups * LANES], vec![0u32; groups * LANES]);
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
                    w[x * LANES + lane] = f[j + x - HALF] + BIAS;
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
///
/// Nearly every batch is four rows against one group, so each SIMD kernel has a version (ONE) that loads
/// a window offset once for all four jobs.
type Masks = [(u64, u64); 4];
type Scorer4 = fn(&[i8], [usize; 4], [&[u8]; 4], (usize, usize), (i8, i8), &mut Masks);

fn one_group(w: &[&[u8]; 4]) -> bool {
    w.iter().all(|x| std::ptr::eq(x.as_ptr(), w[0].as_ptr()))
}

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
    if one_group(&w) {
        score4_neon_inner::<true>(prof, rows, w, span, gates, masks)
    } else {
        score4_neon_inner::<false>(prof, rows, w, span, gates, masks)
    }
}

#[cfg(target_arch = "aarch64")]
#[inline(always)]
fn score4_neon_inner<const ONE: bool>(prof: &[i8], rows: [usize; 4], w: [&[u8]; 4], span: (usize, usize), gates: (i8, i8), masks: &mut Masks) {
    use std::arch::aarch64::*;
    for i in 0..4 {
        assert!((rows[i] + WIN) * STRIDE <= prof.len() && w[i].len() >= WIN * LANES);
    }
    unsafe {
        let (mut cur, mut best) = ([vdupq_n_s8(i8::MIN); 4], [vdupq_n_s8(i8::MIN); 4]);
        let p = rows.map(|r| prof.as_ptr().add(r * STRIDE));
        for o in span.0..span.1 {
            let mut idx = vld1q_u8(w[0].as_ptr().add(o * LANES));
            for i in 0..4 {
                let table = vld1q_s8_x2(p[i].add(o * STRIDE));
                if !ONE {
                    idx = vld1q_u8(w[i].as_ptr().add(o * LANES));
                }
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
unsafe fn score4_sse41_inner<const ONE: bool>(prof: &[i8], rows: [usize; 4], w: [&[u8]; 4], span: (usize, usize), gates: (i8, i8), masks: &mut Masks) {
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
            let mut idx = _mm_loadu_si128(w[0].as_ptr().add(o * LANES + half * 16) as *const __m128i);
            let mut idx_hi = _mm_xor_si128(idx, flip);
            for i in 0..4 {
                let col = p[i].add(o * STRIDE) as *const __m128i;
                if !ONE {
                    idx = _mm_loadu_si128(w[i].as_ptr().add(o * LANES + half * 16) as *const __m128i);
                    idx_hi = _mm_xor_si128(idx, flip);
                }
                let lo = _mm_shuffle_epi8(_mm_loadu_si128(col), idx);
                let hi = _mm_shuffle_epi8(_mm_loadu_si128(col.add(1)), idx_hi);
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
    if one_group(&w) {
        unsafe { score4_sse41_inner::<true>(prof, rows, w, span, gates, masks) }
    } else {
        unsafe { score4_sse41_inner::<false>(prof, rows, w, span, gates, masks) }
    }
}

/// AVX2: one job's 32 lanes per register. `vpshufb` shuffles each 128-bit half with its own table, so both
/// halves get the same one, broadcast from memory, which costs a load and no shuffle.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn score4_avx2_inner<const ONE: bool>(prof: &[i8], rows: [usize; 4], w: [&[u8]; 4], span: (usize, usize), gates: (i8, i8), masks: &mut Masks) {
    use std::arch::x86_64::*;
    for i in 0..4 {
        assert!((rows[i] + WIN) * STRIDE <= prof.len() && w[i].len() >= WIN * LANES);
    }
    let p = rows.map(|r| prof.as_ptr().add(r * STRIDE));
    let flip = _mm256_set1_epi8(i8::MIN);
    let (mut cur, mut best) = ([flip; 4], [flip; 4]);
    for o in span.0..span.1 {
        let mut idx = _mm256_loadu_si256(w[0].as_ptr().add(o * LANES) as *const __m256i);
        let mut idx_hi = _mm256_xor_si256(idx, flip);
        for i in 0..4 {
            let col = p[i].add(o * STRIDE) as *const __m128i;
            if !ONE {
                idx = _mm256_loadu_si256(w[i].as_ptr().add(o * LANES) as *const __m256i);
                idx_hi = _mm256_xor_si256(idx, flip);
            }
            let lo = _mm256_shuffle_epi8(_mm256_broadcastsi128_si256(_mm_loadu_si128(col)), idx);
            let hi = _mm256_shuffle_epi8(_mm256_broadcastsi128_si256(_mm_loadu_si128(col.add(1))), idx_hi);
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
    if one_group(&w) {
        unsafe { score4_avx2_inner::<true>(prof, rows, w, span, gates, masks) }
    } else {
        unsafe { score4_avx2_inner::<false>(prof, rows, w, span, gates, masks) }
    }
}

/// AVX-512 VBMI: `vpermb` is a lookup in the whole 32-byte column, its table straight from memory. It uses
/// an index's low 5 bits, so the BIAS of 0x70 (16 modulo 32) is taken off with one bit. Four instructions a
/// job with AVX2's lanes: under half as many.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx512f,avx512bw,avx512vl,avx512vbmi")]
unsafe fn score4_avx512_inner<const ONE: bool>(prof: &[i8], rows: [usize; 4], w: [&[u8]; 4], span: (usize, usize), gates: (i8, i8), masks: &mut Masks) {
    use std::arch::x86_64::*;
    for i in 0..4 {
        assert!((rows[i] + WIN) * STRIDE <= prof.len() && w[i].len() >= WIN * LANES);
    }
    let p = rows.map(|r| prof.as_ptr().add(r * STRIDE));
    let (zero, unbias) = (_mm256_set1_epi8(i8::MIN), _mm256_set1_epi8(16));
    let (mut cur, mut best) = ([zero; 4], [zero; 4]);
    for o in span.0..span.1 {
        let mut idx = _mm256_xor_si256(_mm256_loadu_si256(w[0].as_ptr().add(o * LANES) as *const __m256i), unbias);
        for i in 0..4 {
            if !ONE {
                idx = _mm256_xor_si256(_mm256_loadu_si256(w[i].as_ptr().add(o * LANES) as *const __m256i), unbias);
            }
            let s = _mm256_permutexvar_epi8(idx, _mm256_loadu_si256(p[i].add(o * STRIDE) as *const __m256i));
            cur[i] = _mm256_adds_epi8(cur[i], s);
            best[i] = _mm256_max_epi8(best[i], cur[i]);
        }
    }
    let (g0, g1) = (_mm256_set1_epi8(gates.0 ^ i8::MIN), _mm256_set1_epi8(gates.1 ^ i8::MIN));
    for i in 0..4 {
        masks[i] = (_mm256_cmpge_epi8_mask(best[i], g0) as u64, _mm256_cmpge_epi8_mask(best[i], g1) as u64);
    }
}

#[cfg(target_arch = "x86_64")]
fn score4_avx512(prof: &[i8], rows: [usize; 4], w: [&[u8]; 4], span: (usize, usize), gates: (i8, i8), masks: &mut Masks) {
    if one_group(&w) {
        unsafe { score4_avx512_inner::<true>(prof, rows, w, span, gates, masks) }
    } else {
        unsafe { score4_avx512_inner::<false>(prof, rows, w, span, gates, masks) }
    }
}

/// Every kernel this CPU has, fastest first. FILTER_SIMD=name forces one.
fn scorers() -> Vec<(&'static str, Scorer4)> {
    let mut all: Vec<(&'static str, Scorer4)> = vec![];
    #[cfg(target_arch = "aarch64")]
    all.push(("neon", score4_neon));
    #[cfg(target_arch = "x86_64")]
    {
        if is_x86_feature_detected!("avx512vbmi") && is_x86_feature_detected!("avx512vl") && is_x86_feature_detected!("avx512bw") {
            all.push(("avx512", score4_avx512));
        }
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
/// and WIN - HALF after; padding and OUTSIDE score i8::MIN, which resets a local score.
fn int8_profile(pf: &Pfam, m: usize, out: &mut Vec<i8>) {
    let (a, b) = (pf.offs[m] as usize, pf.offs[m + 1] as usize);
    let at = out.len();
    out.resize(at + (b - a + WIN) * STRIDE, i8::MIN);
    let out = &mut out[at..];
    for c in 0..b - a {
        for r in 0..21 {
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

// A survivor is model << 48 | ORF << 24 | diagonal key, so it sorts as one integer. The key is column less read position,
// offset to stay positive.
const KEY_BITS: u32 = 24;
const KEY_OFFSET: u32 = 1 << (KEY_BITS - 1); // longest sequence, in residues
const MAX_BLOCK: usize = 1 << (48 - KEY_BITS); // ORFs
fn pack(model: u32, orf: u32, column: u32, pos: u32) -> u64 {
    (model as u64) << 48 | (orf as u64) << KEY_BITS | (column + KEY_OFFSET - pos) as u64
}
/// Model, ORF, column less position.
fn unpack(x: u64) -> (u32, u32, i64) {
    ((x >> 48) as u32, (x >> KEY_BITS) as u32 & (MAX_BLOCK as u32 - 1), (x & (2 * KEY_OFFSET as u64 - 1)) as i64 - KEY_OFFSET as i64)
}

/// FILTER_PER_ORF's checks sit in the join's inner loops, so the join is compiled with and without them.
fn join(pf: &Pfam, ti: &TileIndex, path: &str, bar: f32, gate_bits: f32, block: usize, out: &str) {
    match std::env::var("FILTER_PER_ORF") {
        Ok(p) => join_with::<true>(pf, ti, path, bar, gate_bits, block, out, Some(p)),
        Err(_) => join_with::<false>(pf, ti, path, bar, gate_bits, block, out, None),
    }
}

#[allow(clippy::too_many_arguments)]
fn join_with<const PER_ORF: bool>(pf: &Pfam, ti: &TileIndex, path: &str, bar: f32, gate_bits: f32, block: usize, out: &str,
    per_orf: Option<String>) {
    use std::io::BufRead;
    let mut lines = std::io::BufReader::with_capacity(1 << 20, fs::File::open(path).unwrap()).lines();
    let mut n_orfs = 0usize;
    let gate = (gate_bits * 4.0).round() as i8;
    assert!(pf.models() < 1 << 16 && block <= MAX_BLOCK, "a block is at most {MAX_BLOCK} ORFs");
    // gapped score in a band of +-8 around the best diagonal (step 4); FILTER_BAND=b or =whole to change it. Up to
    // 16 has a SIMD bound ahead of it
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
    // FILTER_DIAGONALS=all: score a pair in a band around every diagonal at the bar, not only its best one. any: the same
    // pairs, but a pair is done at the first diagonal that reaches its bar, so the scores it is written with are that one's
    let diagonals = std::env::var("FILTER_DIAGONALS").unwrap_or_default();
    let (all_diagonals, first_only) = (matches!(&diagonals[..], "all" | "any"), diagonals == "any");
    let covered_gate = (((bar - SLACK) * 4.0).floor().min(127.0) as i8).max(gate);
    // FILTER_MIN_GAPPED=g: write only pairs whose gapped score is at least g (HMMER's bar), so pairs scored
    // at the diagonal bar but below the hand-off are not formatted and written (324 M rows at diagonal 17).
    // Half a printed digit lower, so a score that prints as g is kept and readers' cut on the text decides.
    let min_gapped: f32 = std::env::var("FILTER_MIN_GAPPED").map_or(f32::NEG_INFINITY, |v| v.parse::<f32>().unwrap() - 0.005);
    // FILTER_GAPPED_FILE: a bar per model, one per line, where that is lower. A short repeat at GA scores less
    // than any one bar for the rest: several copies of it add up to GA1
    let min_gapped: Vec<f32> = match std::env::var("FILTER_GAPPED_FILE") {
        Ok(p) => fs::read_to_string(p).unwrap().lines().map(|l| (l.trim().parse::<f32>().unwrap() - 0.005).min(min_gapped)).collect(),
        Err(_) => vec![min_gapped; pf.models()],
    };
    assert_eq!(min_gapped.len(), pf.models(), "FILTER_GAPPED_FILE must hold a line per model");
    let bound = band.and_then(|b| Bound::new(pf, &min_gapped, b));
    let bounded = AtomicU64::new(0); // pairs the bound lets through to the f32 score
    let (engine, score4) = scorer4();
    let (hits, survivors, pairs_n) = (AtomicU64::new(0), AtomicU64::new(0), AtomicU64::new(0));
    let (surv_ns, unique, gap_ns) = (AtomicU64::new(0), AtomicU64::new(0), AtomicU64::new(0)); // surv: sort, dedupe and exact diagonals
    let (mut t_post, mut t_join) = (0.0f64, 0.0f64);
    // FILTER_BINARY=1: per block its ORFs, residues and bytes of pairs (u64), then the pairs, (ORF in the block, model) as
    // two u32 in model order, little endian. For a reader that streams the same FASTA beside it: 8 bytes a pair, not 60
    let binary = std::env::var("FILTER_BINARY").is_ok_and(|v| v == "1");
    // FILTER_PER_ORF=path: a diagnostic, each ORF's work (see PerOrf). Off, nothing is counted or timed
    let mut per_orf = per_orf.map(|p| {
        let mut f = std::io::BufWriter::new(fs::File::create(p).unwrap());
        use std::io::Write;
        writeln!(f, "orf\tname\tlength\thits\tdiagonals\tdiagonal_ns\tbounded\tbound_ns\tgapped\tgapped_ns\twritten").unwrap();
        f
    });
    let timer_ns = if per_orf.is_some() { timer_cost() } else { 0.0 };
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
        assert!(chunk.iter().all(|f| f.len() <= KEY_OFFSET as usize), "a sequence is at most {KEY_OFFSET} residues");
        let t = Instant::now();
        let p = postings(&chunk, ti.shape);
        t_post += t.elapsed().as_secs_f64();
        let real: u64 = p.count.iter().map(|&c| c as u64).sum();
        eprintln!("block of {} ORFs: {} postings in {} groups, lanes filled {:.3}, {:.0} MB", chunk.len(), real,
            p.gstarts.last().unwrap(), real as f64 / (*p.gstarts.last().unwrap() as f64 * LANES as f64), p.windows.len() as f64 / 1e6);
        let t = Instant::now();
        // each tile's pairs are scored and formatted before the next tile's are found, so a block holds its
        // text, not its pairs
        let on = PER_ORF;
        let rows: Vec<(Vec<u8>, Vec<PerOrf>)> = ti.tiles
            .par_iter()
            .map_init(|| (Vec::<i8>::new(), Vec::<u64>::new(), vec![], vec![]), |(prof, surv, h, del), tile| {
                let mut work = if on { vec![PerOrf::default(); chunk.len()] } else { vec![] };
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
                let mut compact_at = SURV_COMPACT;
                let mut drain = |jobs: &[(usize, usize)], masks: &Masks, surv: &mut Vec<u64>| {
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
                            surv.push(pack(m, p.orf[q], c as u32, p.pos[q]));
                        }
                    }
                    // a diagonal survives once per k-mer on it: compact before the list grows large. Again only when it
                    // has doubled, or a tile with more than SURV_COMPACT unique diagonals sorts them at every batch
                    if surv.len() > compact_at {
                        surv.sort_unstable();
                        surv.dedup();
                        compact_at = SURV_COMPACT.max(2 * surv.len());
                    }
                };
                for code in 0..tile.starts.len() - 1 {
                    let (r0, r1) = (tile.starts[code] as usize, tile.starts[code + 1] as usize);
                    let (g0, g1) = (p.gstarts[code] as usize, p.gstarts[code + 1] as usize);
                    if r0 == r1 || g0 == g1 || dropped[code] {
                        continue;
                    }
                    n_hits += (r1 - r0) as u64 * p.count[code] as u64;
                    // Four rows against every group in turn, so a profile row is read once per k-mer. Rows left
                    // over wait for company, whatever its k-mer.
                    let rows = &tile.rows[r0..r1];
                    let whole = rows.len() / 4 * 4;
                    for four in rows[..whole].chunks_exact(4) {
                        let four = [four[0] as usize, four[1] as usize, four[2] as usize, four[3] as usize];
                        for g in g0..g1 {
                            let w = &p.windows[g * WIN * LANES..(g + 1) * WIN * LANES];
                            score4(prof, four, [w; 4], (p.span[g].0 as usize, p.span[g].1 as usize), (gate, covered_gate), &mut masks);
                            // two batches in three have no survivor
                            if masks.iter().fold(0, |any, m| any | m.0) != 0 {
                                drain(&four.map(|row| (row, g)), &masks, surv);
                            }
                        }
                    }
                    for &row in &rows[whole..] {
                        for g in g0..g1 {
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
                    let (m, o, d) = unpack(x); // d: column - position
                    let f = chunk[o as usize];
                    let cols = (pf.offs[m as usize + 1] - pf.offs[m as usize]) as usize;
                    // a point on the diagonal inside both the read and the model
                    let j = if d >= 0 { 0usize } else { (-d) as usize };
                    let col = (j as i64 + d) as usize;
                    if j >= f.len() || col >= cols {
                        continue;
                    }
                    let td = on.then(Instant::now);
                    let s = diagonal(pf, m as usize, col, f, j);
                    if let Some(td) = td {
                        let w = &mut work[o as usize];
                        w.diagonals += 1;
                        w.diagonal_ns += (td.elapsed().as_nanos() as f64 - timer_ns).max(0.0);
                    }
                    if s >= bar {
                        match res.last_mut() {
                            Some(x) if !all_diagonals && x.0 == o && x.1 == m => {
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
                use std::io::Write;
                let (mut text, mut n) = (Vec::<u8>::new(), 0u64);
                // a pair's diagonals come one after the other, and it is written once, with its best scores
                let mut pending: Option<(u32, u32, f32, f32)> = None;
                let write = |text: &mut Vec<u8>, pair: Option<(u32, u32, f32, f32)>| {
                    if let Some((o, m, d, g)) = pair {
                        if binary {
                            text.extend_from_slice(&o.to_le_bytes());
                            text.extend_from_slice(&m.to_le_bytes());
                            return;
                        }
                        writeln!(text, "{}\t{m}\t{d:.2}\t{g:.2}", orfs.name(o as usize)).unwrap();
                    }
                };
                for four in res.chunks(PAIRS) {
                    let tb = on.then(Instant::now);
                    let reach = bound.as_ref().map_or([true; PAIRS], |b| {
                        // a short last chunk repeats its first pair
                        b.may_reach(std::array::from_fn(|i| {
                            let &(o, m, _, diag) = four.get(i).unwrap_or(&four[0]);
                            (m as usize, chunk[o as usize], diag)
                        }))
                    });
                    if let Some(tb) = tb {
                        let share = (tb.elapsed().as_nanos() as f64 - timer_ns).max(0.0) / four.len() as f64;
                        for &(o, ..) in four {
                            work[o as usize].bounded += 1;
                            work[o as usize].bound_ns += share;
                        }
                    }
                    for (&(o, m, d, diag), _) in four.iter().zip(reach).filter(|(_, r)| *r) {
                        if first_only && pending.is_some_and(|p| p.0 == o && p.1 == m) {
                            continue;
                        }
                        let f = chunk[o as usize];
                        n += 1;
                        let a = pf.offs[m as usize] as usize;
                        let cols = pf.offs[m as usize + 1] as usize - a;
                        let te = on.then(Instant::now);
                        let g = match band {
                            Some(b) => gapped_band(&pf.lo, a, cols, f, diag, b, h, del),
                            None => gapped(&pf.lo, a, cols, f, h, del),
                        };
                        if let Some(te) = te {
                            work[o as usize].gapped += 1;
                            work[o as usize].gapped_ns += (te.elapsed().as_nanos() as f64 - timer_ns).max(0.0);
                        }
                        if g >= min_gapped[m as usize] {
                            if on && !pending.is_some_and(|p| p.0 == o && p.1 == m) {
                                work[o as usize].written += 1;
                            }
                            match &mut pending {
                                Some(p) if p.0 == o && p.1 == m => (p.2, p.3) = (p.2.max(d), p.3.max(g)),
                                _ => write(&mut text, pending.replace((o, m, d, g))),
                            }
                        }
                    }
                }
                write(&mut text, pending);
                bounded.fetch_add(n, Relaxed);
                pairs_n.fetch_add(res.len() as u64, Relaxed);
                gap_ns.fetch_add(ts.elapsed().as_nanos() as u64, Relaxed);
                (text, work)
            })
            .collect();
        t_join += t.elapsed().as_secs_f64();
        if let Some(f) = per_orf.as_mut() {
            write_per_orf(f, &orfs, &chunk, &p, ti, &dropped, &rows, n_orfs - orfs.len());
        }
        if binary {
            for x in [orfs.len(), orfs.seqs.len(), rows.iter().map(|r| r.0.len()).sum()] {
                sink.write_all(&(x as u64).to_le_bytes()).unwrap();
            }
        }
        for (r, _) in rows {
            sink.write_all(&r).unwrap();
        }
        if binary {
            // the reader works on a block while the next is searched
            sink.flush().unwrap();
        }
    }
    sink.flush().unwrap();
    if let Some(mut f) = per_orf {
        f.flush().unwrap();
    }
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

/// One ORF's work in `join`, for FILTER_PER_ORF. Counts: diagonals scored whole (the gate's survivors, deduplicated),
/// pairs through the SIMD bound and through the exact gapped score, pairs written. Times in core-ns, less the timer's
/// own cost; the bound scores PAIRS pairs at once, and each gets an equal share. The gate and the survivors' sort are
/// not timed, as that would slow the gate's loop even when this is off: their cost per k-mer hit hardly varies, so a
/// reader shares the log's gate time out by `hits`, and the sort's by `diagonals`.
#[derive(Clone, Default)]
struct PerOrf {
    diagonals: u64,
    diagonal_ns: f64,
    bounded: u64,
    bound_ns: f64,
    gapped: u64,
    gapped_ns: f64,
    written: u64,
}

/// The cost of one Instant::now and elapsed pair, in ns.
fn timer_cost() -> f64 {
    let n = 100_000;
    let t = Instant::now();
    let mut x = 0u128;
    for _ in 0..n {
        x += std::hint::black_box(Instant::now()).elapsed().as_nanos();
    }
    std::hint::black_box(x);
    t.elapsed().as_nanos() as f64 / n as f64
}

/// A block's rows of FILTER_PER_ORF: the tiles' work summed, and each ORF's k-mer hits, as the join counts them.
#[allow(clippy::too_many_arguments)]
fn write_per_orf(f: &mut impl std::io::Write, orfs: &Packed, chunk: &[&[u8]], p: &Postings, ti: &TileIndex, dropped: &[bool],
    rows: &[(Vec<u8>, Vec<PerOrf>)], first: usize) {
    let mut hits = vec![0u64; chunk.len()];
    for code in 0..ti.tiles[0].starts.len() - 1 {
        let n: u64 = ti.tiles.iter().map(|t| (t.starts[code + 1] - t.starts[code]) as u64).sum();
        if n == 0 || dropped[code] {
            continue;
        }
        let (g0, count) = (p.gstarts[code] as usize, p.count[code] as usize);
        for q in g0 * LANES..g0 * LANES + count {
            hits[p.orf[q] as usize] += n;
        }
    }
    for i in 0..chunk.len() {
        let mut w = PerOrf::default();
        for (_, t) in rows {
            let x = &t[i];
            w.diagonals += x.diagonals;
            w.diagonal_ns += x.diagonal_ns;
            w.bounded += x.bounded;
            w.bound_ns += x.bound_ns;
            w.gapped += x.gapped;
            w.gapped_ns += x.gapped_ns;
            w.written += x.written;
        }
        writeln!(f, "{}\t{}\t{}\t{}\t{}\t{:.0}\t{}\t{:.0}\t{}\t{:.0}\t{}", first + i, orfs.name(i), chunk[i].len(), hits[i],
            w.diagonals, w.diagonal_ns, w.bounded, w.bound_ns, w.gapped, w.gapped_ns, w.written).unwrap();
    }
}

// What a seed would do at each threshold, model by model, from the HMM alone.
const READ: usize = 50; // residues in a read of 150 nt
const THRESHOLDS: usize = 25; // steps of a quarter bit from lowest_t
const Q: [f64; 20] = [0.0787945, 0.0151600, 0.0535222, 0.0668298, 0.0397062, 0.0695071, 0.0229198, 0.0590092, 0.0594422, 0.0963728,
    0.0237718, 0.0414386, 0.0482904, 0.0395639, 0.0540978, 0.0683364, 0.0540687, 0.0673417, 0.0114135, 0.0304133]; // prep.py's background
fn lowest_t(shape: Shape) -> f32 { if shape.k() == 3 { 3.0 } else { 4.5 } }
const PER_BIT: f32 = 8.0; // bins of a score histogram
const LOW: f32 = -6.0; // a residue below this is in no word worth having
const ONE: usize = 128; // bins of one column's scores, LOW..10
const TWO: usize = 2 * ONE; // of two columns' sum

/// How each column's residues are drawn in a read that just reaches GA: p(r) ~ Q(r) 2^(beta score(r)), beta set for
/// each READ-column stretch so that its mean score is `need`. Gives every column's beta, by the stretch it is
/// the middle of, and each stretch's weight: one that needs a stronger tilt holds fewer of the reads that pass.
fn tilts(lo: &[f32], cols: usize, need: f64) -> (Vec<f64>, Vec<f64>) {
    let betas: Vec<f64> = (0..20).map(|i| 0.2 + 3.8 * i as f64 / 19.0).collect();
    let n = READ.min(cols);
    let mut beta = vec![*betas.last().unwrap(); cols - n + 1];
    let mut done = vec![false; cols - n + 1];
    let mut before = vec![f64::NAN; cols - n + 1];
    for (b, &x) in betas.iter().enumerate() {
        let mean: Vec<f64> = (0..cols).map(|c| {
            let (mut z, mut e) = (0.0, 0.0);
            for r in 0..20 {
                let l = lo[c * 21 + r] as f64;
                let w = Q[r] * (x * l).exp2();
                z += w;
                e += w * l;
            }
            e / z
        }).collect();
        let mut sum: f64 = mean[..n].iter().sum();
        for w in 0..beta.len() {
            if w > 0 {
                sum += mean[w + n - 1] - mean[w - 1];
            }
            // a stretch's mean score rises with beta
            if !done[w] && sum >= need {
                beta[w] = if b == 0 { x } else { betas[b - 1] + (x - betas[b - 1]) * (need - before[w]) / (sum - before[w]) };
                done[w] = true;
            }
            before[w] = sum;
        }
    }
    let least = beta.iter().copied().fold(f64::MAX, f64::min);
    let weight = beta.iter().map(|b| (-4.0 * (b - least)).exp2()).collect();
    let column = (0..cols).map(|c| beta[c.saturating_sub(n / 2).min(beta.len() - 1)]).collect();
    (column, weight)
}

/// Per shape and T: background hits per residue, and -ln P(a read at GA holds no word), words taken as independent.
/// A word of four is a pair of residues and another pair, a word of three a pair and one, so every shape with a gap in common
/// shares that gap's pair histograms, and a word's score distribution is one histogram against the other's tail.
fn seed_scores(lo: &[f32], cols: usize, need: f64, shapes: &[Shape]) -> Vec<f32> {
    let (beta, weight) = tilts(lo, cols, need);
    let mut p = vec![[0.0f32; 20]; 2 * cols]; // background, then tilted, per column
    let mut bin = vec![[None::<usize>; 20]; cols];
    for c in 0..cols {
        let w: Vec<f64> = (0..20).map(|r| Q[r] * (beta[c] * lo[c * 21 + r] as f64).exp2()).collect();
        let z: f64 = w.iter().sum();
        let q: f64 = Q.iter().sum();
        for r in 0..20 {
            p[2 * c][r] = (Q[r] / q) as f32;
            p[2 * c + 1][r] = (w[r] / z) as f32;
            let b = ((lo[c * 21 + r] - LOW) * PER_BIT).round();
            bin[c][r] = (b >= 0.0).then(|| (b as usize).min(ONE - 1));
        }
    }
    // histograms of a pair `gap` apart from column c, and one column's: [2 * c + e] is a run of bins
    let pairs = |gap: usize| -> Vec<f32> {
        let mut h = vec![0.0f32; cols.saturating_sub(gap) * 2 * TWO];
        for c in 0..cols.saturating_sub(gap) {
            for x in 0..20 {
                let Some(bx) = bin[c][x] else { continue };
                for y in 0..20 {
                    let Some(by) = bin[c + gap][y] else { continue };
                    for e in 0..2 {
                        h[(2 * c + e) * TWO + bx + by] += p[2 * c + e][x] * p[2 * (c + gap) + e][y];
                    }
                }
            }
        }
        h
    };
    let tails = |mut h: Vec<f32>, bins: usize| -> Vec<f32> {
        for run in h.chunks_mut(bins) {
            for b in (0..bins - 1).rev() {
                run[b] += run[b + 1];
            }
        }
        h
    };
    let mut ones = vec![0.0f32; 2 * cols * ONE];
    for c in 0..cols {
        for r in 0..20 {
            if let Some(b) = bin[c][r] {
                ones[2 * c * ONE + b] += p[2 * c][r];
                ones[(2 * c + 1) * ONE + b] += p[2 * c + 1][r];
            }
        }
    }
    let ones = tails(ones, ONE);
    let (mut first, mut second): (Vec<Option<Vec<f32>>>, Vec<Option<Vec<f32>>>) = (vec![None; 16], vec![None; 16]);
    let mut out = vec![0.0f32; shapes.len() * 2 * THRESHOLDS];
    for (shape, out) in shapes.iter().zip(out.chunks_mut(2 * THRESHOLDS)) {
        let at = shape.offsets();
        if cols < shape.span() {
            continue;
        }
        let gap = at[1];
        if first[gap].is_none() {
            first[gap] = Some(pairs(gap));
        }
        let (tail, bins, tail_at) = if at.len() == 4 {
            let gap = at[3] - at[2];
            if second[gap].is_none() {
                second[gap] = Some(tails(pairs(gap), TWO));
            }
            (second[gap].as_ref().unwrap(), TWO, at[2])
        } else {
            (&ones, ONE, at[2])
        };
        let head = first[gap].as_ref().unwrap();
        let lowest = ((lowest_t(*shape) - at.len() as f32 * LOW) * PER_BIT).round() as usize;
        let starts = cols - shape.span() + 1;
        let mut expect = vec![0.0f64; starts * THRESHOLDS]; // words a start column is expected to give a read at GA
        for c in 0..starts {
            for e in 0..2 {
                let h = &head[(2 * c + e) * TWO..][..TWO];
                let t = &tail[(2 * (c + tail_at) + e) * bins..][..bins];
                let top = h.iter().rposition(|&x| x > 0.0).map_or(0, |b| b + 1);
                let reach = |at: usize| -> f32 { (at.saturating_sub(bins - 1)..top).map(|b| h[b] * t[at.saturating_sub(b)]).sum() };
                for i in 0..THRESHOLDS {
                    let at = lowest + i * (PER_BIT as usize / 4);
                    // a bin holds the scores that round to it, half of the first one's below T
                    let v = 0.5 * (reach(at) + reach(at + 1));
                    if e == 0 { out[i] += v } else { expect[c * THRESHOLDS + i] = v as f64 }
                }
            }
        }
        let width = (READ.min(cols) - shape.span().min(READ.min(cols)) + 1).min(starts);
        for i in 0..THRESHOLDS {
            let (mut run, mut none, mut total) = (0.0, 0.0, 0.0);
            for c in 0..starts {
                run += expect[c * THRESHOLDS + i];
                if c >= width {
                    run -= expect[(c - width) * THRESHOLDS + i];
                }
                if c + 1 >= width {
                    let w = weight.get(c + 1 - width).copied().unwrap_or(0.0);
                    total += w;
                    none += w * (-run).exp();
                }
            }
            out[THRESHOLDS + i] = -(none / total + 1e-12).ln() as f32;
        }
    }
    out
}

/// The table choose_seeds.py reads: f32, models x shapes x (hits, quality) x THRESHOLDS. GA: accession, GA1, GA2 per line.
fn seeds(pf: &Pfam, prep: &str, ga: &str, shapes: &[Shape], out: &str) {
    let t0 = Instant::now();
    let ga2: std::collections::HashMap<String, f64> = fs::read_to_string(ga).unwrap().lines().map(|l| {
        let mut c = l.split_whitespace();
        (c.next().unwrap().to_string(), c.nth(1).unwrap().parse().unwrap())
    }).collect();
    let accs = fs::read_to_string(format!("{prep}/accs.txt")).unwrap();
    let accs: Vec<&str> = accs.lines().collect();
    assert_eq!(accs.len(), pf.models());
    let rows: Vec<Vec<f32>> = (0..pf.models()).into_par_iter().map(|m| {
        let (a, b) = (pf.offs[m] as usize, pf.offs[m + 1] as usize);
        let cols = b - a;
        // raw bits a domain at GA2 has: HMMER's score, plus what entering the model and the ORF's length cost it
        let need = ga2[accs[m]] + (cols as f64 * (cols as f64 + 1.0) / 2.0).log2() + 2.0 * (READ as f64 / 2.0 + 1.5).log2() - (READ as f64 + 1.0).log2();
        seed_scores(&pf.lo[a * 21..b * 21], cols, need, shapes)
    }).collect();
    let mut f = std::io::BufWriter::new(fs::File::create(out).unwrap());
    for r in &rows {
        f.write_all(bytes_of(r)).unwrap();
    }
    f.flush().unwrap();
    eprintln!("{} models, {} shapes, {:.0} s", pf.models(), shapes.len(), t0.elapsed().as_secs_f32());
}

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let pf = Pfam::load(&a[2]);
    match a[1].as_str() {
        "build" => build(&pf, Shape::parse(&a[3]), a[4].parse().unwrap(), &a[5]),
        "search" => search(&pf, &load_index(&a[3]), &a[4], a[5].parse().unwrap(), &a[6], a.get(7)),
        "window" => window(&pf, &load_index(&a[3]), &a[4], a[5].parse().unwrap(), &a[6]),
        "tile" => tile(&pf, &load_index(&a[3]), &a[4]),
        "join" => join(&pf, &load_tiles(&a[3]), &a[7], a[4].parse().unwrap(), a[5].parse().unwrap(), a[6].parse().unwrap(), &a[8]),
        "seeds" => seeds(&pf, &a[2], &a[3], &a[5].split(',').map(Shape::parse).collect::<Vec<_>>(), &a[4]),
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
        enumerate(&lo, &[0, 1], 9.0, &|_, _| true, &mut |c, s| seen.push((c, s)));
        assert_eq!(seen, vec![(0, 10.0)]);
        seen.clear();
        enumerate(&lo, &[0, 1], 9.0, &|i, _| i < 2, &mut |c, s| seen.push((c, s)));
        assert!(seen.is_empty());
    }

    #[test]
    fn a_spaced_word_is_the_same_word_in_the_index_and_in_a_read() {
        for (text, k, span) in [("3", 3, 3), ("1011001", 4, 7), ("11001", 3, 5)] {
            let shape = Shape::parse(text);
            assert_eq!((shape.k(), shape.span(), Shape::from_header(shape.header())), (k, span, shape));
            // a model whose consensus is a read, with every other residue hopeless, holds exactly the read's words
            let f: Vec<u8> = (0..30).map(|i| (i * 7 % 20) as u8).collect();
            let lo: Vec<f32> = f.iter().flat_map(|&b| (0..21u8).map(move |r| if r == b { 4.0 } else { -9.0 })).collect();
            let words: Vec<usize> = (0..=f.len() - span)
                .map(|c| {
                    let mut seen = vec![];
                    enumerate(&lo, &shape.offsets().iter().map(|o| c + o).collect::<Vec<_>>(), 4.0 * k as f32, &|_, _| true, &mut |w, _| seen.push(w));
                    assert_eq!(seen.len(), 1);
                    seen[0]
                })
                .collect();
            let p = postings(&[&f], shape);
            let mut in_read = vec![usize::MAX; words.len()];
            for w in 0..p.count.len() {
                for g in p.gstarts[w] as usize..p.gstarts[w + 1] as usize {
                    for lane in (0..LANES).filter(|l| p.orf[g * LANES + l] == 0) {
                        in_read[p.pos[g * LANES + lane] as usize] = w;
                    }
                }
            }
            assert_eq!(in_read, words, "{text}");
        }
    }

    #[test]
    fn every_kernel_equals_scalar() {
        let mut x = 99u64;
        let mut next = || { x ^= x << 13; x ^= x >> 7; x ^= x << 17; x };
        // scores up to 60 a column, so some lanes saturate
        let prof: Vec<i8> = (0..3 * WIN * STRIDE).map(|i| if i % STRIDE > 20 { i8::MIN } else { (next() % 90) as i8 - 30 }).collect();
        for t in 0..400 {
            // residues and padding, or only well-scoring ones
            let ws: Vec<Vec<u8>> = (0..4).map(|_| (0..WIN * LANES).map(|_| (next() % 22) as u8 + BIAS).collect()).collect();
            // four groups, or one for all four jobs
            let w = [0, 1, 2, 3].map(|i| &ws[if t % 2 == 0 { i } else { 0 }][..]);
            let rows = [0, 1, 2, 3].map(|_| (next() % (2 * WIN as u64)) as usize);
            let (span, gates) = (((next() % 10) as usize, WIN - (next() % 10) as usize), ([66, 10, 127, 0][t / 2 % 4], [68, 127, 100, 40][t / 2 % 4]));
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
    fn an_unknown_residue_does_not_end_an_alignment() {
        // a 30-column model that scores 1 bit for one residue a column, and a read of it with an X at 15
        let n = 30;
        let mut lo: Vec<f32> = (0..n * 21).map(|i| if i % 21 == i / 21 % 20 { 1.0 } else { -0.5 }).collect();
        lo.chunks_mut(21).for_each(score_unknown);
        let pf = Pfam { lo, offs: vec![0, n as u32] };
        let mut f: Vec<u8> = (0..n as u8).map(|c| c % 20).collect();
        f[15] = 20;
        let want = 1.0 * (n - 1) as f32 + (Q[15] * 1.5 / Q.iter().sum::<f64>() - 0.5) as f32;
        assert!((diagonal(&pf, 0, 0, &f, 0) - want).abs() < 1e-3);
        assert!((gapped(&pf.lo, 0, n, &f, &mut vec![], &mut vec![]) - want).abs() < 1e-3);
        // the window of the word at 0 holds 26 residues and padding before them: in quarter bits 25 x 4 - 2 through
        // the X, and 15 x 4 if it broke the run
        let p = postings(&[&f], Shape::parse("4"));
        let mut prof = vec![];
        int8_profile(&pf, 0, &mut prof);
        let at = (0..p.orf.len()).find(|&l| p.orf[l] == 0 && p.pos[l] == 0).unwrap();
        let w = &p.windows[at / LANES * WIN * LANES..][..WIN * LANES];
        for (name, kernel) in scorers() {
            let mut got: Masks = [(0, 0); 4];
            kernel(&prof, [0; 4], [w; 4], (0, WIN), (98, 99), &mut got);
            let lane = LANE << (LANE_BITS * (at % LANES));
            assert!(got[0].0 & lane != 0 && got[0].1 & lane == 0, "{name}");
        }
    }

    #[test]
    fn seed_scores_are_what_enumerating_gives() {
        let mut x = 5u64;
        let mut next = || { x ^= x << 13; x ^= x >> 7; x ^= x << 17; x };
        let cols = 14;
        let lo: Vec<f32> = (0..cols * 21).map(|i| if i % 21 == 20 { -100.0 } else { (next() % 1000) as f32 / 1000.0 * 9.0 - 5.0 }).collect();
        let shapes = ["1101001", "1011", "11001", "111"].map(Shape::parse);
        let table = seed_scores(&lo, cols, 20.0, &shapes);
        let z: f64 = Q.iter().sum();
        for (shape, row) in shapes.iter().zip(table.chunks(2 * THRESHOLDS)) {
            for i in [0, 4, 8] {
                let t = lowest_t(*shape) + 0.25 * i as f32;
                let exact = |t: f32| -> f64 {
                    (0..=cols - shape.span()).map(|c| {
                        let mut p = 0.0;
                        let at: Vec<usize> = shape.offsets().iter().map(|o| c + o).collect();
                        enumerate(&lo, &at, t, &|_, _| true, &mut |code, _| {
                            let (mut w, mut code) = (1.0, code);
                            for _ in 0..at.len() {
                                w *= Q[code % 20] / z;
                                code /= 20;
                            }
                            p += w;
                        });
                        p
                    }).sum()
                };
                // a score is rounded to an eighth of a bit per residue
                let (more, less) = (exact(t - 0.2), exact(t + 0.2));
                assert!((less..=more).contains(&(row[i] as f64)), "{shape:?} T {t}: {} not in {less}..{more}", row[i]);
            }
            let quality = &row[THRESHOLDS..];
            assert!(quality.windows(2).all(|w| w[0] >= w[1] - 1e-3) && quality[0] > 0.0);
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
                lo.extend((0..21u8).map(|r| if r == best { 2.0 + (next() % 200) as f32 / 100.0 } else { -(next() as u32 as f32 / u32::MAX as f32) * 12.0 }));
            }
            offs.push(offs.last().unwrap() + n as u32);
            top.push(t);
        }
        let pf = Pfam { lo, offs };
        let (mut kept, mut dropped) = (0, 0);
        // at 40 a gap of more than 16 columns can still pay, which the wide band's last step is for
        for (min, band) in [17.0f32, 21.25, 26.25, 40.0].into_iter().flat_map(|min| [8i64, 16, 3, 12, 7, 15].map(|band| (min, band))) {
            // one model with a lower bar than the others, below what one bar for all may be
            let mins = [min - 0.005, min - 0.005, 13.0 - 0.005];
            let Some(b) = Bound::new(&pf, &mins, band) else { return }; // no SIMD here
            let (mut h, mut del) = (vec![], vec![]);
            for _ in 0..3000 {
                let mut pairs: Vec<(usize, Vec<u8>, i64)> = vec![];
                for _ in 0..PAIRS {
                    let m = (next() % 3) as usize;
                    let (mut f, mut c) = (vec![], (next() % lens[m] as u64) as usize);
                    let d = c as i64 + [0, -band, band, 1 - band, band - 1, band + 1][(next() % 6) as usize] + (next() % 3) as i64 - 1;
                    while f.len() < 1 + (next() % 60) as usize {
                        match next() % 12 {
                            0 => c += (next() % 4) as usize * if next() % 3 == 0 { band as usize / 2 } else { 1 },
                            // residues of the read in a gap
                            1 => {
                                let n = next() % 3 * if next() % 3 == 0 { band as u64 / 2 } else { 1 };
                                f.extend((0..=n).map(|_| (next() % 21) as u8));
                                continue;
                            }
                            2 => f.push((next() % 21) as u8),
                            _ => f.push(top[m].get(c).copied().unwrap_or(20)),
                        }
                        c += 1;
                    }
                    if next() % 4 == 0 {
                        // two clean stretches with one long gap in the model, from one edge of the band towards the other
                        let (c, n) = ((next() % 60) as usize, 8 + (next() % 10) as usize);
                        let gap = 1 + (next() % (2 * band as u64)) as usize;
                        let two = (c..c + n).chain(c + n + gap..c + 2 * n + gap);
                        pairs.push((2, two.map(|c| top[2].get(c).copied().unwrap_or(20)).collect(), c as i64 + band - (next() % 2) as i64));
                        continue;
                    }
                    pairs.push((m, f, d));
                }
                let reach = b.may_reach(std::array::from_fn(|i| (pairs[i].0, &pairs[i].1[..], pairs[i].2)));
                for ((m, f, d), reach) in pairs.iter().zip(reach) {
                    let g = gapped_band(&pf.lo, pf.offs[*m] as usize, lens[*m], f, *d, band, &mut h, &mut del);
                    assert!(reach || g < mins[*m], "dropped {g} at {}, band {band}", mins[*m]);
                    (kept, dropped) = (kept + (g >= mins[*m]) as usize, dropped + !reach as usize);
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
    fn survivors_round_trip_in_sort_order() {
        let mut all = vec![];
        for m in [0, 1, 65_534] {
            for o in [0, 7, MAX_BLOCK as u32 - 1] {
                for (c, j) in [(0, KEY_OFFSET), (0, 65_536), (3, 3), (0, 0), (65_535, 0)] {
                    let x = pack(m, o, c, j);
                    assert_eq!(unpack(x), (m, o, c as i64 - j as i64));
                    all.push((m, o, c as i64 - j as i64, x));
                }
            }
        }
        all.sort();
        assert!(all.windows(2).all(|w| w[0].3 <= w[1].3));
    }

    /// The same domain at the start of a sequence, and past position 65,535 and 131,071: the same scores at every place.
    #[test]
    fn join_finds_a_domain_wherever_it_is() {
        let dir = std::env::temp_dir().join(format!("kfpfam-long-{}", std::process::id()));
        fs::create_dir_all(&dir).unwrap();
        let path = |name: &str| dir.join(name).to_str().unwrap().to_string();
        let n = 60usize;
        let lo: Vec<f32> = (0..n * 21).map(|i| if i % 21 == (i / 21) * 7 % 20 { 4.0 } else { -2.0 }).collect();
        fs::write(path("profile.f32"), bytes_of(&lo)).unwrap();
        fs::write(path("offs.u32"), bytes_of(&[0u32, n as u32])).unwrap();
        let pf = Pfam::load(dir.to_str().unwrap());
        build(&pf, Shape::parse("4"), 7.0, &path("i.idx"));
        tile(&pf, &load_index(&path("i.idx")), &path("i.tidx"));
        let domain: String = (0..n).map(|c| A[c * 7 % 20] as char).collect();
        // filler of one residue, which scores -2 nearly everywhere and is in no word
        let at = [0usize, 100, 65_500, 65_536, 70_000, 140_000, 300_000];
        let faa: String = at.iter().map(|&j| format!(">at{j}\n{}{domain}{}\n", "W".repeat(j), "W".repeat(30))).collect();
        fs::write(path("p.faa"), faa).unwrap();
        join(&pf, &load_tiles(&path("i.tidx")), &path("p.faa"), 17.0, 16.5, 100, &path("out.tsv"));
        let out = fs::read_to_string(path("out.tsv")).unwrap();
        fs::remove_dir_all(&dir).unwrap();
        let rows: Vec<Vec<&str>> = out.lines().map(|l| l.split('\t').collect()).collect();
        assert_eq!(rows.iter().map(|r| r[0].to_string()).collect::<Vec<_>>(), at.map(|j| format!("at{j}")).to_vec());
        assert!(rows.iter().all(|r| r[1..] == rows[0][1..]), "{out}");
    }

    #[test]
    fn sketch_keeps_the_stated_fraction() {
        let kept = (0..1_000_000).filter(|&c| in_sketch(c, 8)).count();
        assert!((kept as f64 / 125_000.0 - 1.0).abs() < 0.02, "{kept}");
        assert!(in_sketch(12345, 1));
    }
}
