# Scaling the filter: where a 10× is, and where it is not

Written after the RawR pipeline measured Mode A at 2,153 cpu-h against nail's 287 on eleven real
libraries (`nar-review-response/response/B11_v6_kfpfam_arms.md`), about 13× the per-read cost this
repo's README advertises. All measurements below are on this machine (M4 Pro, 6 threads, 309,149
ORFs / 9.26 M residues from ERR10889147, index `k4T6.5`) unless stated.

## Where the time goes

```
bar 22, gate 16.5: 32 core-s; join 28.5, gapped 1.6, survivors 2.7
18,882,393,580 hits, 1.5 ns per hit, 44,351,912 survivors (0.0023), 533,697 pairs
```

The join is **89% of the cost**, and it is `Σ over k-mer codes (index rows × query k-mers)`. Two
things follow, both measured:

**Thresholds cannot help.** The hit count is *identical* at bar 17 and bar 22 —
18,882,393,580 both times — because the bar filters after the window is scored. Raising it from 17
to 22 takes 47 core-s to 32 (the saving is all in the gapped stage and the output), and costs .999
→ .995 recall at GA. Worth doing for a GA-only consumer like RawR, which never uses the sub-GA
bins Mode A's bar of 17 was tuned for, but it is 1.5×, not 10×.

**Index density is the only tuning lever, and it is bounded.** At bar 22, GA-only:

| index | hits | core-s | vs T6.5 | recall ≥GA |
|---|---|---|---|---|
| T6.5 (Mode A) | 18.9 G | 32 | 1.0× | .995 |
| T7.5 (Mode B) | 6.39 G | 14 | 2.3× | .984 |
| T8 | 3.67 G | 10 | 3.2× | .965 |

The plan's own sweeps bound the rest: per-model T is worth about 1.2× as a length rule and
overfits if fitted per model; per-model seed shape about 1.5× at equal index size; dropping
promiscuous k-mers about 1.5× at equal recall, 2.7× at 5× the index. Multiplied out that is ~2.7×,
and the plan already concluded the Speed plan's 25× "is not in reach at any index size".

So **10× at equal recall is not available from tuning.** It is available somewhere else.

## The 10×: x86 has no fused kernel

`scorer4()` returns the fused `score4_neon` on aarch64 and falls through to `score4_generic`
everywhere else. There is no `score4_sse41`. On x86 the hot loop therefore:

1. calls the per-group `score_sse41` **four times** instead of once, reloading the 32-byte profile
   column for each job at every offset;
2. runs **one serial `cur = max(adds(cur, s), zero)` chain** per call, so it is latency-bound,
   where the fused kernel keeps four independent chains in flight;
3. computes the lane masks in **scalar** code (`nibbles`, 16 byte compares per job) where
   `score4_neon` does it in SIMD;
4. emulates NEON's single 32-entry `vqtbl2q_s8` with **five** SSE instructions — two `pshufb`, a
   `sub`, a `cmpgt` and a `blendv` — because `pshufb` only indexes 16 bytes.

`FILTER_SIMD=generic` now selects that path while keeping this CPU's per-group kernel, so the loss
is measurable anywhere:

| path | ns per hit | core-s | vs fused |
|---|---|---|---|
| `score4_neon` (fused, 4 jobs) | 1.5 | 32 | 1.0× |
| `FILTER_SIMD=generic` (per-group NEON ×4) | **33.9** | **643** | **22.5× slower** |
| `FILTER_SIMD=scalar` | 61.3 | 1161 | 36× slower |

**The fused kernel is worth 22.5×, and x86 does not have it.** Codon's measured gap is about 8×
per ORF residue, which is smaller — `score_sse41` is still SIMD, and the compiler recovers some of
the hoisting — but the shape of the deficit is this, and it is the single largest factor in why
the pipeline's cost bore no resemblance to the README's.

Nothing here costs recall. Writing `score4_sse41` reproduces `score4_neon`'s structure exactly:
one offset loop, the profile column loaded once, four independent accumulator chains, masks via
`_mm_movemask_epi8`. On AVX2, `_mm256_shuffle_epi8` over 32 lanes doubles throughput again for the
same instruction count per lane; on AVX-512 VBMI, `_mm512_permutexvar_epi8` is a true 64-entry
lookup in one instruction — NEON's `tbl2` but four times wider — which would put x86 ahead of this
laptop rather than behind it.

## One kernel design, three lookups (2026-11-23)

`score4_sse41` and `score4_avx2` are now written. They are **not** two different designs, and
should not become two: everything that makes `score4_neon` fast is architecture-independent and
they share it —

- all four jobs in one pass over the window, so each offset's loop overhead is paid once;
- each job's profile column loaded once per offset;
- independent accumulator chains, so the `adds`-then-`max` latency is hidden;
- lane masks computed in SIMD (`nibbles_sse` reproduces NEON's shift-and-narrow: compare, shift
  each u16 right by 4, `pshufb` the low bytes together, one `movq`).

The only part that differs is the one the hardware forces: **a 32-entry byte lookup for 16
lanes**. NEON does it in one `vqtbl2q_s8`. `pshufb` indexes only 16 bytes, so x86 needs two
shuffles, a compare and a blend. AVX2 doesn't remove that — `vpshufb` still indexes 16 bytes per
128-bit half — but it lets one register carry two jobs, each half with its own table, so the same
work issues half as many instructions. That is roughly 1.5×, not 2×; the limit is `pshufb`'s
16-byte window, not register width.

So what is good on x86 is good on ARM, and vice versa: the structure is the speed. The lookup is
three lines per architecture. Going genuinely wider — 32 lanes on AVX2, 64 on AVX-512 VBMI with
`vpermb`, a true 64-entry lookup in one instruction — means widening `LANES` in the postings
layout, which is shared by every kernel. That is the next step if x86 measurements say the
kernel, rather than memory, is still the limit.

### Verified

Homebrew's Rust here has only the aarch64 standard library, so the x86 code was checked in
`rust:1-slim-trixie` under `--platform linux/amd64`, which reports SSE4.1 and AVX2:

- **It compiles on x86_64** — which it would not have before: the first draft had a `*const i8`
  passed where `*const u8` was expected, invisible on ARM because the code was `cfg`'d out.
- **`score4_equals_scalar` now checks every kernel the CPU has**, not only the one `scorer4()`
  dispatches to, so SSE4.1 and AVX2 are both tested against scalar wherever the suite runs on x86.
  All four (generic, dispatched, sse4.1, avx2) pass.
- **End to end, byte-identical:** 20,000 real ORFs against `k4T8`, the x86 AVX2, SSE4.1 and
  generic kernels each write the same 26,319 pairs as ARM's NEON, `cmp` clean.

Under emulation the fused kernels ran at 13–15 ns per hit against generic's 39, about 3×. That is
direction, not magnitude: emulation translates x86 into ARM, and in particular says nothing about
AVX2 against SSE4.1 (they came out level). **The speed-up has to be measured on Codon.**

`FILTER_SIMD=sse4.1` forces the narrower kernel on an AVX2 machine. `FILTER_SIMD=generic`, the
path before these kernels, was removed once Codon had measured it (below); `f40154c` still has it.

## Codon, and what it left (2026-11-23)

`speed_test.sh` in `nar-review-response/hpc/B11_pfam_arms`: 800,000 ORFs (25.4 M residues) of
SRR23926920, Mode A (T6.5, bar 17, gate 16.5, gapped 21.25), 8 threads, 400,000 ORFs per block.
All three kernels wrote the same 4,498,149 pairs.

| kernel at `f40154c` | core-s | join | gapped | core-µs per residue |
|---|---|---|---|---|
| AVX2 | 347 | 165 | 177 | 13.7 |
| SSE4.1 | 374 | 186 | 183 | 14.7 |
| generic (what every B11 arm ran on) | 1,053 | 857 | 191 | 41 |

The fused kernels are worth 3.0x overall, 5.2x on the join. `generic`'s 41 agrees with the 46
estimated from the B11 trace, so Mode A's 2,153 cpu-h would be about 710: still 2.5x nail's 287.
Three things stood out:

- **The gapped stage was the largest**, 2.5 µs for each of 71.7 M pairs, of which 6% were written.
  Its loop was one scalar dependency chain per row.
- **AVX2 was 1.13x SSE4.1 on the join**, as half its instructions put two jobs into one register.
- **The kernel, not memory, limits the join.** llvm-mca gives 0.125 cycles per lane-offset on
  Cascade Lake; at about 43 lane-offsets per hit that is 1.8 ns at 3 GHz, and 1.5 was measured once
  the survivor stage is taken off.

What changed, with the same pairs out, byte for byte, as `f40154c` (NEON on 309,149 ORFs; AVX2,
SSE4.1 and NEON on 40,000):

- **`Bound`**: an upper bound on the banded gapped score, in bytes. The ±8 band is one 16-byte
  register, its 17th column merged into the 16th, one read residue per step; scores round up and gap
  costs down, so it should not drop a pair the f32 score keeps (it did one; see below), and the f32 score runs only on what it
  passes (16% of pairs here). Four pairs run side by side, as each one is a long dependency chain. 0.13 µs
  per pair on this laptop. It needs `FILTER_MIN_GAPPED`.
- **`gapped_band`** opens gaps from the cell before the gap, not after: opening from a gap never beats
  extending it, so the score is the same and the chain is half as long.
- **Survivors** are one `u64` each, so they sort as integers.
- **The gapped stage runs inside each tile**, so a block no longer holds its pairs.
- **x86 has 32 lanes, one job per register.** A profile column's halves are broadcast from memory,
  which costs a load and no shuffle; windows are stored plus 0x70, so `pshufb`'s zero for an index with
  its top bit set replaces the compare and blend; scores are held less 128, so the saturating add is
  also the floor. The last of these is on NEON too.

llvm-mca, cycles per lane-offset (predicted, not timed):

| | Cascade Lake | Ice Lake | Zen 2 | Zen 3 |
|---|---|---|---|---|
| AVX2 before | .125 | .116 | .146 | .126 |
| AVX2 now | .067 | .067 | .110 | .073 |
| SSE4.1 before | .188 | .168 | .169 | .113 |
| SSE4.1 now | .134 | .134 | .181 | .121 |

On this laptop, 3.2 M ORFs, 8 threads, 400,000 per block: **516 → 378 core-s**, gapped 152 → 59,
survivors 52 → 32, join 345 → 297 with survivors in it, peak 4.0 → 3.3 GB. The gapped figures are not
quite alike: before, the stage was charged for idle threads, and now it is not. Under Rosetta
the gapped stage goes 5.4 → 0.9 core-s on 40,000 ORFs. Rosetta runs the new AVX2 kernel at half
the SSE4.1 one's speed, so nothing about AVX2 can be read from it.

**Block size now matters.** A group is full only if its k-mer has enough postings in the block:
hit-weighted, with lane-offsets per hit (which is the kernel's work):

| ORFs per block | 16 lanes | 32 lanes | 64 lanes |
|---|---|---|---|
| 300,000 | .83, 45.7 | .69, 59.5 | .48, 90.6 |
| 800,000 | .93, 38.3 | .87, 43.3 | .75, 53.6 |
| 1,600,000 | .97, 35.8 | .93, 38.6 | .87, 43.5 |

At 300,000 the 32 lanes give back 30% of what the kernel gains; at 1.6 M, 8%. Memory is what it
costs: 4.9 GB at 1.6 M on this laptop. On Codon it made no difference (below).

### On Codon, and a pair lost (2026-11-23)

`ac87746`, same ORFs and settings, on an EPYC 9555 (Zen 5), the two builds one after the other:

| | core-s | join | gapped | survivors (in join) | peak |
|---|---|---|---|---|---|
| `f40154c`, AVX2 | 361 | 167 | 188 | 33 | 4.9 GB |
| `ac87746`, AVX2 | **157** | 123 | 28 | 23 | 3.7 GB |
| `ac87746`, SSE4.1 | 206 | 172 | 29 | 24 | 3.8 GB |

2.3x, and 6.6x on the gapped stage. But it wrote 4,498,148 pairs, **one short**: `Bound` was not an
upper bound. The band's 17th column is merged into the 16th lane, and a gap from the 17th to the 16th
moves on one residue without scoring it, while a cell that stays in a lane adds that residue's score,
which a conserved column can make worse than the gap. The last lane now keeps its own gapped value as
well as passing it on. The unit test had weak penalties and few edge diagonals and missed it; with both
it fails on `ac87746`, and the fix passes 48 M random pairs. It lets 0.02% more pairs through.

### The fix, AVX-512, and block size (2026-11-23)

The bound fixed, a window offset loaded once for a batch (nearly every batch is four rows against one
group), and `score4_avx512`: AVX-512 VBMI's `vpermb` looks up the whole 32-byte column in one
instruction. Same node, ORFs and settings; every run gave `f40154c`'s 4,498,149 pairs.

| | core-s | join | gapped | survivors (in join) |
|---|---|---|---|---|
| `f40154c`, AVX2 | 354 | 162 | 187 | 31 |
| AVX-512 (now the default where the CPU has it) | **125** | 92 | 28 | 21 |
| AVX2 | 137 | 104 | 28 | 21 |
| SSE4.1 | 194 | 162 | 28 | 22 |

Two predictions were wrong:

- **AVX-512 is 1.13x AVX2 on the join, not 3x.** llvm-mca assumes every load hits L1. With the
  arithmetic this cheap the join waits on something else. The profile's rows were the guess, as they
  outgrow L1 and a tile outgrows L2; the next section tests that, and it isn't so.
- **Block size doesn't matter.** 1.6 M ORFs at AVX-512: 262, 270 and 262 core-s in blocks of 400,000,
  800,000 and 1.6 M, for 3.9, 4.8 and 6.5 GB. Fuller groups save the kernel what longer survivor lists
  then cost (47 → 64 core-s). `kfpfam_block` stays at 400,000.

### Where it stops paying (2026-11-23)

`fd1f74c`: four rows meet every group of a k-mer in turn, where one group met every row, so a profile row
is read once per k-mer and tile; and a batch without a survivor skips the survivor scan. Same node, ORFs
and settings, AVX-512, the same pairs from every run:

| | core-s | join | gapped | survivors (in join) |
|---|---|---|---|---|
| Mode A, one group against every row | 121 | 89 | 27 | 19 |
| Mode A, four rows against every group | 121 | 90 | 27 | 20 |
| — tiles of 32,768 rows (315, against about 157) | 160 | 122 | 33 | 27 |
| — tiles of 16,384 rows (632) | 175 | 140 | 30 | 28 |
| Mode B, one group against every row | 52 | 45 | 2.3 | 13 |
| Mode B, four rows against every group | 56 | 49 | 2.4 | 15 |

The loop order makes no difference there (it stays: 7% on the join of the laptop). Mode B's 52 against 56
is one run each, and its survivor stage, which didn't change, is slower by as much in proportion.
Smaller tiles don't help: some of the slowdown is the node (the gapped and survivor stages, the same
work, are 12 to 40% slower in those runs) but not all of it, and L2 is therefore not what the join waits
on. `FILTER_TILE_ROWS` is gone again. What it does wait on is not known: three guesses here were wrong,
so the next step would be `perf stat` on a node rather than a fourth.

From the two modes, the join without survivors is about 0.6 ns per hit plus 14 to 18 core-s whatever the
hits: an eighth of Mode A's total, and a third of Mode B's. The per-hit part is within 1.6x of llvm-mca's figure
for the kernel alone.

In cpu-h, scaling Mode A's 2,153 by core-s on these ORFs: **Mode A about 250, Mode B about 110** (Mode B's
arm took 585). Both are against 53,275 models, and projections from one chunk; nail's 287 was against
30,134, so it needs measuring again on 39.0 before they are compared.

The lever left is fewer residues. The pipeline's ORF filter kept 62% of them for 0.18% of reads, and
the filter is linear in residues, so about 155 and 70 cpu-h if that holds for kfpfam. Not run with it.

### The bound at ±16 (2026-11-23)

`FILTER_BAND` up to 16 now has a bound too: 32 lanes, one AVX2 register on x86 and two on NEON. Before, anything but 8 went
straight to the f32 score. The 16-lane gap step takes its cost off in two halves, as 16 extensions don't fit in a byte. A
band between gets the next bound up, since a wider band only raises a score.

311,680 merged-fragment ORFs, Pfam 38.2, k = 3 with a T per model, M4 Pro, 4 threads, `FILTER_MIN_GAPPED` 21.25; the pairs written are the same byte for byte:

| | core-s | gapped | scored gapped |
|---|---|---|---|
| bar 15, gate 14.5, ±16 | 153 → **90** | 82.8 → 20.4 | 68.8 M → 3.9 M |
| bar 16, gate 15.5, ±16 | 103 → **74** | 41.7 → 12.0 | 33.9 M → 3.1 M |
| bar 15, gate 14.5, ±12 | 134 → 89 | 64.1 → 19.0 | 68.8 M → 3.9 M |
| Mode A, ±8 | 63 → 62 | 5.2 → 5.1 | unchanged |

The unit test covers bands 3, 8, 12 and 16, at a bar high enough (26 bits and up) for a gap of more than 16 columns to pay, and fails on each of twelve deliberate
faults in the step, the edges and the joins between registers (nine on ARM, three in the AVX2 code). It passes 25.6 M random pairs on ARM. Under Rosetta the AVX2 build passes it and
writes ARM's pairs on 4,000 ORFs. Not timed on x86.

### Spaced seeds (2026-11-23)

`filter build PREP SHAPE T` takes a shape, 1011001 say, in place of k: the residues a word takes. The shape is in the index's header
(a k-mer's is k, as before), and only the index and the cutting of reads into words change: not the kernel. Indexes for a k-mer
are byte for byte what they were. `search`, `window` and `sketch` read k-mers only.

`filter seeds PREP GA.tsv OUT SHAPES` tabulates, for each model, shape and T, background hits and -ln P(a read that just reaches GA
holds no word), to choose a T per model. `FILTER_DIAGONALS=all` scores a pair's band around every diagonal at the bar.

For recall they are better than contiguous words at equal cost, and catch repeats with a period (Collagen's
Gly-X-Y) that a 3-mer cannot: raw-reads-analysis-pipeline, `docs/plans/kfpfam-recall.md`.

## Ranked plan

| # | Change | Expected | Recall cost | Effort |
|---|---|---|---|---|
| 1 | `score4_sse41`, fused like `score4_neon` | most of 8× on x86 | **none** | **done**: 5.2x on the join on Codon |
| 2a | AVX2, two jobs per register | 1.13× beyond #1 on Codon | **none** | **done**, replaced by 2b |
| 2b | 32 lanes, one job per AVX2 register | 1.3–1.9× beyond #2a (predicted) | **none** | **done**: join 167 → 123 core-s on Codon at 400,000 per block |
| 2c | SIMD bound before the gapped score | 6.6× on that stage on Codon | **none** | **done** |
| 2c′ | The same at a band of ±16 | 4× on that stage on ARM | **none** | **done**, not timed on x86 |
| 2d | AVX-512 VBMI | 1.13× on the join on Codon | **none** | **done** |
| 2e | Four rows against every group; smaller tiles | none on Codon; smaller tiles slower | **none** | loop kept, tiles dropped |
| 2f | The pipeline's ORF filter in front | 1.6× | 0.18% of reads with `hmmsearch` | config only; not run with kfpfam |
| 3 | GA-only settings for RawR: bar 22, min_gapped ~24, not bar 17 | 1.5× | .999 → .995 ≥GA | config only |
| 4 | T7.5 for RawR instead of T6.5 | 2.3× | .995 → .984 ≥GA | config only |
| 5 | Per-model T as a length rule | 1.2× | none at equal hits | implemented (`FILTER_T_FILE`), needs the rule |
| 6 | Drop promiscuous k-mers | 1.5× | ~.001 | implemented (`FILTER_DROP`), needs a default |

#1 to #2c cost no recall, so #3 and #4 wait on what they measure on Codon. The x86 figure, not this
laptop's, is the one that matters.

Two things *not* worth pursuing on the evidence: more threads (Codon already ran at 95.7%
efficiency on 16 cores, against nail's 50.7%), and narrowing the 48-residue window — the `span`
tracking in `join` already restricts the scored range per batch, and consecutive jobs in the inner
loop share a group, so the span is already tight.

## The 183 families: not promiscuity, and not this repo

RawR found 183 families that nail finds and both kfpfam modes returned **zero** of. The hypothesis
that their k-mers were all promiscuous and dropped is **ruled out**, twice over:

- `build` has no promiscuity cap at all. Its `keep` closure is `|_, _| true`, so a k-mer is never
  excluded for being common. `FILTER_DROP` drops promiscuous k-mers at *search* time, and only
  when set — the pipeline does not set it.
- All 183 have ample index entries at T 6.5: **minimum 1,352, median 13,280, maximum 151,931**
  (the same enumeration `build` does, over `profile.f32`). The sparsest is PF07981 Plasmod_MYXSPDY,
  17 columns and still 1,352 entries.

They are also unremarkable as models — median length 142 aa against 120 for the rest, median GA
27.0 against 27.0, mean information 0.69 against 0.73 — and locally **Mode A emits candidates for
all 183**, and for 30,131 of 30,134 models overall. So the loss is in the prep or index built on
Codon by `hpc/pfam-arms/kfpfam_setup.sh`, not in the method: confirmed 2026-11-23, they came from another
Pfam release, with 53,275 models. See
[benchmark representativeness](benchmark-representativeness.md) for the diagnostic.

The only genuine blind spot here is three models too short to carry a 4-mer above any bar:
PF02044 Bombesin (14 columns), PF08261 Carcinustatin (8), PF03373 Octapeptide (8).
`representativeness.py models` reports them, and a rise above three is a build failure.

**Unknown residues.** X scored -100 and ended every alignment; HMMER scores it at the column's mean over the background and carries
on. The filter now does too, when it loads the profile, so an old prep needs no rebuild (seeds still hold no X). In bird gut's B11 reads
12 of the 21 pairs the best mode lost held an X. 7,796 V3 pairs with one residue masked: 7 lost before, 1 after.

**`FILTER_GAPPED_FILE`**: a gapped bar per model, one per line, where that is lower than `FILTER_MIN_GAPPED`. A short repeat reaches GA1 as
several copies, each under any bar that suits the rest. The bars are the raw score of a domain at GA2 less 4 bits; on V3's ORFs they add 6% to candidates at
bar 15 and 10-12% at 17 and 20, and nothing to the filter's time.

**A band of 7 or 15 has a tighter bound.** Bands of 8 and 16 are a column wider than the register, so the bound merges the last two columns, and
on whole proteins that lane climbs by chance: the bound passed 45% of pairs at bar 20 (band 8) where 6% pass at 7, and 48% against 18% at 16 and 15 (3,013
Pfam models, 50,000 reversed MGnify proteins of 165 residues). On V3's ORFs it matters less: bar 15 at 15 and `any` takes 90 CPU-s against 95, with the same
8,533 of 8,533 pairs; bar 17 at 7 takes 63 against 63. The bound for 8 and 16 is as it was.

**`FILTER_DIAGONALS=any`** writes the pairs `all` does, and leaves a pair at the first diagonal that reaches its bar: 188,485 scored in full against 336,092 on
proteins that hold many true pairs. Its scores are that diagonal's, not the pair's best.

## Whole proteins: where the time goes, and a window for the gapped stage (2026-11-26)

24,000 MGnify proteins (6.2 M residues, median 191) and the same residues cut into pieces, 7,534 Pfam 38.2 models
(every fourth), M4 Pro, 4 threads, band 15, `any`. The new stage timers (`stages, core-s:` in the log) split the join.
Candidates in, CPU-s, and cost per candidate; `fast` (θ 5, bar 18, gate 17.5, gapped 10 under each model's weakest
domain):

| stage | 50-residue pieces | whole proteins | per candidate |
| --- | --- | --- | --- |
| window gate (k-mer hits) | 2.40 G, 3.3 s | 2.75 G, 4.2 s | 1.37 → 1.53 ns |
| whole diagonal (unique diagonals) | 3.20 M, 0.5 s | 4.74 M, 1.8 s | 163 → 370 ns |
| gapped bound (pairs at the bar) | 2.13 M, 0.5 s | 3.27 M, 3.6 s | 0.24 → 1.11 µs |
| gapped exact (past the bound) | 95 k, 0.1 s | 274 k, 1.3 s | 1.4 → 4.6 µs |
| HMMER (pairs written) | 43.9 k, 2.4 s | 53.6 k, 9.2 s | 56 → 171 µs |

`top` (θ 9, bar 15, gate 14.5, gapped 21.25): gate 9.7 → 11.2 G hits, 12.2 → 21.8 s; diagonals 31 → 44 M, 4.8 → 19.1 s
(154 → 431 ns); bound 20 → 30 M pairs, 4.9 → 33.5 s (0.24 → 1.13 µs); exact 1.0 → 2.7 M, 1.5 → 15.6 s; HMMER 0.45 → 1.02 M
pairs, 7.1 → 41.4 s. So it is mostly the cost per candidate that grows, roughly with the stretch scored (the whole
diagonal, then the band over the model's length), and only partly the count (1.2-1.5×, the exact stage 2.6-2.9×).
The gate's per-hit cost hardly moves.

`FILTER_GAPPED_WINDOW=w` (off by default; the same pairs as before without it) gives the bound and the exact score only
the residues from w before the diagonal's best ungapped segment to w after it. Whole proteins, pairs missed of 8,723
that `hmmsearch` passes (HMMER with `--forward-slack 8`):

| w | `fast` filter CPU-s | pairs to HMMER | missed | `top` filter CPU-s | pairs to HMMER | missed |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| whole diagonal | 10.6 | 53,567 | 69 | 85.1 | 1,019,878 | 0 |
| 64 | 8.8 | 51,678 | 77 | 65.2 | 1,006,323 | 0 |
| 32 | 7.8 | 47,710 | 122 | 55.8 | 949,650 | 4 |
| 16 | 7.2 | 42,827 | 225 | 49.7 | 847,505 | 24 |
| 8 | 6.9 | 39,140 | 356 | 46.2 | 765,576 | 40 |
| 0 | 6.6 | 36,214 | 577 | 43.0 | 715,172 | 47 |

At 64 `top` loses nothing and its filter takes 0.77× the time; `fast` pays 8 pairs for 0.83×. Windows narrower than
that lose recall faster than they save, and barely change HMMER's work (its pairs fall less than its cost per pair
rises). The gate's figures vary by about 3 s between runs of the same work, so read differences under that as noise.
What is left of `top`'s gap to the pieces is the whole-diagonal stage (19 s at 431 ns a diagonal), which a window
would need the k-mer's position for: survivors keep only (model, protein, diagonal).

### A window for the whole-diagonal stage too (2026-11-26)

Survivors now keep the read position of each diagonal's first and last k-mer through the gate (a u128 per survivor, compacted
to at most two a diagonal; the pairs written are byte for byte the same, and sorting `top`'s survivors takes 1.4 CPU-s against
0.6). `FILTER_DIAG_WINDOW=w` scores the diagonal on the read from w before the first to w after the last, not on the whole
stretch that meets the model; its best segment then places `FILTER_GAPPED_WINDOW`. Same proteins and models as above; HMMER with
`--forward-slack 8`; pairs missed of 8,723:

| mode | diagonal window | gapped window | diagonal stage | pairs at the bar | gapped bound + exact | filter | HMMER | missed |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `top` | whole | whole | 18.9 | 29.6 M | 49.4 | 86.5 | 41.6 | 0 |
| `top` | whole | 64 | 18.9 | 29.6 M | 29.4 | 66.5 | 41.5 | 0 |
| `top` | 32 | 64 | 8.1 | 29.5 M | 29.0 | 55.4 | 41.4 | 0 |
| `top` | 16 | 64 | 5.7 | 26.8 M | 25.7 | 49.5 | 40.2 | 0 |
| **`top`** | **8** | **64** | 4.3 | 20.4 M | 19.1 | **41.5** | **37.1** | **0** |
| `top` | 8 | 32 | 4.4 | 20.4 M | 12.2 | 34.7 | 36.2 | 4 |
| `top` | 4 | 64 | 3.6 | 16.0 M | 15.0 | 36.6 | 34.6 | 4 |
| `top` | 0 | 64 | 2.4 | 9.8 M | 9.7 | 30.0 | 31.1 | 13 |
| `fast` | whole | whole | 1.8 | 3.27 M | 5.0 | 10.9 | 9.4 | 69 |
| `fast` | 32 | whole | 0.8 | 3.26 M | 5.0 | 9.9 | 9.4 | 69 |
| `fast` | 32 | 64 | 0.8 | 3.26 M | 2.9 | 7.8 | 9.4 | 77 |
| `fast` | 16 | 64 | 0.6 | 2.82 M | 2.4 | 7.1 | 9.2 | 102 |

The diagonal window does two things: each diagonal costs less (428 → 98 ns at 8), and a score over a stretch near the k-mers is a
sharper test than one over the whole diagonal, which can collect chance segments far away, so fewer pairs reach the bar
(29.6 → 20.4 M) and every later stage, HMMER included, has less to do. `top` at 8 and 64 is 0.61× its CPU with filter and HMMER
together, 128 → 79 CPU-s, with no pair lost; on these proteins that is 2.3× `hmmsearch` (179) where it was 1.4×. `fast`'s
diagonals were cheap already: 32 and whole is free (69 lost, as before), 32 and 64 saves a little more for 8 pairs.
What is left of `top`'s filter is the window gate itself, 17-18 CPU-s of 41.5.

### Long queries or long models? (2026-11-26)

The same 7,534 models cut into 24,367 pieces of 50 columns (25 to 74; T and a flat gapped bar of 21.25 from the parent; a
timing test only, as pieces have no thresholds), `top`'s settings, whole diagonals and bands, filter CPU-s:

| | models whole (median 121 columns) | models cut to 50 |
| --- | ---: | ---: |
| 50-residue queries | 25.8 (diagonal 161 ns, bound 0.24 µs) | 20.7 (140 ns, 0.25 µs) |
| whole proteins (median 191) | 89.1 (427 ns, 1.13 µs) | 31.8 (145 ns, 0.40 µs) |

The diagonal and gapped stages cost what the overlap of query and model costs, the shorter of the two, so shortening either one
brings them back: whole proteins against 50-column models cost 1.2× what 50-residue pieces do, not 3.5×. The windows of the
sections above do the same without cutting anything (41.5). HMMER is different: a candidate's Forward, Backward and domain
definition run over the whole query whatever the model, and `hmmsearch` pays for every column too, so model length cancels
out of HMMER's speed-up and query length does not (`fast`: 56 µs a candidate on 50-residue pieces, 171 on whole proteins,
for about as many candidates).

### One windowed score or two? (2026-11-26)

With `FILTER_DIAG_WINDOW=8` a diagonal of one k-mer is scored on 24 residues, inside the gate's 48, so the gate's outer residues
only decide which hits reach a test that no longer looks at them. Builds with `WIN` 32 (`HALF` 12) and 24 (`HALF` 8), their own
tiles, the same proteins; total is filter and HMMER (`--forward-slack 8`), CPU-s, pairs missed of 8,723:

| mode | gate window | gate | diagonal, gapped windows | gate CPU-s | diagonals | at the bar | filter | total | missed |
| --- | ---: | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `top` | 48 | 14.5 | 8, 64 | 17.3 | 44.2 M | 20.4 M | 41.4 | 78.4 | 0 |
| `top` | 32 | 14.5 | 8, 64 | 14.2 | 36.1 M | 18.7 M | 33.7 | 68.9 | 3 |
| `top` | 24 | 14.5 | 8, 64 | 12.8 | 28.5 M | 17.8 M | 29.4 | 63.6 | 5 |
| `top` | 24 | 14 | 8, 64 | 15.1 | 40.9 M | 17.9 M | 32.7 | 67.0 | 5 |
| `top` | 24 | 14.5 | 0, 64 | 11.4 | 28.5 M | 7.3 M | 18.5 | 43.6 | 64 |
| `fast` | 48 | 17.5 | 32, 64 | 4.0 | 4.7 M | 3.26 M | 7.9 | 17.3 | 77 |
| `fast` | 32 | 16.5 | 32, 64 | 3.5 | 6.9 M | 2.72 M | 6.8 | 16.0 | 101 |
| `fast` | 24 | 17.5 | 32, 64 | 3.0 | 2.5 M | 1.79 M | 4.5 | 13.2 | 207 |

- The gate costs 0.74× at half the window, not 0.5×: its 11.2 G hits each cost something whatever the window.
- A narrower gate is also a stricter one (fewer survivors), and the few pairs it loses do not come back with a lower gate,
  which only lets more through (14 against 14.5 costs 3.3 CPU-s and keeps the same 5 lost).
- Merging the two into one stage would cost more, not less. The gate is int8, per hit, in SIMD; the diagonal is exact f32, once
  per diagonal (0.25-0.4% of hits) and over every k-mer on it; it halves what reaches the bound (44 → 20 M diagonals at the bar).
  Without it the bound would see every survivor, and int8 rounding alone moved 5% of pairs across the bar before.
- So two stages stay, with a narrower gate as a choice: `top` 32 is 0.88× its CPU for 3 pairs (0.03%), 24 is 0.81× for 5.
  `WIN` is a constant of the postings' layout and the tiles: as an option it would need the tiles built for it.

`FILTER_WINDOW=len:half` at `filter tile` now sets the gate's window (default 48:22, as before); the tiles carry it and `join`
reads it from them, so one binary serves both. A 32:12 index gives the pairs of a build with `WIN` 32, byte for byte, and a
default one those of before. `top` in hmm-annotate is 32:12 with `FILTER_DIAG_WINDOW=8` and `FILTER_GAPPED_WINDOW=64`
(on 50-residue pieces too, 0 of 11,835 pairs lost, 22.1 → 18.5 CPU-s for the filter).

### Models trimmed to 50 columns, and where `top`'s time goes now (2026-11-26)

The same 24,000 proteins and 7,534 models; `top` at 32:12, 8, 64; `hmmsearch` 176 CPU-s; HMMER with `--forward-slack 8`;
8,723 pairs to find. Each stage, what reaches it, and its CPU-s:

| stage | in | CPU-s |
| --- | ---: | ---: |
| window gate | 11.2 G k-mer hits | 13.6 |
| survivors sorted, whole diagonals | 85.0 M survivors, 36.1 M diagonals | 1.6 + 3.5 |
| gapped bound | 18.7 M pairs at the diagonal bar | 11.5 |
| gapped exact | 1.64 M past the bound | 6.0 |
| HMMER MSV | 752 k pairs written | 4.2 |
| HMMER bias filter and Viterbi | MSV's passes | 12.9 |
| HMMER Forward | Viterbi's | 9.9 |
| HMMER domain definition | 12.7 k within 8 bits | 7.6 |

69 CPU-s, 2.5× `hmmsearch`, 3 missed. The false pairs sit on the chance tail just over the bar (median gapped score
22.8 at a bar of 21.25), about half as many for each 1.75 bits more.

**One 50-column stretch per model** (the 50 with the highest summed best log-odds; models of 75 or fewer whole; parent T and
bars). A stand-in for a real choice of trim. 1.21 M columns become 0.38 M, index entries 348 M → 136 M; filter 10.3 CPU-s
(gate 6.5), HMMER 13.3 on 217 k pairs, 23.6 in all (7.5×), but 1,360 pairs missed. By how much of the stretch the true domain
covers: none, 629 of 630 missed; under half, 291 of 528; half to 80%, 92 of 759; 80-99%, 81 of 1,448; all of it, 267 of 5,358.
The last are weak domains whose 50 columns score 10-21 in the filter, under its bar. Choosing each model's stretch with
hindsight, to cover the most of these very pairs, still leaves 331 pairs that overlap it by fewer than 10 columns, 600 by fewer
than 25: MGnify's hits are often part of a domain, at different places. One stretch per model cannot keep `top`'s recall;
pieces covering the model can, and cost what the windows above cost.

**A second gapped bar on the whole diagonal** (`FILTER_SECOND_GAPPED=s`, uncommitted; not for the repeat models of
`FILTER_GAPPED_FILE`). The windows lower the weakest true pairs' scores (best of their diagonals: 22.7 at the least; on the
whole diagonal's band, 26.4), so a pair that passes the window is scored again over its whole diagonal and needs s there:

| | filter | pairs to HMMER | HMMER | total | × `hmmsearch` | missed |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `top` | 34.6 | 752 k | 34.6 | 69.2 | 2.5 | 3 |
| s 24 | 42 | 253 k | 20.3 | 62 | 2.8 | 4 |
| s 25 | 43 | 183 k | 16.9 | 60 | 2.9 | 4 |
| s 26 | 43 | 139 k | 14.3 | 57 | 3.1 | 4 |
| s 25, θ 8 | 38 | 177 k | 16.5 | 55 | 3.2 | 6 |
| s 25, θ 7 | 34 | 170 k | 16.3 | 50 | 3.5 | 7 |
| s 25, diagonal bar 16 | 34 | 163 k | 16.3 | 50 | 3.5 | 9 |
| s 25, diagonal bar 17 | 29 | 142 k | 15.6 | 45 | 3.9 | 13 |

HMMER halves for a quarter of the pairs: the pairs left are the longer proteins (false ones 534 residues on average).

**Cheaper HMMER, tried and dropped.** A Viterbi filter set like `--forward-slack` (a P-value that stops a protein 12 bits
under the threshold) saves 12 CPU-s and misses 182 more pairs: Viterbi is far under the final score for weak, gapped
domains. HMMER first on a window of the protein (the diagonals' stretch of the model, 20 or 40 residues more each side; then
the whole protein for those within 8 bits) reads a third of the residues, but costs 24.8 against 34.6 and misses 455
(repeats: EGF, cadherin, whose copies add up outside the window). On the second bar's pairs it saves 2.4 CPU-s.

**Floors.** The gate alone is 13.6 CPU-s, 13× `hmmsearch` if all else were free; HMMER on the true pairs alone, 3 CPU-s. 10×
needs the gate's 11 G hits made several times cheaper or fewer, which lower θ does only by missing pairs.

### How concentrated a model's evidence is (2026-11-26)

Per domain, the filter's scores along HMMER's alignment (match log-odds; gaps 6, and 1 a residue more) give the best local
score inside any w columns, a lower bound on what a model cut to those columns scores. Pfam 38.2's 412,351 seed domains
at GA, and the slice's 11,622 true domains (8,723 pairs); bar 21.25.

- **Each domain is concentrated.** For all but 368, 107, 34 and 2 of 8,723 pairs some 20, 30, 50 or 75 columns alone
  reach the bar; their best 50 columns hold a median 75% of the domain's score (10th percentile 46%).
- **But not in the same place.** Regions chosen on the seeds alone (each disjoint from the last, up to four) miss many:
  one of 50 columns keeps 6,358 pairs, four keep 8,271. Seed domains are strong and whole; MGnify's are weak and partial.
- **Regions chosen on the slice's other half** (with the seeds; greedy set cover, overlapping allowed, a training domain
  must clear the bar by a margin where it can), tested on the held-out half:

| w | margin | regions a model | columns | test pairs lost | of which no w columns can keep |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 30 | 0 | 1.13 | 0.21 | 108 of 4,055 | 53 |
| 30 | 4 | 1.22 | 0.23 | 79 | 53 |
| 50 | 0 | 1.07 | 0.32 | 61 of 3,515 | 14 |
| 50 | 8 | 1.13 | 0.34 | 40 | 14 |
| 50 | 12 | 1.17 | 0.35 | 37 | 14 |

(test pairs of models longer than w; shorter models stay whole.)

The filter with those regions as models (parent T and bars), HMMER on the whole proteins against the whole models:

| | pieces | columns | filter | HMMER | total | × `hmmsearch` | missed, all | missed, held-out half |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `top` | 7,534 | 1.21 M | 34.6 | 34.6 | 69 | 2.5 | 3 | |
| 50 columns, margin 8 | 8,400 | 0.41 M | 11 | 14.1 | 25 | 7.0 | 64 | 34 of 4,307 |
| 75 columns, margin 8 | 7,957 | 0.56 M | 16 | 18.7 | 35 | 5.1 | 25 | 16 of 4,307 |

Of the held-out misses, about 40% have no region of that width at all (weak domains spread over the model) and the rest
want a region the training half never showed. More training proteins would shrink the second part, not the first.

### Back to the reads' filter (2026-11-26)

kfpfam's use is short reads; MGnify proteins are annotated with `hmmsearch`, and `top` on them stays as a benchmark of the
method (its misses show where recall goes). On 309,149 read ORFs of ERR10889147 (9.3 M residues, the same 7,534 models,
`hmmsearch` 357 CPU-s, 2,142 pairs), `top`'s read settings, filter CPU-s (two runs each), all pairs found in every row:

| filter | CPU-s |
| --- | ---: |
| 210c93c, before the whole-protein work | 22.4 |
| 7a6cce4: u32 positions, survivors compacted again only when doubled, writes checked | 22.3 |
| c69307e: stage timers, `FILTER_GAPPED_WINDOW` (off) | 26.2 |
| cdaa97e: u128 survivors, `FILTER_DIAG_WINDOW` (off) | 27.1 |
| 63a4cf6: the gate's window from the tiles | 27.4 |
| 63a4cf6 with 32:12, 8, 64 | 22.1 |
| 63a4cf6 with 48:22, 8, 64 | 25.7 |

The windows only won back what their own code cost on reads, so the filter is 7a6cce4's again: the window options,
`FILTER_WINDOW` and the stage timers are gone (63a4cf6 has them), and an index tiled with `FILTER_WINDOW` other than the
default does not load. The sections above stand as the record of what whole proteins need. `--forward-slack 8` in HMMER
(hm-annotate) does help reads: 4.14 → 3.37 CPU-s, nothing lost.

Short models are the other side of this: the filter's cost follows the shorter of query and model, and a model of 50 columns
or fewer has no region to choose, so short families such as weefams suit it as they are.
