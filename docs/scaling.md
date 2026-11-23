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
  costs down, so it can't drop a pair the f32 score keeps, and the f32 score runs only on what it
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
costs: 4.9 GB at 1.6 M on this laptop. `kfpfam_block` in the pipeline is 400,000 until
`speed_test.sh` has timed 800,000 and 1.6 M on Codon.

**Next, if the join is still the limit:** AVX-512 VBMI's `vpermb` is a 32-entry lookup in one
instruction, about half the µops per lane again. `speed_test.sh` prints whether the node has it.
Nothing is written for it, as it can't be run here.

## Ranked plan

| # | Change | Expected | Recall cost | Effort |
|---|---|---|---|---|
| 1 | `score4_sse41`, fused like `score4_neon` | most of 8× on x86 | **none** | **done**: 5.2x on the join on Codon |
| 2a | AVX2, two jobs per register | 1.13× beyond #1 on Codon | **none** | **done**, replaced by 2b |
| 2b | 32 lanes, one job per AVX2 register | 1.3–1.9× beyond #2a (predicted) | **none** | **done**, verified; time on Codon, with larger blocks |
| 2c | SIMD bound before the gapped score | 6× on that stage under Rosetta | **none** | **done**, verified; time on Codon |
| 2d | AVX-512 VBMI | ~2× on the join | **none** | if Codon's nodes have it |
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
Codon by `hpc/pfam-arms/kfpfam_setup.sh`, not in the method. See
[benchmark representativeness](benchmark-representativeness.md) for the diagnostic.

The only genuine blind spot here is three models too short to carry a 4-mer above any bar:
PF02044 Bombesin (14 columns), PF08261 Carcinustatin (8), PF03373 Octapeptide (8).
`representativeness.py models` reports them, and a rise above three is a build failure.
