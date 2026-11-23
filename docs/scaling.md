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

## Ranked plan

| # | Change | Expected | Recall cost | Effort |
|---|---|---|---|---|
| 1 | `score4_sse41`, fused like `score4_neon` | most of 8× on x86 | **none** | a day; `score4_equals_scalar` already tests it |
| 2 | AVX2 (32 lanes) and an AVX-512 VBMI path | 2–4× beyond #1 | **none** | days; runtime feature detection already present |
| 3 | GA-only settings for RawR: bar 22, min_gapped ~24, not bar 17 | 1.5× | .999 → .995 ≥GA | config only |
| 4 | T7.5 for RawR instead of T6.5 | 2.3× | .995 → .984 ≥GA | config only |
| 5 | Per-model T as a length rule | 1.2× | none at equal hits | implemented (`FILTER_T_FILE`), needs the rule |
| 6 | Drop promiscuous k-mers | 1.5× | ~.001 | implemented (`FILTER_DROP`), needs a default |

**#1 + #3 is the 10×**, and only #3 touches recall. Measure #1 on Codon with `measure/speed.sh`
before and after; the x86 figure, not this laptop's, is the one that matters.

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
