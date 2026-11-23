# Why the local benchmark said .999 recall and 1.20 core-ms, and the pipeline said otherwise

The RawR pipeline ran kfpfam Modes A and B against nail on eleven real metagenome libraries at full
depth (102.9 Gbase, ~500 M merged fragments; `nar-review-response/response/B11_v6_kfpfam_arms.md`).
Two headline numbers from this repo's README did not survive:

| | local (this repo) | RawR, eleven libraries | gap |
|---|---|---|---|
| Mode A recall ≥GA | .999 | 183 families found **zero** times | — |
| Mode A cost | 1.20 core-ms per read | ~16 core-ms per read | **13×** |
| Mode B cost | 0.37 core-ms per read | ~4.3 core-ms per read | **12×** |
| nail cost | 1.24 core-ms per read | ~2.1 core-ms per read | 1.7× |

Neither gap is a measurement error in the pipeline. Both are blind spots in how this repo measures.
Reproduce the coverage figures with `measure/representativeness.py`.

## 1. The recall figure cannot see a missing model

`measure/index_score.py` computes recall as the share of (read, model) truth pairs at or above GA
that the filter emits as candidates, where truth is `hmmsearch --cut_ga`-equivalent over all models
on the same reads. That is the right quantity. The problem is the denominator: it comes from 50,000
reads of one human gut run (ERR10889147).

```
local truth: 8,431 (read, model) pairs at or above GA, over 1,975 families
  that is 6.6% of the 30,134 models in the index; 28,159 models are never tested at GA
  one pair is 0.012% of recall, so a defect costing fewer than about 42 pairs reads as 0.995 or better
```

**The benchmark exercises 6.6% of Pfam.** Weighted by what real data actually annotates, it is
better but still partial: the families it covers hold 63.1% of the pipeline's 135.6 M annotated
reads, so **36.9% of real annotated reads come from 18,842 families the benchmark never tests**.
The largest are ordinary, abundant families — PF03466, PF00577, PF01041, PF06458, PF07715 — not
exotica.

That is exactly the size of hole the pipeline fell into. Of the 183 families the pipeline found
zero of, only **15** appear above GA in the local truth, contributing **39 of 8,431 pairs**. Had the
local index missed all 183, recall would have read **0.9954** instead of 1.000 — indistinguishable
from the .999 in the README.

A single pooled ratio also hides all-or-nothing failures. 183 families at zero and 183 families at
95% give the same recall if the pair counts match, but only the first is a bug.

### Those 183 are not a kfpfam defect

Worth stating plainly, because the pipeline report could be read that way. Locally, Mode A emits
candidates for **all 183**, and for 30,131 of 30,134 models overall. They are all present in
`data/prep/models.tsv`, and they are unremarkable — median length 142 aa, median GA 27.0, mean
information 0.69, against 120 aa / 27.0 / 0.73 for the rest. So the loss is in the prep or index
**built on Codon by `hpc/pfam-arms/kfpfam_setup.sh`**, not in the method. Pin it down there:

```bash
grep -c '^NAME ' pfam_unchunked/Pfam-A.hmm && wc -l < kfpfam/prep/models.tsv
```

Those should be N and N+1. If they are, the HMM and the prep agree and the loss is in
`filter build`/`tile`; if the prep is short, `prep.py` or the joined HMM is.

The only genuine local blind spot the new invariant finds is three models too short for a 4-mer
index at any bar — PF02044 Bombesin (14 aa), PF08261 Carcinustatin (8 aa), PF03373 Octapeptide
(8 aa). Worth a note in the README; not worth fixing.

## 2. The cost figure is a laptop figure

Three things were ruled out before blaming hardware.

**Thread scaling is fine.** On Codon `KFPFAM_FILTER` used 15.3 of its 16 cores (95.7% efficiency);
nail used 8.1 of 16 (50.7%). kfpfam's cpu-hours are not inflated by poor parallelism — it genuinely
does the work. In wall-clock terms Mode A took 141 h against nail's 36 h.

**Cost is linear in ORF residues, not superlinear in ORF length.** At 3 M residues held constant,
Mode A's wall time is flat across ORF length: 4.28 s at 30 aa, 4.06 s at 45 aa, 4.98 s at 68 aa.
What does rise is the candidate count — 320k, 585k, 844k pairs — which is downstream HMMER's cost,
not the filter's.

**So the difference is per-core throughput.** Two points on this machine (6 threads, M4 Pro, NEON):
3 M residues in 4.63 s, 12 M in 12.68 s, giving a marginal **5.4 core-µs per ORF residue** and
about 2 s of fixed index loading. Codon's Mode A spent 2,153 cpu-h on an estimated 1.5–1.9 × 10¹¹
ORF residues, i.e. roughly **46 core-µs per residue — about 8× this machine's**, ±25% on the
residue estimate. nail's equivalent gap is 1.7×, because MMseqs2 is heavily tuned on x86.

The candidates are clock speed and, more likely given the size of the gap, that the **SSE4.1 path
is much less optimised than the NEON path** (the Dockerfile already notes x86 takes SSE4.1). An
order-of-magnitude gap is not explained by clock alone.

There is a second, smaller distortion: the local reads are ~124 nt (ORFs max out at 41 aa, mean
29.9), while the pipeline searches ~206 nt merged fragments whose ORFs reach ~68 aa. Per *read*
comparisons therefore mix hardware with read length; per *residue* is the fair unit, and the
per-read column above overstates the gap by roughly the length ratio.

## What to change

1. **Measure cost on the target hardware.** `measure/speed.sh` runs anywhere. Run it on an x86 node
   and quote both numbers, or quote core-µs per ORF residue rather than core-ms per read. Until
   then no cost claim in the README transfers off this laptop. If the x86 figure is ~8× the Apple
   one, the SSE4.1 path is the thing to work on, and that is a bigger win than any retuning.
2. **Add the model-coverage invariant to every engine run.** It is one line and it would have caught
   this class of bug:
   ```bash
   python3 measure/representativeness.py models data/prep data/speed/NAME/real.tsv
   ```
   Every model that can be seeded at all must appear at least once across a large ORF set. Treat a
   rise above the known three as a failure, not a metric.
3. **Report recall per family as well as pooled.** Pooled recall over 1,975 families cannot show an
   all-or-nothing loss in one of them. The minimum is: number of families with a truth pair at GA
   whose recall is 0, and the read-weighted share they carry.
4. **Widen the truth set beyond one gut run.** Two cheap changes, in order of value:
   - Weight by what real data annotates. The families holding the other 36.9% of real reads are
     known; `measure/sets.py fragments` already builds family-balanced fragments from proteins, so
     the machinery exists to cover them deliberately rather than by whatever the sample contained.
   - Add a second and third biome. The pipeline's own set (six Zymo mocks plus five biome runs) is
     the obvious source, and the mocks come with a truth of their own.
5. **Match the read geometry.** `measure/rawr_orfs.sh` already uses the pipeline's exact translation
   settings, which is right, but the benchmark's reads are ~124 nt. Build the standard ORF set from
   *merged* fragments of the length the pipeline actually produces (~200 nt), so the ORF length
   distribution and the candidate counts per read are the real ones.
6. **Quote recall at the read level too.** The pipeline's unit is a read that ends up annotated,
   after `PARSEHMMSEARCHCOVERAGE` keeps essentially one family per read (1.002 domains per annotated
   read measured across all arms). Candidate-pair recall and annotated-read agreement are different
   quantities, and only the second predicts what a pipeline user sees.
