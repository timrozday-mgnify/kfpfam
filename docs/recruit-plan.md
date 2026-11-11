# Recruiting weak reads to confirmed targets

Goal: reads whose hit to a target is too weak to call at GA (too short, or too divergent) are recruited to
targets the sample already confirms. A model of how many weak reads a present target should give at its
depth sets the prior. Each recruited read goes to one target, or keeps an ambiguity set, under a threshold.
The output is candidate (read, target) lists for the follow-up sequence-alignment project.

Exploration: 2026-11-07, 50,000 ERR10889147 reads, hmmsearch truth (Z = 20 M, E ≤ 1000), reversed reads as
the null. Scripts: `measure/recruit_explore.py`, `measure/recruit_trunc.py`, `measure/recruit_model.py`;
outputs in `data/recruit/`.

## What the exploration found

**Restricting to confirmed targets removes most of the noise.** 1,975 targets have at least one read at GA
(8,106 reads). On reads with no GA call, weak pairs to those targets, real against reversed reads:

| E (Z = 20 M) | ≤ 0.1 | ≤ 1 | ≤ 10 | ≤ 100 | ≤ 1000 |
| --- | --- | --- | --- | --- | --- |
| real | 1,469 | 2,849 | 4,618 | 7,264 | 8,191 |
| reversed | 0 | 5 | 75 | 667 | 1,081 |
| FDR | 0 | .002 | .02 | .09 | .13 |

Against the other 28,000 targets, E ≤ 10 gives 3,936 real pairs against 1,458 reversed (FDR .37). The
confirmed set is about 15 times smaller, and the null shrinks with it.

**HMMER's E-values are optimistic on short ORFs.** Of the reversed pairs at E ≤ 1000, 7% have E ≤ 10. That is
7 times the expected rate. The null must be measured per sample, not taken from E.

**Depth model.** Excess weak reads per GA read (real minus reversed, E ≤ 10) by the target's GA depth:

| GA reads | 1 | 2–3 | 4–10 | 11–30 | > 30 |
| --- | --- | --- | --- | --- | --- |
| excess weak reads per GA read | 1.18 | .70 | .47 | .39 | .44 |

Above depth 3 the ratio is about 0.45 and flat. It is also flat with model length (.43–.59). The high ratio
at depth 1 is selection: a target seen once at GA is more often one whose strain sits far from the family, so
most of its reads are weak. So the prior is the target's GA count times about 0.45, with an allowance for
divergence.

**Ambiguity is rare and stays within clans.** Of 4,351 real weak reads (E ≤ 10), 95% hit one confirmed target.
Of the 197 that hit several, 181 stay within one clan. The best-minus-second score gap has a median of
3 bits, and 34% of gaps are under 2 bits.

**Model.** Two steps per read:

1. P(true), the probability that the read comes from any confirmed target. It uses a likelihood ratio by
   E bin, with the null from reversed reads and the true distribution from real minus reversed. The prior
   pi0 is reversed weak reads over real weak reads (0.15 here).
2. P(t | true) over the read's own targets, proportional to prior(t) × LR. The read goes to one target if the
   best P ≥ TAU. Otherwise it gets the smallest set whose P sums to TAU.

**Real reads.** 7,027 weak reads; P(true) cut, kept reads, reversed reads kept and FDR, at TAU 0.95 with the
depth prior:

| P(true) ≥ | kept: one target / set | reversed kept | FDR |
| --- | --- | --- | --- |
| 0.5 | 6,640 / 387 | 1,054 | .15 |
| 0.7 | 6,012 / 372 | 652 | .10 |
| 0.9 | 4,113 / 238 | 74 | .017 |

At FDR 1.7% this adds 4,351 reads to the 8,106 GA reads, which is +54%.

**Known origin.** Windows of 10–25 aa were cut from reads called at GA, plus as many windows from reversed
reads, and searched against the confirmed targets with pyhmmer (8 s for 66,000 windows × 1,975 models).

| window | any pair, E ≤ 1000 | of those at P(true) ≥ 0.9, TAU 0.8, depth prior: one right / set holds origin / wrong | reversed windows kept |
| --- | --- | --- | --- |
| 10 aa | .043 | .168 / .003 / .003 | 0 of 8,122 |
| 15 aa | .293 | .475 / .017 / .001 | 0 |
| 20 aa | .621 | .672 / .022 / .001 | 0 |
| 25 aa | .845 | .805 / .031 / .002 | 2 of 8,105 |

At P(true) ≥ 0.7 the kept share rises (20 aa: .909 / .030 / .003), and 13 of 8,116 reversed 20-aa windows are
kept. The depth prior cuts wrong assignments by 20–40% compared with a flat prior (25 aa at cut 0.5: .008
against .010; 10 aa: .044 against .069) and makes fewer sets. Below 15 aa HMMER barely reports anything, so
recruitment helps 20 aa and up. Windows cut from reads called at GA are the easy case: real weak reads are
also divergent, so treat these numbers as an upper bound on assignment accuracy.

**Targets with no GA read.** Unconfirmed targets with at least k weak reads at E ≤ 1, real against reversed:
k ≥ 2: 184 against 10; k ≥ 3: 100 against 4. Weak evidence summed over reads can confirm a target. This is
the route for the depth-1 targets where nail and Mode B lost most (see the depth breakdown in `plan.md`).

**Alignment position (`measure/recruit_reads.py`, `data/recruit/reads.txt`).** Each weak pair's strongest
domain is classed by where it sits:
- **cut:** the alignment runs to an ORF end while the HMM continues, so the read stops inside the domain;
- **inner:** the ORF continues past the alignment but the HMM is only partly covered;
- **whole:** the whole HMM is covered.

Weak pairs to confirmed targets, real / reversed, FDR in brackets:

| E | cut | inner |
| --- | --- | --- |
| ≤ 0.1 | 1,372 / 0 (0) | 95 / 0 (0) |
| 0.1–1 | 1,237 / 5 (.00) | 140 / 0 (0) |
| 1–10 | 1,476 / 42 (.03) | 292 / 28 (.10) |
| 10–100 | 1,846 / 292 (.16) | 799 / 300 (.38) |
| 100–1000 | 575 / 183 (.32) | 351 / 231 (.66) |

At equal E, 'cut' pairs have a third to a half of the FDR of 'inner' ones, and 95% of GA pairs are 'cut'.
Position adds to E rather than replacing it, because FDR still rises with E within each class. A weak hit
explained by the read ending inside the domain is credible. A read that spans more and still aligns only
partly is mostly noise.

**Reads with a hit, of 50,000:**

| | any pair, E ≤ 1000 | first hit at GA | + recruited, E ≤ 1 | + E ≤ 10 | + E ≤ 100 |
| --- | --- | --- | --- | --- | --- |
| hmmsearch, all ORFs | 27,710 | 8,106 | 2,785 (FDR .002) | 4,351 (.017) | 6,384 (.10) |
| Mode A, then HMMER | 26,164 | 8,085 | 2,716 (.002) | 4,175 (.018) | 5,938 (.11) |

12,577 reads (Mode A: 11,596) have weak pairs only to unconfirmed targets.

## Implementation plan

Inputs come from the Mode A hand-off: pyhmmer pairs at E ≤ 1000, for the real reads and the reversed reads,
with score, E, HMM and ORF coordinates. Everything after that runs in seconds.

1. **Hand-off output.** `measure/handoff.py` → a `kfpfam recruit-input` step that writes every pyhmmer pair at
   E ≤ 1000 with best-domain score, E, hmm_from/to, ali_from/to and ORF length, for real and reversed reads.
   Check: on the 50,000 reads its pairs match hmmsearch's at GA and weak bins (Mode A recall .97–.89).
2. **Reversed-read null, per sample.** Run the reversed reads through the same filter and hand-off. To save
   cost, test a 10–20% sample of reads, scaled. Check: LR by bin from the sample within 20% of the full set.
   Stratify the LR by ORF length (< 20, 20–35, > 35 aa) if that changes it, since the E-value error is a
   short-sequence effect.
3. **Confirmed targets.** At least one read at GA by default; the minimum GA reads is a parameter.
4. **Step 1, P(true) and an FDR cut.** Per-read P(true) from the best pair's LR, by E bin × position class
   (cut / inner / whole), each fitted on real against reversed reads. Then test a position-aware score
   feature: the observed score against the score a true member is expected to reach over the same HMM
   columns. The expected score comes from the model's per-column expected bits, or from the target's GA
   reads over those columns. A short read on low-information columns should not be penalised for its low
   score. The cut is set by an FDR
   target (reversed reads kept / real reads kept) rather than a fixed P, with a default of 2%.
5. **Step 2, assignment.** Prior = GA reads × rho, with rho fitted (about 0.45). Also test one EM update that
   adds recruited mass to the prior. The ambiguity threshold is TAU (default 0.8). Output is one target or
   a set, with P for each member. Report sets by clan too.
6. **Second route: confirm targets by summed weak evidence.** For unconfirmed targets, sum log LR over their
   weak reads. Confirm a target where the real sum beats the reversed-read distribution at the FDR target.
   Then run step 5 once more with the larger set.
7. **Depth check per target.** Expected weak reads = rho × GA reads. Flag targets far above it (a divergent
   strain, worth the alignment step) or far below it (repeats, or a GA call worth checking). Use a negative
   binomial with dispersion fitted across targets.
8. **Output for the alignment project.**
   - Per read: read, ORF, status (GA / one / set / none), P(true), and for each target: P, score, E, HMM
     and ORF coordinates.
   - Per target: GA reads, recruited reads (sum of P), ambiguous mass, depth flag.
   - FASTA per target of candidate reads, with weights.
9. **Evaluation.**
   - (a) Cut windows, as above, in a test.
   - (b) Simulated samples with known origin: MGnify proteins with Pfam hits, sampled to set depths (1–100)
     and mutated to 40–90% identity, as reads of 100 and 150 bp. These measure assignment accuracy on
     divergent reads, which cut windows cannot, and test the depth prior where depth is known.
   - (c) Reversed-read FDR, measured against its target on two or three other biomes.
   - (d) Report recruited reads, sets, and confirmed targets added, on 1.07 M reads.

Open questions:
- Should a set that spans clans be passed to alignment or dropped?
- Is the FDR target per read, or per target for step 6?
- Should GA reads with a near-tie second target also get a set? Today a GA call is final.

## Built (2026-11-07)

The code is in `recruit/`, with usage in the README:
- `pairs.py`: HMMER hand-off writing every pair with coordinates.
- `recruit.py`: `reverse` makes the null, `run` recruits.
- `test_recruit.py`: tests of the logic.

The evaluation scripts are `measure/recruit_sim.py` and `measure/recruit_eval.py`. Outputs are in
`data/recruit/{r50k,sim,r1m}/`.

The user's decisions:
- Sets that span clans go to alignment.
- The FDR is per read.
- Every read keeps its set, and `top` is the single best target.

**Hand-off.** `pairs.py` on Mode A's candidates (50,000 reads, 12 s) gives the same scores and
coordinates as hmmsearch on every shared pair. It finds .997 of hmmsearch's GA pairs.

**Likelihood ratio: logistic, not bins.** Binned features gave P(true) only a few values, so the FDR moved
in steps. Position, for example, could not show its worth at a 2% target. A weighted logistic fit of real
against reversed pairs gives a smooth P(true). The fit uses log E with hinges, position and its slope in
E, score against the expected score over the columns, and log ORF length. Since
h / f0 = odds × Nn / Nr, LR = π (odds − 1) / (1 − π). Reads recruited on the 50,000 reads, with no target
confirmation:

| features | FDR target 1% | 2% | 5% |
| --- | --- | --- | --- |
| E | 3,705 | 4,286 | 5,294 |
| E, position | 3,738 | 4,450 | 5,326 |
| E, position, score ratio | 3,883 | 4,472 | 5,382 |
| E, position, score ratio, ORF length (default) | 4,154 | 4,746 | 5,488 |

So position helps, but less than E: it adds 4% at 2%. All four features together add 11–12% at 1–2%.

**A fraction of the null is enough per read.** With 20% of the reversed reads, 4,881 reads are recruited
against 4,752 with all of them, and every one of the 4,752 is among them. 10% gives 5,024: about 5,000
null reads is too few for a stable estimate. The guide is about 10,000 null reads.

**Target confirmation counts each HMM region once.** Summing log LR over all of a target's reads let deep
proteins with repeats push reversed-read targets to sums over 100. On the simulated sample the cut then
rose to 120, and only 6 targets were added. Counting the best read per 15-column window of the HMM fixed
this. The simulation adds 172 targets, 94.8% of them true, which is the 5% target. On the 50,000 real
reads it adds 297 targets.

**Simulation with known origin.** 3,000 random MGnify proteins, 1,144 of them with a Pfam domain at GA.
Each protein gets an identity of 100–50% (BLOSUM62 substitutions) and 1–64 reads of 100 or 150 bp from
either strand, with six-frame ORFs: 37,541 reads in all. Mode A then pyhmmer, then `recruit.py` with the
defaults:
- **GA calls:** 2,457 reads, .987 right.
- **Recruited:** 2,802 reads. With one target, .940 right; in sets, .85 have the right top and .04 more
  hold it elsewhere in the set.

Reads overlapping a GA domain, called at GA / recruited right:

| identity | 1.0 | .9 | .8 | .7 | .6 | .5 |
| --- | --- | --- | --- | --- | --- | --- |
| GA | .396 | .274 | .219 | .079 | .044 | .014 |
| + recruited | .239 | .248 | .232 | .186 | .101 | .074 |

Recruitment adds more than GA calls give from 80% identity down, and more than doubles what is found at
70% and below. By protein depth, recruited reads rise from .09 at depth 1 to .19 at depth 31–64, because
the depth prior and target confirmation both need other reads.

**The reversed-read FDR is a lower bound.**

| FDR target | .002 | .005 | .01 | .02 | .05 |
| --- | --- | --- | --- | --- | --- |
| estimated | .001 | .005 | .010 | .019 | .050 |
| true (no true target in the set) | .016 | .038 | .044 | .056 | .069 |

At 100% identity the 54 false recruits split into two kinds. 25 are ORFs in other frames of real coding
DNA; all 23 false GA calls are of this kind too. 29 are the correct frame of a real protein hitting a
model the protein has no hit to at E ≤ 1. Reversing a protein reproduces neither, and false reads of
this kind come in runs from one region. The alignment step is what removes them. A coding-frame score
(the FGS frame, or the hexamer score from `measure/orf_rule.py`) could take out the first kind.

**TAU and the prior update.** Sets are rare on the simulation: 54 of 2,802 at TAU 0.8, 116 at 0.95.
Accuracy hardly changes with TAU (0.5–0.95) or the number of EM rounds (0–3). The defaults are TAU 0.8
and one round.

**1.07 M reads** (ERR10889147, Mode A, 6 threads):
- **Null cost:** the whole reversed set costs another filter run (124 s wall, 3.4 GB) and pyhmmer
  (107 s). `recruit.py` takes 26 s and 3.0 GB.
- **Pairs:** 998,202 real and 415,875 reversed at E ≤ 1000.
- **Targets:** 3,695 confirmed at GA, plus 345 from weak evidence.
- **Reads:** 172,605 called at GA, plus 108,401 recruited (+63%) at an estimated 2% FDR. Of the recruited
  reads, 4,191 keep a set.
- **Depth check:** rho is 0.58 recruited reads per GA read. 336 targets are flagged high and 6 low.
- **Null fraction:** 2% of the null (about 21,000 reads) recruits 108,879 reads, and 5% recruits 107,284.
  Both are within 1.5% of the whole null, so the rule of about 10,000 null reads holds at this depth.

**Loci by read position (2026-11-07, the user's request).** A read can hold two real coding frames, one
for each gene it crosses. `measure/real_orfs.py` now writes `OUT_PREFIX.coords.tsv`: each ORF's strand,
its first base on that strand, and the read length. Regenerating ERR10889147 gave byte-identical ORF
FASTAs, so the coordinates fit the existing pairs. With `--orf-coords`, `recruit.py` places each alignment
on the forward read. Pairs whose alignments share at least half the shorter one compete, in any frame.
Separate stretches of a read (two genes, or two domains) become separate loci. `loci.tsv` gains `read_from`
and `read_to`.
- **1.07 M reads:** the 350 reads with GA calls on two ORFs split into 346 whose alignments cover the same
  bases in different frames (never one clan: one frame is wrong) and 4 on separate stretches. With
  placement, kept loci go from 281,330 to 281,369. Reads with two or more kept loci go from 324 to 361,
  with 48 loci now on a second ORF where there were none.
- **Simulation from contigs:** `recruit_sim.py reads` now chains proteins into contigs of 2–6 genes, on
  random strands with 0–150 bp spacers, at a coverage of 0.3–30. That gave 95,922 reads, of which 1,951
  cross two genes by at least 30 bp each and 40 have GA domains in both. The truth is per gene segment, and
  `recruit_eval.py` judges a locus against the segments it covers. Placement changes almost nothing here:
  only 9 of the 180 domain-bearing segments in two-gene reads get a pair to their domain, because 102 of
  them are under 60 bp. Results with placement:
  - GA calls .999 right; recruited .960 right (one target) and .944 (sets).
  - True FDR .037 at an estimated .020.
  - Segments with a GA domain found: .29 overall, .55 at 100% identity, .06 at 50%.
- **Conclusion:** two genes in one read rarely give two hits, because the second segment is too short.
  Frames competing over the same bases is the common case, and placement now tells the two apart. It
  matters for the planned frame feature and wrong-frame null: the comparison is between frames over the
  same stretch, not across the whole read.

Not done:
- Other biomes (step 11), which needs the HPC Nextflow benchmark and a choice of samples.
- The coding-frame feature above.
- Per-read null fractions with target confirmation: confirming targets needs the whole null, and is
  skipped when the fraction is below 1.
