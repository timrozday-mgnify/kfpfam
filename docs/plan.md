# kfpfam: a k-mer pre-filter for Pfam profile HMM search

A pre-filter that cuts the read × model comparisons `hmmsearch` needs to profile a metagenome with Pfam-A. Training data is the HMMs (or Pfam's own alignments); MGnify proteins are for evaluation only.

Status: exploration done on the laptop (2026-11-05); nothing built. Scripts and raw tables: `measure/explore/`.

## Verdict

**It should work, but not as first framed.** Three of the original ideas changed after measuring:

1. **Don't sample k-mers from the HMM; enumerate them.** At the same index size, enumerating by log-odds covers 2–3× as much of the model as sampling from emissions, with fewer background hits.
2. **k-mer hits alone can't give precision.** At thresholds that keep recall, one seed passes 9,000–26,000 of 30,134 models per read.
3. **The precision step should be an ungapped score along the seed's diagonal, not k-mer pairs or a rank rule.** The profile is already in memory. On the same seeds, this beats two-hit on both candidates and recall.

Measured on real data, relative to what `hmmsearch` reports on the same reads: **97.6–98.2% of hits (E ≤ 1e-3) and 98.8–99.1% of hits at GA, at about 13 candidate models per read**. That is about 2,200× fewer comparisons than 30,134 per read.

Comparisons stop being the cost at that point: the seed hits become the cost. I estimate a 10–50× wall-clock gain, and that is **not measured** (see Cost). The prototype has to measure it.

MMseqs2 already implements nearly this design and must be the first baseline (phase 1). It may make building our own unnecessary.

## Prior art

No paper found samples k-mers from Pfam emissions as a pre-filter (absence from my searches, not proof). The closest work selects words from the profile by score, or takes them from real family sequences.

| Work | What it does | Reported result | Relevance |
| --- | --- | --- | --- |
| [HMMERHEAD](https://pmc.ncbi.nlm.nih.gov/articles/PMC2931519/) (Johnson, Eddy, Portugaly 2010) | Ungapped 4-mers from match emission scores above a threshold; DFA scan; ungapped extension; two extended hits; gapped Viterbi | 99.7% of comparisons removed; about 20× (Forward); 2–4% of true positives lost | The direct ancestor, and short k like the results below. Built for full-length sequences against HMMER2, so two hits were affordable |
| [Sun & Buhler 2007](https://read.qxmd.com/read/17237102/designing-patterns-for-profile-hmm-search), [later version](https://manuscript.isc.ac/Inventory/49/1031736.htm) | Patterns designed from the HMM, maximising sensitivity under a false-positive bound | ≥ 90% sensitivity at 10–40×; profile filters about 100% at about 20× | The framing used here: budget false positives, maximise coverage |
| [MMseqs2 profile search](https://pmc.ncbi.nlm.nih.gov/articles/PMC6691333) | For each profile position, indexes all k-mers above a score (k = 5 here); two hits on a diagonal, ungapped, then alignment | 87% of HMMER3's hits at 95% precision; about 300× faster (full-length sequences, Pfam 31) | Same index as enumerating above T. **The main baseline** |
| [nail](https://pmc.ncbi.nlm.nih.gov/articles/PMC10862755/) (Roddy, Rich, Wheeler 2024) | MMseqs2 seeds (k = 6, relaxed `--k-score 80`), then sparse Forward/Backward | About 26× faster than HMMER3; most of the recall gap closed; the rest attributed to the k-mer stage | Evidence that seeding limits recall. Second baseline |
| [UProC](https://pmc.ncbi.nlm.nih.gov/articles/PMC4410661/) (Meinicke 2015) | 18-mers from Pfam-A **full alignment sequences**, nearest-neighbour words, per-residue scores | About 1000× faster than HMMER; TPR 85% against 50% (HMP, 100 bp); HMMER better from 150 bp on GOS and GNHM | Precedent for training on Pfam's sequences. A classifier, not a pre-filter. Its site says recent Pfam exceeds its limits |
| [NEAR](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12261438/) (Olson et al. 2025) | Learned per-residue embeddings, neighbour sharing | At least 5× faster than HMMER3's own pre-filter | The non-k-mer alternative |
| [HMMER3](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3197634/) (Eddy 2011) | MSV/SSV SIMD filter, P ≤ 0.02 | Reference | Default F1 passes about 2% of models, about 600 per sequence. NEAR's authors put MSV at about 70% of run time |
| [kMermaid](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12507277/) (2025), [Brinda et al. 2015](https://arxiv.org/abs/1502.06256) | Amino-acid k-mers to protein clusters; spaced seeds in classifiers | Spaced seeds improve the sensitivity/selectivity trade-off | Context for spaced shapes |

From memory, links not checked: gapped BLAST's two-hit method (Altschul 1997) is the origin of diagonal pairs; RPS-BLAST indexes PSSM neighbourhood words; HHblits pre-filters on column-state profiles; PRINTS fingerprints are ordered conserved motifs.

**Short reads are a known weak spot of `hmmsearch` itself.** UProC found HMMER at about 50% on 100 bp reads; [Techa-Angkoon, Sun & Lei 2017](https://pmc.ncbi.nlm.nih.gov/articles/PMC5657049) found half of RNA-Seq reads missed and nearly 70% of read pairs unaligned at GA. Measured below: 52% at GA on 150 bp, 29% on 100 bp. So the truth has to be defined carefully (see Decisions).

## What was measured

### Data

- **Models:** Pfam-A 38.2 (30,134 models; 4.75 M match states; median length 120).
- **Truth proteins:** 40,000 random proteins from `data/mgnify/mgy_clusters_head.faa`, full Pfam-A `hmmsearch --cut_ga` (26,014 domains, 4,867 families; 15,096 proteins with a hit).
- **Domain fragments:** 20,000 each of 50 aa (150 bp) and 33 aa (100 bp), overlapping a domain by ≥ 30 / 25 aa; each domain equally likely, so abundance-weighted.
- **Real reads:** 10,000 150-bp reads, `data/fmh-funprofiler/metagenome_example.fastq`, six frames, split at stops, ORFs ≥ 15 aa (mean 6.9 ORFs, 222 residues per read). They are simulated from reference genomes, so not independent of Pfam's sources.
- **Seed scoring:** a k-mer is "in the index" for a model when its log-odds (bits, against HMMER's background) reaches T at some window. So everything was scored by brute force across all of Pfam without building an index.

### The HMMs are deliberately flat

Mean relative entropy is **0.63 bits per column** (median 0.44; p90 1.30), set by `hmmbuild`'s entropy weighting (target 0.59). This is good for remote homology and bad for sampling: the mean column has about 11.8 effective letters, so 5-mers drawn from emissions rarely repeat and rarely discriminate.

### Sampling against enumerating (k = 5, 600 windows)

| Entries per window | Sampled from emissions: coverage / background | Top by log-odds: coverage / background |
| --- | --- | --- |
| 100 | 1.6% / 7.1e-5 | 4.8% / 4.3e-5 |
| 1,000 | 6.5% / 6.8e-4 | 13.3% / 4.1e-4 |
| 10,000 (8,600 distinct) | 20.7% / 5.9e-3 | 31.6% / 3.9e-3 |

Sampling is a noisy version of enumeration and costs more. It's only needed where exact enumeration is impossible, which it isn't here (branch and bound over a window, at most 20^k).

### Real sequences hit far more than the emissions predict

At k = 5, T = 10 a window covers only 3.9% of its emissions on average, yet 90% of the fragments `hmmsearch` detects hold such a seed. A read that clears `hmmsearch`'s bar must have high local score density, and a seed measures local score with the same function. This coupling is why the approach works: **recall relative to `hmmsearch` is much easier than recall of homology.**

Against all domain fragments, including those `hmmsearch` itself misses, k = 5, T = 10 reaches only 0.75.

### What `hmmsearch` reports on fragments of its own domains

| Fragment | True model reported, E ≤ 10 | E ≤ 1e-3 | ≥ GA |
| --- | --- | --- | --- |
| 50 aa (150 bp) | 0.80 | 0.64 | 0.52 |
| 33 aa (100 bp) | 0.62 | 0.41 | 0.29 |

On the real reads: 37% have a hit at E ≤ 1e-3 and 30% at GA, with about 1.1 models per hit read. So about 0.4 comparisons per read find anything, against 30,134 made.

### Seeds and read scores

Real 150-bp reads. Recall is over 3,000 read × model hits that `hmmsearch` reports on those reads' ORFs (E ≤ 1e-3 / ≥ GA). Candidates are mean models per read over all frames, whole Pfam, from 60 reads (30 for T = 7 and T = 9, marked †). "Diagonal ≥ X" is the best ungapped local score along any seeded diagonal.

| Seed | Index entries | Seed hits per read (est.) | One seed: candidates | Two hits, ±3: candidates, recall E-3 | Diagonal ≥ 22: candidates, recall E-3 / GA | Diagonal ≥ 25: candidates, recall E-3 / GA |
| --- | --- | --- | --- | --- | --- | --- |
| k4 T7 † | ~3.1e8 | 2.9e5 | 29,319 | 3,585, .93 | 13.4, **.976 / .988** | 2.7, .955 / .978 |
| k4 T8 | 1.1e8 | 9.2e4 | 23,728 | 582, .77 | 10.5, .920 / .955 | 2.1, .899 / .943 |
| k5 T8 | 3.2e9 | 1.4e5 | 26,041 | 830, .94 | 13.5, **.982 / .991** | 2.5, .961 / .982 |
| k5 T9 † | ~1.2e9 | 4.7e4 | 17,687 | 145, .80 | 10.4, .948 / .974 | 2.2, .928 / .964 |
| k5 T10 | 4.7e8 | 1.6e4 | 8,833 | 26, .63 | 8.3, .875 / .923 | 1.8, .855 / .911 |
| k5 T12 | 6.5e7 | 1.6e3 | 1,394 | 1.8, .30 | 4.3, .643 / .701 | 0.9, .630 / .692 |
| k6 T10 | 1.3e10 | 2.3e4 | 11,006 | 35, .81 | 10.5, .952 / .975 | 2.1, .932 / .966 |
| k6 T12 | 1.9e9 | 2.7e3 | 2,053 | 2.4, .50 | 6.5, .823 / .877 | 1.4, .804 / .865 |
| k7 T12 | 5.6e10 | 4.0e3 | 2,666 | 2.6, .65 | 8.3, .912 / .950 | 1.8, .892 / .938 |
| k7 T14 | 8.4e9 | 4.8e2 | 419 | 0.5, .39 | 5.0, .766 / .824 | 1.1, .750 / .815 |
| spaced `s6a` T10 | ~k6 T10 | ~k6 T10 | 11,763 | 34, .78 | 11.7, .965 / .980 | 2.2, .943 / .969 |
| `s6a` T12 | ~k6 T12 | ~k6 T12 | 2,106 | 1.7, .45 | 6.8, .837 / .887 | 1.5, .818 / .875 |

Entries and hits are from exact per-window tables on 400 models, scaled to Pfam; the T = 7 and T = 9 rows interpolate between neighbouring thresholds, and `s6a`'s are assumed equal to k6. `s6a` is the shape 11101011. Fragment recall tells the same story: at 33 aa, k5 T8 diagonal ≥ 22 gives .988 / .996; at 50 aa, .973 / .990 (`recall3.main.tsv`).

What this shows:

- **Short k is right, as suspected.** Index size per unit recall favours k = 4–5. k4 T7 matches k5 T8 at a tenth of the entries. Longer k trades index for fewer query hits (k7 T12 has 1/35 of k5 T8's hits for 18× the entries and less recall).
- **Diagonal score beats two-hit, at equal seeds.** k5 T8: 13.5 candidates at .98 against 830 at .94. A 50-aa read rarely has room for two non-overlapping seeds that each clear T.
- **Fixed spaced shapes help a little.** `s6a` T10 is .965 against k6's .952 at about equal size. A pair at a fixed distance is exactly a spaced seed, so this is also the measure of fixed-distance k-mer pairs indexed as one feature. Per-model column choice is untested and may do better (phase 3).
- **A flat threshold is unfair to some families.** Per family (1,400 with ≥ 3 detected fragments, 50 aa): one k5 T10 seed: mean .91, p10 .67, 0.5% lose everything. One k5 T8 seed: mean 1.00, p10 1.00. Two k5 T10 seeds: mean .63, p10 .18, 7.3% lose everything. All models do have seeds: 0.03% have none at k5 T12.

### The per-read rule (only models within 10 bits of the read's best)

From 150 fragments, 98 of them detected by `hmmsearch`, and 60 reads as above:

| Seed | Rule | Candidates mean (median) | True model kept, of detected |
| --- | --- | --- | --- |
| k5 T8 | Seeded, any score | 26,041 (26,438) | 1.00 |
| k5 T8 | Within 10 bits of best diagonal | 693 (70) | .94 |
| k5 T8 | Diagonal ≥ 22 | 13.5 (10) | .99 |
| k5 T8 | Diagonal ≥ 22 and within 10 bits | 8.7 (5) | not measured together; ≤ .94 |
| k5 T10 | Within 10 bits | 492 (117) | .88 |
| k5 T10 | Diagonal ≥ 22 and within 10 bits | 5.5 (4) | not measured together; ≤ .88 |

The rank rule **does work as a precision idea but loses to an absolute bar**: a read with only noise has a best score near 28, so everything within 10 bits passes. It adds about 1.5× on top of the bar, at about 5% of true hits (the best-minus-true gap has p95 of 12 bits at k5 T8). Those losses may be clan relatives outscoring the true family; unchecked, as no clan file was used. Scored on seed scores alone, a rank rule has no dynamic range. It may suit post-`hmmsearch` resolution better. It may also earn its place when made clan-aware (phase 4).

### Cost

| Item | Value | Source |
| --- | --- | --- |
| `hmmsearch` Pfam-A, 150-bp reads | 48 ms/read/core; **133 core-hours per 10 M reads** | Measured (`--cpu 3`, 10 k reads, 2.05e10 cells/s/core, 280 MB) |
| Candidates at diagonal ≥ 22 | About 13 models per read, 2,200× fewer | Measured, small sample |
| `hmmsearch` on those candidates | About 13 × 160 × 30 cells, negligible | Calculated, plus per-model overhead (phase 6) |
| Index | 4 B per entry: k4 T7 1.2 GB, k5 T9 5 GB, k5 T8 13 GB; profile as int8 95 MB | Calculated |
| Pre-filter at k5 T9 | 47 k hits × ~150 ns = **~7 ms/read** | **Estimated**, ns per hit is a guess |
| Pre-filter at k5 T10 | **~2.4 ms/read**, with recall .875 | **Estimated** |

So at high recall the seed-hit stage sets the time. Expect 10–50× overall (frame filtering, shorter-ORF cuts and X-drop rescoring all help; unfiltered reads carry 222 residues, most in wrong frames). Hundreds-fold would mean trading recall, as MMseqs2's 87% at 300× did. Measuring it is phase 5's purpose.

### Limits of these numbers

- Candidate counts come from 60 reads of one simulated set; 3,652 of 30,134 families were sampled by fragments, so rare families are unassessed.
- E-values used small Z (40 k–70 k), so a real run's E ≤ 1e-3 is a stricter bar and recall should be a bit higher.
- No sequencing errors, no frameshifts, no eukaryotic or viral sample, no low-complexity masking.
- This HMM file's Pfam release may not be MGnify's. Diagonal scores have no per-model calibration or length correction.

## Design this points to

```
reads → six frames, ORFs (frame filter from kfp's core)
      → k-mer lookup in {k-mer → (model, column)} with entries = all k-mers scoring ≥ T at that window
      → dedupe hits to diagonals; ungapped score on the profile, X-drop
      → (read, model) pairs with score ≥ bar (maybe calibrated per model)
      → hmmsearch per model on its reads, with -Z / --domZ set to the unfiltered search space
```

The output is then a subset of an unfiltered `hmmsearch`, with its scores and thresholds unchanged. The filter only decides which comparisons run.

## Plan

Every step goes to a Progress log at the end of this file. Anything big ships as Nextflow (SLURM + Singularity) for the user's HPC runs. Local runs: one job at a time, capped threads.

### Phase 1: baselines, before building anything

Do this first because it could end the project.

- `mmseqs` target-profile search of reads against Pfam, from the same Pfam release's alignments, sweeping `--k-score`, `-s`, and short-ORF settings (default extraction minimum lengths would drop these ORFs). `nail` too.
- Score both as pre-filters by the same recall and candidates measures, on the same reads.
- **Gate:** if MMseqs2 reaches ≥ 97% of `hmmsearch`'s GA hits at ≥ 20× overall speed, use it, wrap it in phase 6's hand-off and stop. If not, record what it misses and continue.

### Phase 2: truth at scale (HPC)

- Workflow `workflows/truth.nf`: proteins sampled from the full MGnify release with `mgy_proteins_pfam`, reads simulated with errors at 100 and 150 bp, background ORFs and at least one real metagenome per biome type. Unfiltered `hmmsearch` over all as truth, fixed `-Z`.
- Confirm the Pfam release behind `mgy_proteins_pfam` and use that HMM file. Earlier plans left this open.
- Enough reads that most families have ≥ 20 detected reads. Reuse `recall3.py` and `brute.py` scoring.
- **Output:** per-family and overall recall tables, the reference for all later steps.

### Phase 3: choose the seeds

No index needed; scoring by brute force stays valid on samples. In order of expected return:

1. **Per-model thresholds.** A budget per model, or a coverage target, against global T (per-family p10 is the target). The window's T as well: cap entries from very conserved windows, which take many entries for one read population.
2. **Real words as well as enumerated ones.** Add k-mers from `Pfam-A.seed`, then `Pfam-A.full`, never MGnify. This is UProC's source. It tests whether site correlation carried by real sequences adds recall that a flat, site-independent HMM can't. Also try thresholds from a non-flattened rebuild (`hmmbuild --enone` on seed) while scoring with the released HMM.
3. **Model-specific spaced seeds.** Column choice: the top-information columns within a span of 15–20, drawn from a library of 20–50 shapes a read can afford to look up. Sun & Buhler's approach in k-mer form, and the user's creative pairs at their strongest. Compare at equal index size with k4/k5.
4. **Pairs indexed as one feature** with gap tolerance (the `pairs.py` generator idea, with enumerated rather than observed k-mers). Only if 3 shows promise: entries grow with the square of each part's list.
5. **Reduced alphabets** (`alphabets.py`), low-complexity masking of reads and models.

**Output:** a chosen seed set per budget (index ≤ 2 / 8 / 16 GB), with recall and estimated hits.

### Phase 4: the read score

- Global diagonal bar against per-model calibration (the HMM file's `STATS LOCAL MSV` gives μ and λ, so a P-value is close to free), read-length correction.
- Read pairs: a model's best score over both mates, or the sum for concordant hits.
- A clan-aware rank rule (Pfam-A.clans.tsv), both as a filter and for final family resolution.
- Allowing one indel on the diagonal, if phase 3's losses show it matters.

### Phase 5: Rust prototype, and the speed gate

- `pfam_seed` builder (branch and bound per window; `Pfam-A.hmm` parser) and query (translation, FASTQ and frame filtering from kfp's `crates/core`, as a dependency or copied; direct-indexed table at k ≤ 5; sorted batch rescoring for cache).
- Measure reads/s/thread and peak memory against `hmmsearch` on the same machine. Rerun the chosen configuration over phase 2's truth.
- **Gate:** ≥ 97% of hits at GA, per-family p10 ≥ 0.85, ≥ 20× overall, index ≤ 16 GB. Below 10× even at relaxed recall, stop and report.

### Phase 6: hand-off to HMMER

- `pyhmmer` (per-model sequence blocks) against per-model FASTA + `hmmsearch`. `pyhmmer` means one new dependency; choose by measured overhead at 30 k models.
- `-Z`, `--domZ`, and `--seed` fixed so results match. **Check:** on a sample, filtered output equals unfiltered output restricted to candidates, score for score.
- `hmmscan` is not needed: both give the same scores, and grouping by model is the cheaper direction.

### Phase 7: does the profile survive?

The real question for profiling is not read recall. On phase 2's data:

- Family abundance against unfiltered: correlation, Bray-Curtis, per-family bias with depth, family detection at low depth.
- Bias by read position in the domain: seeds favour conserved blocks, so uneven read loss within a family may shift abundances.
- Benchmark against `kfp`'s Pfam roll-up, for context.

## Decisions

- **Truth (default until changed).** `hmmsearch` on the same read at GA, as this is what a pipeline would report. Domains a read overlaps (protein-level truth) is a far higher bar, which `hmmsearch` itself reaches only about half the time at 150 bp. A filter plus direct score could go beyond `hmmsearch` there (UProC's result), but that's a classifier with its own precision problem, and a different project.
- **Where it lives.** This repo (the user's choice, 2026-11-05), separate from `kfp`, whose scope is target abundances.
- **Tools run from Docker** for local runs (`docker/Dockerfile`). HPC runs use the same tools as Singularity through Nextflow.

## Reproducing

Run from `data/explore/` (holds `pfam.npz` and `truth/`), calling scripts by path, for example `python ../../measure/explore/recall3.py`:

1. `parse.py Pfam-A.hmm.gz pfam.npz`; `analytics.py pfam.npz` (`analytics.out`: entries, coverage and background per k and T).
2. `sample.py` then `hmmsearch --cut_ga` on the proteins (10 min at `--cpu 4`); `fragments.py`, then `hmmsearch -E 10` on the fragments; `fragtruth.py`.
3. `reads2orfs.py`, `hmmsearch -E 10` on the ORFs; `readtruth.py`.
4. `recall2.py` (by family traits; sampling against enumerating), `recall3.py` (`recall3.main.tsv`, `recall3.extra.tsv`), `floor.py`.
5. `mkreads.py`, `run_brute.sh` (whole-Pfam candidates, about 30 s per read, single thread), `summarize.py`.

## Progress log

- **Exploration (2026-11-05).** All of the above, on the laptop. Nothing committed to a design; phase 1 is next.
- **Measurement scan and explorer (2026-11-05, the user's request; repo moved to `kfpfam`).** `measure/scan` (Rust) records every seed of a sequence with what any policy could ask of it, across all of Pfam and with no index. `web/` evaluates policies from those seeds: enumerated, top-N, sampled or seed-alignment k-mers; hash and window sparsity; one seed, two seeds or a diagonal bar; a rank rule; linear adaptation of T and the bar on Pfam traits; budget and coverage targets. Truth: 23,439 hits in 7,269 families, from all 188,471 proteins in `mgy_clusters_head.faa` and the real reads, at Z = 20 million. A larger Z makes E ≤ 1e-3 a stricter truth, so recall at that bar reads higher than in the tables above; GA is unaffected. `checks/test_policy.mjs` holds the engine to the Python scan's GA recall within 0.05.
- **Phase 1, in progress (2026-11-05).** Tools run natively from Docker; the user's shell defaults to `linux/amd64`, so the platform must be forced. On 10,000 real reads at Z = 20 million: MMseqs2 prefilter on `Pfam-A.seed` profiles reaches GA recall 0.33 / 0.63 / 0.80 / 0.91 at 0.7 / 9.5 / 107 / 425 candidates per read (`s1`, `s4`, `s5.7` k5, `s7.5` k5). `s5.7` and `s7.5` at k = 6 ran out of memory at 14 GB. nail (Pfam HMMs, defaults): seeds GA recall 0.980 and E ≤ 1e-3 recall 0.992 at 2.4 candidates per read. Its time was about the same at 1,000 and 10,000 reads, so a fixed cost dominates.
- **nail at scale (2026-11-05).** MGnify fragments: GA recall 0.986 and E ≤ 1e-3 recall 0.993 at 1.3 seeded models per fragment; family p10 1.0 (3,644 families with 3+ hits). The same fragments among 1 M read-shaped sets of ORFs: GA 0.984, E ≤ 1e-3 0.992, so recall holds at that size. Wall time at 6 threads: 144 s at 100 k reads, 398 s at 1 M, 3.2 GB peak. That is about 2 min fixed plus 0.28 ms per read, or about 1.7 core-ms if all 6 threads stay busy, against 48 for `hmmsearch`: about 28×. **On this evidence nail passes phase 1's gate** (≥ 97%, ≥ 20×). Limits: the scale sets are one true-frame protein window plus shuffled ORFs per read, with 1.0 seeded models per read against 2.4 on the real reads, so real reads may cost more; recall at scale was tested on fragments, not reads; no sequencing errors.
- **Anchor window rule and family view (2026-11-05, the user's request).** The explorer gains a pseudo-alignment rule after k-mer selection: a chain of 2 or 3 non-overlapping anchors, in read order, each within ±w diagonals of the last, with no distance hashing, alone or combined with the diagonal bar. `scan bg` now also keeps weak-diagonal anchors that have a partner within ±6 diagonals, above a floor per seed (5-mers: 9.5 bits); below it, background candidates read low. On real reads at GA, 5-mers T 9.5: diagonal ≥ 22 gives recall 0.961 at 12.2 candidates; two anchors on one diagonal 0.772 at 10.3; window ±3 0.809 at 64; ±6 0.819 at 138; diagonal ≥ 25 or anchors ±3 with chain ≥ 24 gives 0.951 at 6.8, against 0.948 at 2.6 for diagonal ≥ 25 alone. So the window rescues about 4 points, plausibly indels, at about 6× the candidates, and adds nothing clear over moving the bar. Reads this short rarely hold two strong anchors. The family view shows one Pfam along its columns: the distribution of the score that selects k-mers with the policy's cut-off, k-mers kept per window, real seeds and anchors, column information, clustering figures, and each window's top k-mers. Same-spend-per-column (250) beat one global T at about equal index size: 0.984 against 0.975.
- **Emissions and co-occurrence (2026-11-05).** Match emissions are independent given the path, so pairs drawn or scored from the HMM carry no co-occurrence: a pair's log-odds is the sum of its parts, which the diagonal score already counts. Real co-occurrence can only come from sequences. `measure/cooccur.py` tests it on `Pfam-A.seed`; results pending.
