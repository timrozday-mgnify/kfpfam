# How kfpfam works now

As of `9e7bc16` (2026-11-24), in the setting called `top`. [method.md](method.md) describes the earlier contiguous 4-mer
modes (A and B) and what was tried and dropped; this page describes what runs today.

kfpfam is a prefilter. Given six-frame ORFs from reads, it picks the (ORF, model) pairs worth scoring. HMMER then scores those pairs with
the real HMMs at Pfam's gathering thresholds (GA). So kfpfam can lose a hit, but it can't change a score or add a false call.

On the slides it is called KAP, *K-mer-Anchored Profile filter*. "K-mer-set approximation of HMM profiles" describes the index
only. The search does not approximate a profile by its k-mers: they anchor the read on a model, and the model's own match scores then
score it.

## In one page

| Step | What happens | What is left, per ORF (`top`) |
|---|---|---|
| 1. Seed | Spaced 4-residue words shared by the ORF and a model's index | 278,000 hits |
| 2. Window gate | Ungapped int8 score of the 48 residues round the word, ≥ 14.5 bits | 1,225 survivors, 740 distinct diagonals |
| 3. Diagonal bar | Exact ungapped score of the whole diagonal, ≥ 15 bits | 469 |
| 4. Gapped bar | A cheap SIMD upper bound first, then an exact banded score, ≥ 21.25 bits (lower for short repeats) | 50 pass the bound, 7.2 pass the bar |
| 5. HMMER | `pyhmmer`, GA cutoffs, Z fixed | 0.04 |

Counts are per ORF over 16.3 M ORFs of 11 libraries, Pfam 39.0 (53,275 models). They predate the two fixes at the end of this page, which add about 6% to `top`'s
candidates.

## Before the search: building an index

Done once per Pfam release, from `Pfam-A.hmm` alone.

**1. Profiles (`measure/prep.py`).** Each match column becomes 20 log-odds scores in bits, emission over HMMER's background. Insert and
transition probabilities are not used. Output: `profile.f32`, `offs.u32`, `models.tsv`. An unknown residue (X) is scored by the filter as
the column's mean over the background, as HMMER has it.

**2. The seed.** A word is four residues taken from a stretch of eight, at positions 1, 3, 4 and 8: `10110001`. This is a Golomb ruler: no distance between two
positions repeats, so two overlapping words of one read share at most one residue and fail more nearly independently. Contiguous
words (`1111`) clump, and lose more at the same cost. It also catches period-three models such as collagen.

**3. A threshold per model (`filter seeds`, then `choose_seeds.py`).** A word is indexed at a column if its log-odds sum reaches T.

- For each model the filter works out how likely a 50-residue read that just reaches GA is to hold *no* indexed word, as a function of T. The read is drawn from the model's scores, tilted towards
  its high-scoring residues until the mean equals what GA needs, for each 50-column stretch.
- Each model gets the highest T with −ln P(no word), times a factor of 0.936 for words of one read overlapping, at least θ. `top` is θ = 9, so the model's estimate of a miss is e⁻⁹ ≈ 0.01%.
- T is therefore low for flat models and high for conserved ones: median 5.5 bits, 5–95% 4.7 to 6.4. One T for every model would cost more for the same loss.

**4. Enumeration (`filter build`).** For every window of every model, branch and bound lists all 4-residue words scoring ≥ T. They
are kept as (word → model, column). At 39.0 this is about 2.1 billion entries.

**5. Tiling (`filter tile`).** Models are grouped into tiles, whose int8 profiles fit in 2 MB (L2), and within a tile entries are grouped by word. An entry is then a u16 row of
the tile's profile, 2 bytes. 39.0 makes 157 tiles.

**6. Gapped bars (`recall_score.py bars`).** A model that is short and has a low GA2 can pass GA with few raw bits, as several copies add up. Each model gets its own gapped
bar where that is below the mode's: the weakest raw score GA2 allows, less 4 bits. 158 models are affected at 39.0.

## The search: `filter join`

ORFs are read in blocks of 400,000, so memory doesn't grow with the sample.

1. **Postings.** Every word of every ORF becomes a posting: the 48 residues round it, its ORF and position. Postings with the same word are sorted by where the ORF lies in the window and packed in groups of 32 lanes (16 on ARM).
2. **Join from the model's side.** Tiles run in parallel (rayon). For each word, each index entry of the tile is scored against each group of postings: 16 or 32 reads at a time against one column of one
   model. Four entries go together, as independent chains hide instruction latency. There are kernels for NEON, SSE4.1, AVX2 and AVX-512 VBMI, chosen at run time, and all give the scalar version's answer.
3. **Window gate.** The score is a running local maximum in int8 quarter bits. A window that already holds its whole ORF must reach the diagonal bar less 3 bits instead, since it is the whole diagonal to
   within rounding. Anything under 14.5 bits (or that) stops here. 99.6% of hits do.
4. **Deduplicate and rescore.** Survivors are sorted as 64-bit (model, ORF, diagonal) keys, repeats dropped, and the whole diagonal is scored in `f32`. Quarter-bit rounding would move about 5% of pairs across the bar, so the
   exact score decides.
5. **Every diagonal at the bar.** `top` scores a gapped band around *every* diagonal that reaches 15 bits, not only the best one (`FILTER_DIAGONALS=all`). The right alignment's diagonal isn't always the best ungapped one.
6. **Gapped score.** Affine local alignment in bits, opening 6 and extending 1, in a band of ±16 columns (`FILTER_BAND=16`). A SIMD upper bound goes first, four pairs at once. It can only overestimate, so it drops
   nine in ten pairs without losing any.
7. **Output.** Pairs at their gapped bar are written as (ORF, model index, diagonal score, gapped score). With `FILTER_BINARY=1` they are two u32,
   the ORF in its block and the model, after a header per block: 8 bytes a pair where the text takes about 50.

**HMMER.** In the pipeline the filter and HMMER are one task, `KFPFAM_HMMER`, joined by a pipe. `hmmer_pairs.py --stream` reads a block's pairs and that
block's ORFs, gives `pyhmmer` each model's candidates, and runs `--cut_ga` with Z set to the whole database while the next block is filtered. It writes a domtbl, as
`hmmsearch` would, and if asked the alignments as Parquet, as they are found. Nothing is written between the two. HMMER runs once: the filter's own gapped score is a cheap one, only to choose pairs.

## `top` against the other settings

| | `r97` | `r99` | `r999` | **`top`** |
|---|---|---|---|---|
| θ for the seed thresholds | 2.5 | 3.4 | 6.4 | **9** |
| Window gate, bits | 21.5 | 19.5 | 16.5 | **14.5** |
| Diagonal bar, bits | 22 | 20 | 17 | **15** |
| Band, and diagonals scored | 8, best | 8, best | 8, best | **16, all** |
| Gapped bar, bits | 28.25 | 26.25 | 23.25 | **21.25** |
| Median k-mers per model (38.2) | 2,200 | 4,200 | 16,000 | **34,000** |
| Index, 39.0 | 0.4 GB | 0.6 GB | 2.1 GB | **4.3 GB** |
| Filter peak memory | 2.3 GB | 2.6 GB | 4.6 GB | **7.8 GB** |
| Candidates per ORF | 0.15 | 0.33 | 1.9 | **7.2** |
| Pairs lost against `hmmsearch` | 2.4% | 0.80% | 0.041% | **0.002%** (13 of 715,654) |
| CPU-h per 50 M read pairs, filter + HMMER | 6 | 9 | 27 | **90** |

Pairs lost and CPU are from 11 libraries on Codon, Pfam 39.0 (`nar-review-response/response/B11_prefilter.md`). `top` costs about the same CPU as nail and keeps far more. `r999` is a third of its cost
at 0.04% lost. Every setting needs its own index.

## What it runs on

- **Filter:** Rust, about 2,100 lines in one file with its tests, and one dependency (rayon). About 1 MB.
- **Preparation and choosing thresholds:** Python and NumPy.
- **Final calls:** `pyhmmer`, in the same image and task.
- **Packaging:** a Docker image converted to Singularity on Codon (61 MB before pyhmmer and pyarrow were added). In RawR it is `--pfam_search_engine kfpfam`, process `KFPFAM_HMMER`, with `top`'s band and diagonals in `ext.args`.
- **Build cost at 39.0 (estimated on an M4 Pro from a tenth of 38.2's models):** reading the HMMs about 30 s, the threshold table about 70 core-s, `build` about 66 s on one thread, `tile` about 5 s. `tile` holds
  the index twice, so allow about 13 GB.

## What it loses, and what it doesn't promise

- **No filter of this kind is provably lossless.** A 30-bit alignment over 50 residues can hold no word above any T that still filters.
- **At `top`**, 13 of 715,654 pairs were lost on real reads, no family with 20 or more pairs falls below 95%, and it passes the parity test (at most 0.01% of pairs). It was tuned knowing an earlier run's misses on these libraries, so it's not an independent test.
- **Known misses:** a model shorter than the seed's 8 columns (WWamide, PF08258) is never indexed. Others are seeds just under T, and alignments outside the band.
- **Two fixes since the recall test:** X is scored as HMMER scores it, and short repeats have their own gapped bar. `B11`'s arms and the 715,654-pair figures include both; the per-ORF counts above do not.
- **Not measured:** Pfam 37.4 (to be run), reads over about 150 bp in the pipeline, and a release-matched `hmmsearch` at full depth.

## Running `top` by hand

```bash
# once per release; Pfam-A.hmm.gz and an accession, GA1, GA2 table (recall_test.sh writes ga.tsv)
python3 measure/prep.py Pfam-A.hmm.gz Pfam-A.clans.tsv.gz prep
filter seeds prep ga.tsv seeds.f32 10110001
python3 choose_seeds.py seeds.f32 10110001 9 g9.T.txt
FILTER_T_FILE=g9.T.txt filter build prep 10110001 0 g9.idx
filter tile prep g9.idx g9.tidx
python3 recall_score.py bars --prep prep --ga ga.tsv --out prep/gapped.txt > weak_models.tsv

# per sample: diagonal bar 15, window gate 14.5, blocks of 400,000 ORFs, one line per sequence
FILTER_MIN_GAPPED=21.25 FILTER_GAPPED_FILE=prep/gapped.txt FILTER_BAND=16 FILTER_DIAGONALS=all \
  filter join prep g9.tidx 15 14.5 400000 orfs.faa pairs.tsv
```

`choose_seeds.py` and `recall_score.py` are in `nar-review-response/hpc/B11_pfam_arms`, not in this repository yet.
