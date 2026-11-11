# kfpfam: method

kfpfam is a pre-filter for Pfam-A annotation of metagenome reads. It decides which (ORF, model) pairs
`hmmsearch` should score, then hands only those to HMMER. HMMER still makes every call, so scores,
E-values and gathering thresholds are those of an unfiltered `hmmsearch`; the filter only changes which
comparisons run.

The details, every measurement and its date are in [the plan](plan.md)'s Progress log. This page
summarises where things stand. All numbers are from one human gut metagenome, `ERR10889147`
(1.07 M reads, 6.63 M six-frame ORFs, mean 30 aa), on an Apple M4 Pro, unless stated.

## In short

| | Mode A: nail's cost, higher recall | Mode B: nail's recall, less cost | nail defaults |
| --- | --- | --- | --- |
| Index | k4 T6.5, 466 M entries, 0.95 GB | k4 T7.5, 168 M entries, 0.40 GB | |
| Diagonal bar, window gate | 17, 16.5 bits | 20, 16.5 bits | |
| HMMER on pairs with gapped score ≥ | 21.25 bits | 26.25 bits | |
| CPU, core-ms per read (filter + HMMER) | **1.20** (0.875 + 0.326) | **0.37** (0.28 + 0.09) | 1.24 |
| Wall, 1.07 M reads, 6 threads | 237 s | 71 s | 349 s |
| Peak memory | 3.5 GB filter, 5.7 GB HMMER (both flat in sample size) | 5.5 GB | 3.0 GB |
| Recall at GA | .999 | .987 | .987 |
| Recall at GA − 5–GA, −10–−5, −20–−10 | .989 .963 .952 | .898 .655 .385 | .85 .54 .33 |

Recall is against `hmmsearch` on all ORFs of 50,000 of the reads (Z = 20 M, E ≤ 1000), binned by the
pair's best domain score minus the family's GA. Weak bins carry noise: on reversed reads, precision is
about 96%, 65% and 36% in the three weak bins.

## The method

```
ORFs (six frames, ≥ 15 aa)
  → block of ~400,000–1,000,000 ORFs: postings, one per 4-mer occurrence, each with the 48 residues
    around it, grouped 16 at a time by 4-mer
  → join from the model's side, tile by tile: every index entry (model column, 4-mer) scores the
    4-mer's postings against the 48 profile columns around it, 16 lanes per NEON instruction
  → window gate: best local score in the window ≥ gate
  → exact ungapped diagonal (f32) on survivors ≥ diagonal bar
  → banded gapped score (±8 columns, open 6, extend 1 bits)
  → pyhmmer on each model's candidates with gapped score ≥ the HMMER bar, Z fixed
```

**Index.** For every window of k = 4 columns in every model, all 4-mers whose log-odds sum (bits,
against HMMER's background) reaches T, found by branch and bound. Short k with a low T is the point:
it catches the weak, short matches 150-bp reads give, and an enumerated index covers 2–3× as much of a
model as words sampled from its emissions. The index is stored in tiles: runs of whole models whose
padded int8 profile has at most 65,535 rows (2 MB, so it stays in L2). Within a tile, entries are
grouped by 4-mer, and an entry is a u16 row of the tile's profile. `filter build`, then `filter tile`.

**Postings.** Each block of ORFs is turned into postings: for every 4-mer occurrence, the 48 residues
centred on it (padding outside the ORF), its ORF and position. A 4-mer's postings are sorted by where
the ORF sits in the window and laid out 16 lanes at a time, offset-major, with a per-group record of the
offsets that hold a residue in some lane and of the lanes whose window holds the whole ORF.

**Join.** Tiles run in parallel. Within a tile, each 4-mer's postings list is read from memory once and
scored against every entry of the tile that holds the 4-mer. The kernel is NEON `tbl`: each step looks
up 16 residues in a profile column's 32-byte row, then a saturating add and two maxes give a running
local best. Four (row, group) jobs share one batch to hide latency; a batch is one group against four
rows, so they share the group's offset span. A hit costs about 1.4 ns scored and 1.6 ns all in (k4 T7);
the scalar version is kept as the oracle and gives identical output.

**Window gate, then the diagonal.** The window score is int8 in quarter bits. Survivors pass a gate of
14.5–16.5 bits. A window that holds its whole ORF already scores the whole diagonal to within rounding,
so it must reach the bar less 3 bits instead. Survivors are deduplicated by (model, ORF, diagonal) and
rescored in f32 along the whole diagonal, as rounding to quarter bits moves about 5% of pairs across the
bar. Pairs at the diagonal bar get a gapped score in a band of ±8 columns around their best diagonal.

**HMMER hand-off.** pyhmmer in fork worker processes (`recruit/pairs.py` `search_each`), each model's optimized profile against its candidate ORFs,
Z fixed, E ≤ 1000 or `--cut_ga`. Best-domain scores agree with `hmmsearch` within 0.05 bits on every
shared pair. HMMs are read from a binary file (`Pfam-A.h3m`, 886 MB, 1 s to load against 224 s as
text). For GA-only output the HMMER bar can rise to gapped ≥ 28.25 with GA recall .9921 against .9923
at 24.25, at a fifth of the candidates.

**Parameters.** T sets index size and hits; the diagonal bar sets which pairs are rescored with gaps;
the gate trades survivors against pairs lost (0.18% of the reference's pairs at 14.5 bits, 0.81% at
16.5, recall by bin unchanged); the HMMER bar sets candidates per read. Block size trades memory for
speed (400,000 ORFs: 3.06 GB; 1 M: about 5 GB, slightly faster). Lowering the diagonal bar from 22 to 20
is the cheapest recall: gaps lift pairs over the gapped bar.

## Compared with nail

[nail](https://pmc.ncbi.nlm.nih.gov/articles/PMC10862755/) (Roddy, Rich and Wheeler 2024) seeds with
MMseqs2 and then runs sparse Forward/Backward in a "cloud" around each seed. Its seeding is the
pre-filter, and on this sample it is most of its time: 307 of 349 s wall at 1.07 M reads.

| | nail's pre-filter (MMseqs2 seeding) | kfpfam's pre-filter |
| --- | --- | --- |
| Words | MMseqs2's similar-k-mer neighbourhood of each profile position, k = 6 with a relaxed k-score (from the paper) | Every 4-mer (or other k) scoring ≥ T in each window of the HMM, enumerated exactly |
| Hit rule | MMseqs2's diagonal rule, then ungapped and gapped alignment to give seeds | One word hit, then a 48-residue window score, the exact ungapped diagonal and a banded gapped score |
| Join | k-mer-ordered index lookups per query | Model-side join over a block's postings, tiled to L2, 16 lanes per SIMD op |
| Candidates per read | 1.5 seeded models (defaults), 64.5 at `-S 1e-2` | 19.3 (Mode A) or 1.2 (Mode B) |
| Recall of candidates (GA; −5–GA, −10–−5, −20–−10) | defaults .987; .85 .54 .33 (reported); `-S 1e-2` seeds .989; .908 .832 .760 | Mode A .999; .989 .963 .952; Mode B .987; .898 .655 .385 |
| Pre-filter cost | seeding 1.7 core-ms per read (wall × threads) at 1.07 M reads; nail as a whole 7.4 at 50,000 reads, as about 2 min is fixed | 0.28–0.875 CPU core-ms per read; about 1–2 s fixed (index and HMM load); about the same per read at 50,000 reads |
| After the pre-filter | nail's own sparse Forward/Backward and scores (30 s here) | HMMER itself on candidates (0.09–0.33 core-ms per read) |
| Memory and disk | 3.0 GB; no index beyond MMseqs2's temporary files, which run to many GB | 3–5.5 GB; index 0.25–1.5 GB on disk |

Why kfpfam's pre-filter keeps more weak hits at less cost:

- **One short seed, scored.** A 150-bp read gives about 50 aa; it rarely holds two independent
  k-mer hits on one diagonal, so two-hit rules lose weak hits (k5 T8: two hits ±3 kept .94 at 830
  candidates per read; a diagonal bar kept .98 at 13.5). kfpfam takes every single hit and scores it.
- **Shorter words at a lower threshold.** 4-mers at T 6.5–7.5 bits against MMseqs2's 6-mers: a weak
  match on a 50-aa read is more likely to hold a short high-scoring word than a long one, and the
  cost of the extra hits is paid in the join, not in candidates.
- **Cost per hit, not hits.** k4 T7 makes about 220,000 hits per read; kfpfam spends about 1.6 ns on
  each, against about 130 ns in the first, k-mer-ordered version. The join reads each postings list
  once per tile and keeps the profile in cache.
- **No fixed start-up cost.** nail's 2 minutes make it about 7× dearer per read at 50,000 reads.

What nail does that kfpfam doesn't: it aligns and scores in one tool, and its sparse Forward/Backward
is fast once seeded. Its scores differ slightly from HMMER's: on 50,000 merged fragments, 1.67% of GA
pairs were called by one and not the other, so the raw-reads pipeline keeps HMMER for the calls.

## Features tried and not kept

Each was measured on the same truth; the plan's Progress log has the numbers and dates.

**Seeds and the index**

| Feature | Result | Why not kept |
| --- | --- | --- |
| k-mers sampled from emissions | At equal index size, half to a third of the coverage of enumerated words, more background | Enumeration is exact and cheap |
| Longer k (k5 T9) as default | About 6× fewer hits, but GA recall .988 → .974 | The user's goal is recall at nail's speed |
| Fixed spaced seeds | `s6a` .965 against k6's .952 at about equal size; 3-g-3 (two 3-residue blocks 1–2 apart) 0.12× hits at k4 T7's recall | Needs about 50× the index, and long keys make the join's lists too short to fill SIMD lanes |
| Spaced seeds with a drift window (±1–2 residues in the gap) | The window never helped, in any structure class or at any hit level | Its chance hits cost more T than the indels it rescues |
| k or T by secondary structure (GOR IV/V on emissions, DSSP) | Per-class choice 0.8× hits at 2.5× index; indel-rate terciles did as well (0.74×) | Hits fall alike in helix, strand and coil as k grows; structure adds nothing over conservation |
| Per-profile T, fitted to each model's hits | In-sample 0.69× hits at no loss; held out it loses (GA .913) | Most models have one or two hits per sample |
| Per-model hit cap (a length correction) | About 1.2× fewer hits at equal recall; at 0.8× hits, GA .992 → .990 | Trades GA recall; kept as an option |
| Per-profile seed shape (k4–k7, spaced) | Brute force: 0.75× hits at 0.93× index, no recall loss | In the SIMD join, k5 lists fill 35% of a 16-lane group (k6 14%) against 92% for k4, so it costs 1.08–1.23× k4 T7's lanes |
| Promiscuous k-mers: longest index lists | 14% of hits for −.006 to −.013 in the weak bins | Lists are not skewed (mean 1,754, p99 4,579); kept as an option |
| Promiscuous k-mers: list length × background | Worse than the plain list at every share | Common-composition k-mers carry real signal |
| Promiscuous k-mers: frequent in this sample | 30% of hits within the gate when scored on the reads that built it | A list from one block of 75 reads doesn't carry to the next; needs large blocks |
| Same spend per column | .984 against .975 at equal index size | Not adopted: a small gain that changes the k4 T7 index the user chose; untested since the join was built |
| Whole-sample k-mer sketch | At most 1.2–1.5× fewer models at 90% of present models | Lumping a sample makes chance matches reach 30–34 bits at any k |
| Words from `Pfam-A.seed` co-occurrence | 1.6× over expected in true hits, 0.9 in background; as a filter 58–66% recall | Too weak as a filter; maybe as evidence on weak hits |

**Read-side rules**

| Feature | Result | Why not kept |
| --- | --- | --- |
| Gate on the read's own secondary structure | 0.56× hits, no better than raising T (0.58×) | Same trade as T |
| Predicting informative zones from the read (GOR-style window, 4-mer tables, seed rate per hit) | At best matches raising T at light cuts, loses beyond | Likely true seeds are the promiscuous k-mers (r .55–.77 with hits), and a read-side rule can't see the column |
| Scoring repeated windows once | 5.5% of (k-mer, window) postings repeat in a block | Not worth the bookkeeping |
| ORF rule: FGS frame plus others with hexamer score ≥ −5 | 0.55× hits; diagonal ≥ 22 recall .989; .915 .761 .648 against .992; .936 .875 .852 | Costs weak hits; it is the raw-reads pipeline's rule, so an option, and Mode B's alternative |

**Precision rules**

| Feature | Result | Why not kept |
| --- | --- | --- |
| Two hits on a diagonal | .94 at 830 candidates against .98 at 13.5 for a diagonal bar (k5 T8) | Short reads rarely hold two seeds |
| Anchor chains (2–3 anchors within ±w diagonals) | About 4 points of recall at about 6× the candidates | No better than moving the bar |
| Rank rule (models within 10 bits of the read's best) | 693 candidates at .94 | A noise-only read's best is about 28 bits, so everything passes |
| Pre-computed upper bounds per hit | Remove 0.7–1.3% of hits | Columns' best scores are too high to bound a short read |
| Hit counts in place of scores | Two hits within ±3 at k4 T7 gave .93 at 3,585 candidates | Scoring is cheaper per unit of recall |
| Whole-matrix gapped score | 9.5 core-s per 50,000 reads against 0.6 for ±8 columns, recall within .001 | Replaced by the band |

**Engineering**

| Feature | Result | Why not kept |
| --- | --- | --- |
| k-mer-ordered index, scalar, whole diagonal per hit | 130 ns per hit, 29 core-ms per read | Each hit jumps to a new model's profile |
| Window clamped into the ORF | Needs 32 residues at 1.08% surviving, but no fixed SIMD layout | Centred 48-residue windows: 0.94% surviving |
| Quarter-bit integer diagonals | Move 4.6% of pairs across the bar at 22 | The exact diagonal stays in f32; only the gate is int8 |
| Model-side join without tiles | 7.4 ns per hit with NEON, the kernel alone 0.9 | Each entry fetched its postings list from memory: tiles fixed it (1.4 ns) |
| Software prefetch of postings lists | 4.9 → 3.4 ns per hit | Superseded by tiles |
| Lane bitmasks in scalar code | No change | Superseded by masks from NEON compares |
| One `hmmsearch` call per model in pyhmmer | 5.7 core-ms per read, thread start-up | One `Pipeline` per worker thread: 0.25–0.56; fork workers instead: 2.7× less (threads lose most of their time to the interpreter lock) |
| Gathering cut-offs for speed | Same cost as E ≤ 1000 | No gain; used only for output |
| HMMs as text | 224 s to load | Binary file, 1 s |
| AVX2 | Not written | No x86 toolchain on the laptop; the SSE4.1 path is written but not yet built |

## After HMMER: recruiting weak reads

HMMER's pairs on the filter's candidates go on to [recruitment](recruit.md). Reads too weak for GA are
assigned to families that the sample confirms, with ambiguity kept as sets. Each kept hit's alignment is
written for the Pfam-anchored search of MGnify proteins. On the 1.07 M reads, recruitment adds 108,750
loci to 172,619 at GA, with an estimated per-read FDR of 2%. On a simulated sample the true FDR was
3.7%.

## Limits

- One sample, one machine. Other biomes, x86 (SSE4.1 untested, AVX2 not written) and the full
  read set's truth are step 7 of the plan: a Nextflow benchmark for the HPC.
- Truth is `hmmsearch` on the same ORFs. Reads where `hmmsearch` itself misses a domain don't count.
- Recall is per (ORF, model) pair; weak bins include noise (see reversed-read precision above).
- nail ran in Docker; pyhmmer in Docker; the filter natively.
- The diagonal bar is not calibrated per model.
