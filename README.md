# kfpfam

Pfam annotation of metagenomes with k-mers: a pre-filter that cuts the read × model comparisons
`hmmsearch` needs, so a metagenome can be profiled with Pfam-A at a fraction of the cost. `hmmsearch`
still makes every call, so its scores and thresholds are unchanged.

Early stage: a prototype filter, measurements and an explorer. [The method](docs/method.md) summarises where it stands; see [the plan](docs/plan.md) for the literature, the measurements so far
and the phases. In short, on 150 bp reads, k-mers enumerated from the HMMs plus an ungapped score along
the seed's diagonal keep about 98% of `hmmsearch`'s hits at about 13 candidate models per read, out
of 30,134. Speed is in the Filter section below.

Training uses the HMMs and Pfam's own alignments only. MGnify proteins are for evaluation.

## Tools

`hmmsearch`, MMseqs2 and nail run from one image:

```bash
docker build --platform linux/arm64 -t kfpfam-tools docker/  # on Apple silicon, or an amd64 image runs emulated
```

Local data goes in `data/` (not tracked). Pfam-A comes from
[the Pfam FTP site](https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/).

## Filter

`filter/` is the pre-filter in Rust: an index of every 4-mer scoring at least 7 bits at a column, joined
to blocks of reads from the model's side with SIMD (NEON, or SSE4.1 on x86), then an exact diagonal and
a banded gapped score on what passes. On one human gut metagenome (1.07 M reads, M4 Pro, 6 threads) it
takes about 0.7 core-ms per read, and with HMMER on its candidates about 0.9 (Mode A's old settings), against 1.2–2.0 for nail,
with higher recall below GA. See the plan's Speed plan and Progress log.

```bash
(cd filter && cargo build --release && cargo test --release)
filter/target/release/filter build data/prep 4 7 data/index/k4T7.idx
filter/target/release/filter tile data/prep data/index/k4T7.idx data/index/k4T7.tidx
# bar 22 bits on the diagonal, window gate 14.5 bits, blocks of 500,000 ORFs; out: ORF, model, diagonal, gapped
filter/target/release/filter join data/prep data/index/k4T7.tidx 22 14.5 500000 orfs.faa pairs.tsv
measure/speed.sh NAME 4 join data/index/k4T7.tidx 22 14.5 500000  # timed and scored on the 50,000-read truth
```

Two settings against nail on the same machine, from the plan's Progress log (1.07 M reads; recall on 50,000;
Mode A re-tuned 2026-11-11 after HMMER moved to worker processes, Mode B's HMMER cost is from before):

| Mode | Index | Diagonal bar, gate | HMMER on gapped ≥ | core-ms per read | Recall ≥GA; −5–GA, −10–−5, −20–−10 |
| --- | --- | --- | --- | --- | --- |
| A: nail's cost, higher recall | `build data/prep 4 6.5` | 17, 16.5 | 21.25 | 1.20 | .999; .989 .963 .952 |
| B: nail's recall, less cost | `build data/prep 4 7.5` | 20, 16.5 | 26.25 | 0.37 | .987; .898 .655 .385 |
| nail defaults | | | | 1.24 | .987; .85 .54 .33 |

## Recruiting weak reads

`recruit/` takes reads whose HMMER hit is too weak to call at GA and assigns them to targets the sample
confirms. Each read goes to one target, or to a set of targets when it is ambiguous. The output is a
list of candidates for a later alignment step. Reversed ORFs give the null. See
[the recruitment plan](docs/recruit-plan.md) for the model and its measurements.

```bash
python3 recruit/recruit.py reverse orfs.faa decoy.faa            # the null: reversed ORFs
# Mode A; FILTER_MIN_GAPPED writes only pairs HMMER will see (diagonal 17 alone is ~300 M rows per 1 M reads)
export FILTER_MIN_GAPPED=21.25
filter/target/release/filter join data/prep data/index/k4T6.5.tidx 17 16.5 400000 orfs.faa real.filter.tsv
filter/target/release/filter join data/prep data/index/k4T6.5.tidx 17 16.5 400000 decoy.faa decoy.filter.tsv
# HMMER on the candidates (pyhmmer container): every pair at E <= 1000 with coordinates
python3 recruit/pairs.py data/pfam/Pfam-A.h3m data/prep/models.tsv orfs.faa real.filter.tsv 21.25 6 real.pairs.tsv
python3 recruit/pairs.py data/pfam/Pfam-A.h3m data/prep/models.tsv decoy.faa decoy.filter.tsv 21.25 6 decoy.pairs.tsv
python3 recruit/recruit.py run --real real.pairs.tsv --null decoy.pairs.tsv --models data/prep/models.tsv \
  --prep data/prep --orfs orfs.faa --orf-coords orfs.coords.tsv --out recruit_out
  # --fdr 0.02 per read, --tau 0.8, --target-fdr 0.05; orfs.coords.tsv from measure/real_orfs.py
python3 recruit/test_recruit.py
```

- **`loci.tsv`:** one row per locus: the alignments over one stretch of a read, in any frame. A read
  crossing two genes, or one with two domains, has two loci. Without `--orf-coords`, only alignments that
  don't overlap on one ORF split a read.
  - `status`: `ga`, `one`, `set` or `none`.
  - `p_true`, `q`: the read's probability of coming from any confirmed target, and its q-value.
  - `top`, `top_p`: the best target, for when no later analysis will use the set.
  - `set`, `clans`: the set's targets with their probabilities, and their clans.
  - The top pair's score, E-value, coordinates (HMM, ORF and `read_from`/`read_to` on the read), position
    class and score ratio.
- **`hits.tsv`:** one row per (locus, target in its set), in pfam-reannotate's hits columns with the
  ORF as `query`, including the CIGAR. This is the input to the Pfam-anchored MGnify search. The columns
  are in [docs/recruit.md](docs/recruit.md).
- **`targets.tsv`:** per target.
  - How it was confirmed (`ga`, or `weak` from summed weak evidence).
  - GA reads, and recruited reads as both counts and probability mass.
  - The depth check against the expected recruited reads.
- **`candidates.faa`:** each kept ORF once per target in its set. The header carries `target`, `status`,
  `weight` and `top`, so a later step can keep only `top=1`.

A fraction of the reversed reads is enough for the per-read step: about 10,000 reads, or
`recruit.py reverse orfs.faa decoy.faa FRACTION` with `--null-fraction FRACTION`. Confirming targets from
weak evidence needs the whole null, and is skipped when the fraction is below 1. The reversed-read FDR is
a lower bound: on the simulated sample the true FDR was 3–4 points higher (see the plan).

## Explorer

`web/` compares k-mer policies: how k-mers are drawn from the HMMs, sparsity, the read hit rule and the
information bar, and how each adapts per Pfam. Every number is recomputed in the browser from seeds
measured across all of Pfam. Its data is not tracked (136 MB). To rebuild it, with Pfam-A, its clans and
seed alignments in `data/pfam/`, the exploration sets in `data/explore/` and MGnify proteins:

```bash
(cd measure/scan && cargo build --release)
python measure/prep.py data/pfam/Pfam-A.hmm.gz data/pfam/Pfam-A.clans.tsv.gz data/prep
ln -s ../pfam/Pfam-A.seed.gz data/prep/seed.gz
PY=python measure/run_truth.sh
python measure/export.py data/scan data/prep data/scan/abundance.tsv web/data
# Optional: the family view's secondary structure tracks. See the plan's Progress log for data/ss.
PYTHONPATH=measure python measure/ss_accuracy.py data/ss data/prep data/ss/out
PYTHONPATH=measure python measure/ss_kmers.py web data/prep data/ss/out web/data
python3 -m http.server 8791 --bind 127.0.0.1 -d web  # then open http://127.0.0.1:8791/
```

`web/data/baselines.json` holds the baseline points, from `measure/baseline_score.py`.

```bash
node checks/test_policy.mjs  # policy engine against the separate Python scan
node checks/compare_rules.mjs  # read hit rules side by side
```

Related: [kmer-functional-profiler](https://github.com/timrozday-mgnify/kmer-functional-profiler),
[kfp-explore](https://github.com/timrozday-mgnify/kfp-explore).

Licence: GPL-3.0-or-later.
