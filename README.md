# kfpfam

Pfam annotation of metagenomes with k-mers: a pre-filter that cuts the read × model comparisons
`hmmsearch` needs, so a metagenome can be profiled with Pfam-A at a fraction of the cost. `hmmsearch`
still makes every call, so its scores and thresholds are unchanged.

Early stage: no tool yet, only measurements and an explorer. See [the plan](docs/plan.md) for the literature, the measurements so far
and the phases. In short, on 150 bp reads, k-mers enumerated from the HMMs plus an ungapped score along
the seed's diagonal keep about 98% of `hmmsearch`'s hits at about 13 candidate models per read, out
of 30,134. Wall-clock speed-up is not measured yet.

Training uses the HMMs and Pfam's own alignments only. MGnify proteins are for evaluation.

## Tools

`hmmsearch`, MMseqs2 and nail run from one image:

```bash
docker build --platform linux/arm64 -t kfpfam-tools docker/  # on Apple silicon, or an amd64 image runs emulated
```

Local data goes in `data/` (not tracked). Pfam-A comes from
[the Pfam FTP site](https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/).

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
