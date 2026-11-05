# kfpfam

Pfam annotation of metagenomes with k-mers: a pre-filter that cuts the read × model comparisons
`hmmsearch` needs, so a metagenome can be profiled with Pfam-A at a fraction of the cost. `hmmsearch`
still makes every call, so its scores and thresholds are unchanged.

Early stage: nothing is built yet. See [the plan](docs/plan.md) for the literature, the measurements so far
and the phases. In short, on 150 bp reads, k-mers enumerated from the HMMs plus an ungapped score along
the seed's diagonal keep about 98% of `hmmsearch`'s hits at about 13 candidate models per read, out
of 30,134. Wall-clock speed-up is not measured yet.

Training uses the HMMs and Pfam's own alignments only. MGnify proteins are for evaluation.

## Tools

`hmmsearch`, MMseqs2 and nail run from one image:

```bash
docker build -t kfpfam-tools docker/
```

Local data goes in `data/` (not tracked). Pfam-A comes from
[the Pfam FTP site](https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/).

Related: [kmer-functional-profiler](https://github.com/timrozday-mgnify/kmer-functional-profiler),
[kfp-explore](https://github.com/timrozday-mgnify/kfp-explore).

Licence: GPL-3.0-or-later.
