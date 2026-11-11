# Recruiting weak reads to Pfam targets

`recruit/` assigns reads with HMMER hits too weak to call at GA to Pfam families the sample confirms. It
keeps a set of targets where a read is ambiguous, and writes each kept hit's alignment. Its output is the
input to the Pfam-anchored search of MGnify proteins, which keys k-mers by HMM column and lives in
kmer-functional-profiler's `docs/pfam-anchor-plan.md`. The exploration and measurements behind each
choice are in [the recruitment plan](recruit-plan.md).

Built 2026-11-07, on one human gut metagenome (ERR10889147) and a simulated sample. Nothing has run on
the HPC or on other biomes.

## Where it sits

```text
reads -> six-frame ORFs + coordinates            measure/real_orfs.py
      -> k-mer filter, Mode A                    filter join (k4 T6.5, diagonal 20, gate 15.5)
      -> HMMER on candidates, every pair         recruit/pairs.py (gapped >= 22.25, Z = 20 M, E <= 1000)
         (the same for the reversed ORFs)        recruit/recruit.py reverse
      -> recruitment                             recruit/recruit.py run
      -> hits.tsv: query, family, alignment      -> Pfam-anchored MGnify search
```

Commands are in the README's "Recruiting weak reads" section. `python3 recruit/test_recruit.py` runs the
tests.

## Method

1. **Every HMMER pair.** `pairs.py` runs pyhmmer on the filter's candidate (ORF, model) pairs. For each
   pair it writes the strongest domain: score, E-values, HMM, ORF and envelope coordinates, and a CIGAR.
   On 50,000 reads, its scores and coordinates are identical to `hmmsearch` on all ORFs for every shared
   pair, and it finds .997 of `hmmsearch`'s GA pairs.
2. **The null.** The same steps run on the reversed ORFs. HMMER's E-values are about 7× too optimistic on
   short reversed ORFs, so the null is measured per sample, not taken from E.
3. **Confirmed targets.** These are families with at least one read at GA (`--min-ga`). Weak hits only
   recruit to these. That cuts chance hits about 15-fold: at E ≤ 10, the FDR is .02 against these
   families and .37 against the rest.
4. **Loci.** With ORF coordinates (`--orf-coords`), each alignment is placed on the forward read.
   - Alignments sharing at least half of the shorter one compete as one locus, whatever their frame.
   - Separate stretches of a read (two genes, or two domains) are separate loci.
   - On the real sample, 346 of 350 reads with GA calls on two ORFs had both over the same bases in
     different frames, so one frame was wrong. Only 4 were separate genes.
5. **P(true): the chance a locus comes from any confirmed target.**
   - A weighted logistic fit of real against reversed pairs to confirmed targets. Its inputs are log E
     with hinges, the alignment's position, score against expected score, and ORF length.
     - Position is 'cut' (the read ends inside the domain, which explains a short hit), 'inner' or
       'whole'.
     - Expected score is what a member would reach over the same HMM columns.
   - The fitted odds convert to a likelihood ratio, LR = π (odds − 1) / (1 − π), where π is the null's
     share of pairs.
   - The locus's best pair gives P(true), with the null's share of reads as the prior.
   - Against E alone, the four inputs recruit 11–12% more reads at a 1–2% FDR target. Position adds about
     4%.
6. **FDR cut, per read.** Reversed reads are scored by the same model. The cut is the lowest P(true) at
   which scaled reversed reads kept ÷ real reads kept ≤ `--fdr` (default 0.02). Every locus has a
   q-value, and loci below the cut stay in the output with status `none`.
7. **Which target.** P(t | true) over the locus's targets is proportional to prior(t) × LR.
   - prior(t) is the family's GA reads plus its recruited mass, from one EM round (`--em-iters`).
   - A locus goes to its top target if P ≥ `--tau` (0.8). Otherwise it keeps the smallest set summing to
     TAU, across clans if need be.
   - GA loci get a set the same way, so a GA call with a near tie in another frame or family keeps both.
8. **Targets from weak evidence.** For a family with no GA read, its weak reads' log LR is summed, counting
   the best read in each 15-column window of the HMM only once. Repeats in deep proteins otherwise pile
   reversed reads onto one family. Families whose sum beats the reversed reads' at `--target-fdr` (0.05)
   join the confirmed set, and steps 5–7 run again. This needs the whole null.
9. **Depth check.** Recruited mass per family is compared with rho × GA reads under a negative binomial,
   flagging families with far more or far fewer weak reads than expected. rho was 0.50–0.58 recruited
   reads per GA read.

## Outputs

All files are TSV unless named otherwise.

### `hits.tsv`: for the anchored search

One row per (locus, target in its set).

| Column | Meaning |
| --- | --- |
| `query` | ORF name, `R<read>_<n>`; the sequence is in the ORF FASTA |
| `read`, `locus` | the read, and which stretch of it (0, 1, …) |
| `status` | `ga` (a pair at GA), `one`, `set` (ambiguous), `none` (below the FDR cut) |
| `top` | 1 for the locus's best target: keep `top = 1` for one call per locus |
| `weight` | P(t \| true) × P(true); P(t \| true) alone for `ga` |
| `p_true`, `q` | the locus's probability of coming from any confirmed target, and its q-value (0 for `ga`) |
| `pfam_accession` | integer, as pfam-reannotate stores it (`PF00004` → 4) |
| `seq_evalue`, `score` | sequence E-value at Z = 20 M; best domain's bit score |
| `c_evalue`, `i_evalue` | domain conditional and independent E-values |
| `hmm_from`, `hmm_to` | HMM match columns, 1-based, inclusive |
| `env_from`, `env_to`, `ali_from`, `ali_to` | envelope and alignment on the ORF, 1-based, inclusive |
| `cigar` | over the alignment: M match, I residue in the ORF only, D model column only, as `reannotate.py` writes |
| `above_ga` | domain score ≥ the family's GA (pfam-reannotate also asks the sequence score) |
| `read_from`, `read_to` | the alignment on the forward read, bases, 0-based, half-open; −1 without coordinates |
| `orf_len` | ORF length, aa |

With `ali_from`, the CIGAR and the ORF sequence, the aligned residues can be rebuilt. With `hmm_from`, the
CIGAR gives the HMM column of each residue, which is what column-keyed k-mers need.

This was checked with the anchored plan's own mapping, `column_residues` in kmer-functional-profiler's
`workflows/pfam-anchor/columns.py`, on the 50,000 reads' 16,362 hits:
- Every hit maps consistently: the match columns equal the CIGAR's M count and end within `hmm_to` and
  `ali_to`.
- The ORF residue equals the HMM's consensus at its column for .48 of GA hits' residues, and .38 of
  recruited ones.

### The other outputs

- **`loci.tsv`:** one row per locus: status, P(true), q, top, the set with its probabilities and clans,
  and the top pair's coordinates, position class and score ratio.
- **`targets.tsv`:** per family: confirmed by `ga` or `weak`, GA reads, recruited reads and mass,
  expected recruited reads and the depth flag.
- **`candidates.faa`:** each kept ORF once per target in its set, with `target`, `status`, `weight` and
  `top` in the header.
- **`summary.json`:** the fit's coefficients, the cuts, and counts.
- **`pairs.tsv`** (from `pairs.py`): every pair at E ≤ 1000, with the same alignment columns. This is
  the full set, for a recall-first use that doesn't want recruitment's cut to confirmed families.

## Measured

| | 50,000 reads | 1.07 M reads | Simulation |
| --- | --- | --- | --- |
| Reads (ORFs) | 50,000 (618 k) | 1,070,914 (6.63 M) | 95,922 (224 k) |
| Families at GA, from weak evidence | 1,968, 297 | 3,695, 345 | 515, 145 |
| Loci at GA | 8,064 | 172,619 | 5,470 |
| Recruited loci: one target, set | 5,198, 158 | 104,571, 4,179 | 6,755, 178 |
| Reversed reads kept (estimated FDR) | 106 (.020) | 2,168 (.020) | 137 (.020) |
| `recruit.py run` | 3 s | 26 s, 3.0 GB | 4 s |

The 1.07 M-read null costs another filter run (124 s wall on 6 threads, 3.4 GB) and pyhmmer (107 s).
Filter and pyhmmer on the real reads cost the same again.

**Simulation, known origin.** 760 contigs of 2–6 MGnify proteins, 1,144 of the 3,000 proteins with a GA
domain. Each contig has an identity of 100–50% (BLOSUM62 substitutions) and a coverage of 0.3–30. Genes
lie on random strands with 0–150 bp spacers, and reads are 100 or 150 bp from either strand.

| | right | right target in the set | clan only | false |
| --- | --- | --- | --- | --- |
| GA loci | .999 | .001 | 0 | .001 |
| Recruited, one target | .960 | | .003 | .037 |
| Recruited, set | .916 | .028 | .011 | .045 |

Gene segments with a GA domain found, by identity: 1.0 .55, .9 .43, .8 .29, .7 .26, .6 .15, .5 .06. An
earlier one-protein-per-read simulation found that recruitment adds more than GA calls give from 80%
identity down, and more than doubles what is found at 70% and below.

**Null size.** Per-read results need about 10,000 null reads. At 1.07 M reads, 2% of the null gives
recruited counts within 1.5% of the whole null; at 50,000 reads, 20% is within 3%.

## Limits

- **The reversed-read FDR is a lower bound.** On the simulation the true FDR was .037 at an estimated
  .020. In the earlier simulation it was .056 at .019 and .016 at .001, so the true FDR runs about 2–4
  points above the estimate. Two kinds of hit escape the null:
  - ORFs in other frames of real coding DNA, which gave 1.3% of false GA calls in the first simulation
    (0.1% in the contig one, where frames compete within a locus);
  - real proteins hitting a family they aren't annotated with at E ≤ 1.

  These false hits come in runs from one region. The alignment step against MGnify is what removes
  them. A coding-frame score would remove the first kind (see the plan).
- **Recruitment is a presence cut by family.** Only families with a GA read, or summed weak evidence,
  receive weak reads. A recall-first use should take `pairs.tsv` and use `hits.tsv` to rank.
- **Reach.** Only reads with an HMMER pair at E ≤ 1000 are here. Reads under about 15 aa of domain
  rarely hit: windows of 10, 15, 20 and 25 aa cut from GA reads got a pair .04, .29, .62 and .85 of the
  time. ORFs are at least 15 aa.
- **One alignment per pair.** The ORF's strongest domain for each family. A second domain of the same
  family on the same ORF is not written.
- **Thresholds and versions.** Pfam 38.2 (Pfam-A, the same file as kmer-functional-profiler). Z = 20 M
  for every E-value. The filter's Mode A sits in front, losing .003 of GA pairs and .03–.11 of the
  weak bins against `hmmsearch` on all ORFs.
- **Tested on one sample.** Other biomes, and checking the FDR target on them, are still to do.

## Files

- `recruit/pairs.py`: the HMMER step.
- `recruit/recruit.py`: `reverse` and `run`.
- `recruit/test_recruit.py`: the tests.
- `measure/real_orfs.py`: ORFs and `coords.tsv`.
- `measure/recruit_sim.py`, `measure/recruit_eval.py`: the simulation and its scoring.
- `measure/recruit_{explore,trunc,model,reads}.py`: the exploration.
- `docs/recruit-plan.md`: the plan, its measurements and the build log.
