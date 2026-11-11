"""Reads with a hit, by first call (GA) and post-process recruitment; and whether the alignment position
separates real weak pairs from reversed ones.

Position classes of a pair's strongest domain (HMM columns hmm_from..hmm_to of M, ORF ali_from..ali_to of L):
  'cut'   - the alignment runs to an ORF end where the HMM continues (ali_from <= 3 and hmm_from > 3, or
            ali_to >= L - 2 and hmm_to < M - 2): the read stops inside the domain, which explains a short hit;
  'inner' - the ORF continues past both alignment ends but the HMM does too: the read spans more and still
            aligns only partly;
  'whole' - the HMM is covered end to end (within 3 columns).
Usage: recruit_reads.py REAL.domtbl DECOY.domtbl MODELS.tsv MODE_A.tsv BAR NREADS
"""

import csv
import sys
from collections import Counter, defaultdict

Z = 20_000_000
E_TOPS = (1e-1, 1, 10, 100, 1000)


def pos_class(L, M, hf, ht, af, at):
    if hf <= 3 and ht >= M - 2:
        return "whole"
    if (af <= 3 and hf > 3) or (at >= L - 2 and ht < M - 2):
        return "cut"
    return "inner"


def pairs(path, ga):
    best = {}
    for line in open(path):
        if line[0] != "#":
            c = line.split()
            a = c[4].split(".")[0]
            if a in ga:
                k = (int(c[0][1:].split("_")[0]), a)
                s = float(c[13])
                if k not in best or s > best[k][0]:
                    best[k] = (s, float(c[6]), pos_class(int(c[2]), int(c[5]), *map(int, c[15:17]), *map(int, c[17:19])))
    return best


def main():
    real, decoy, models, mode_a, bar, nreads = sys.argv[1:]
    rows = list(csv.DictReader(open(models), delimiter="\t"))
    ga = {r["acc"]: float(r["ga"]) for r in rows}
    accs = [r["acc"] for r in rows]
    R, D = pairs(real, ga), pairs(decoy, ga)
    A = {(int(q[1:].split("_")[0]), accs[int(m)]) for q, m, _, g in map(str.split, open(mode_a)) if float(g) >= float(bar)}

    for label, keep in (("hmmsearch on all ORFs", lambda k: True), ("Mode A candidates, then HMMER", lambda k: k in A)):
        Rk = {k: v for k, v in R.items() if keep(k)}
        conf = {a for (q, a), v in Rk.items() if v[0] >= ga[a]}
        first = {q for (q, a), v in Rk.items() if v[0] >= ga[a]}
        weak = defaultdict(list)
        for (q, a), v in Rk.items():
            if q not in first and a in conf:
                weak[q].append(v)
        other = {q for (q, a) in Rk if q not in first and q not in weak}
        print(f"\n{label}: {nreads} reads; with any pair at E <= 1000: {len({q for q, _ in Rk})}")
        print(f"  first hit (GA call): {len(first)}; confirmed targets {len(conf)}")
        for e in E_TOPS:
            n = sum(any(v[1] <= e for v in vs) for vs in weak.values())
            nd = len({q for (q, a), v in D.items() if a in conf and v[0] < ga[a] and v[1] <= e})
            print(f"  + recruited at best weak E <= {e:g}: {n} reads (reversed {nd}, FDR {nd / max(n, 1):.3f}); total {len(first) + n}")
        print(f"  weak pairs only to unconfirmed targets: {len(other)}")

    # position: real weak pairs to confirmed targets (reads with no GA call) against reversed, by E bin
    conf = {a for (q, a), v in R.items() if v[0] >= ga[a]}
    first = {q for (q, a), v in R.items() if v[0] >= ga[a]}
    print("\nposition of weak pairs to confirmed targets, real / reversed (share of real in brackets), by E bin")
    lo = 0
    for e in E_TOPS:
        r = Counter(v[2] for (q, a), v in R.items() if a in conf and q not in first and lo < v[1] <= e)
        d = Counter(v[2] for (q, a), v in D.items() if a in conf and v[0] < ga[a] and lo < v[1] <= e)
        tot = sum(r.values())
        print(f"  E {lo:g}-{e:g}: " + "  ".join(f"{c} {r[c]}/{d[c]} ({r[c] / max(tot, 1):.2f}) FDR {d[c] / max(r[c], 1):.2f}" for c in ("cut", "inner", "whole")))
        lo = e
    r = Counter(v[2] for (q, a), v in R.items() if v[0] >= ga[a])
    print("  GA pairs for reference: " + "  ".join(f"{c} {r[c]}" for c in ("cut", "inner", "whole")))


if __name__ == "__main__":
    main()
