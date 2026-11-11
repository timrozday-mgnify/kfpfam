"""Score recruit.py output on the simulated sample (measure/recruit_sim.py) against the reads' origin.

Truth is per gene segment of a read: the Pfam domains at GA overlapping it by >= 10 aa ('ga'), plus
domains below GA ('other'). A locus is judged against the segments its alignment overlaps on the read
(by at least half the alignment), or against every segment of the read when it has no read position.
A kept locus is right if its top is a true target, 'in set' if a set member is, 'clan' if only a clan
sibling is, false otherwise. Recall is per segment with a GA domain: found when a kept locus over it is
right or holds a true target.
Usage: recruit_eval.py TRUTH.tsv OUT_DIR MODELS.tsv
"""

import csv
import json
import sys
from collections import Counter, defaultdict


def main():
    truth_path, out, models = sys.argv[1:]
    clan = {r["acc"]: r["clan"] or r["acc"] for r in csv.DictReader(open(models), delimiter="\t")}
    segs = defaultdict(list)
    for r in csv.DictReader(open(truth_path), delimiter="\t"):
        r["ok"] = set(filter(None, (r["ga"] + "," + r["other"]).split(",")))
        segs[r["read"]].append(r)
    S = json.load(open(f"{out}/summary.json"))
    targets = {r["acc"]: r for r in csv.DictReader(open(f"{out}/targets.tsv"), delimiter="\t")}
    loci = list(csv.DictReader(open(f"{out}/loci.tsv"), delimiter="\t"))

    def over(locus):
        f, t = int(locus["read_from"]), int(locus["read_to"])
        if f < 0:
            return segs[locus["read"]]
        return [s for s in segs[locus["read"]] if min(t, int(s["seg_to"])) - max(f, int(s["seg_from"])) >= 0.5 * (t - f)]

    found = set()  # (read, segment index) found right
    verdicts = defaultdict(Counter)
    for lc in loci:
        if lc["status"] == "none":
            continue
        ss = over(lc)
        ok = set().union(*(s["ok"] for s in ss)) if ss else set()
        st = {x.split(":")[0] for x in lc["set"].split(",")}
        if lc["top"] in ok:
            v = "right"
        elif st & ok:
            v = "in set"
        elif {clan[a] for a in st} & {clan[a] for a in ok}:
            v = "clan"
        else:
            v = "false"
        verdicts[lc["status"]][v] += 1
        if v in ("right", "in set"):
            for s in ss:
                if s["ok"] & st:
                    found.add((lc["read"], id(s)))
    print(f"reads {len(segs)}, loci {len(loci)}")
    for status in ("ga", "one", "set"):
        c = verdicts[status]
        n = sum(c.values())
        print(f"  {status:4} {n:6}: " + "  ".join(f"{k} {c[k] / max(n, 1):.3f}" for k in ("right", "in set", "clan", "false")))
    rc = verdicts["one"] + verdicts["set"]
    n = sum(rc.values())
    print(f"recruited loci {n}: estimated FDR {S['null_reads_kept_scaled'] / max(S['reads_recruited'], 1):.3f}, "
          f"true FDR {rc['false'] / max(n, 1):.3f}")

    units = [(k, s) for k, ss in segs.items() for s in ss if s["ga"]]
    print(f"\ngene segments with a GA domain: {len(units)}, found {sum((k, id(s)) in found for k, s in units) / len(units):.3f}")
    for name, f, bins in (("identity", lambda s: float(s["identity"]), (1.0, 0.9, 0.8, 0.7, 0.6, 0.5)),
                          ("coverage", lambda s: float(s["coverage"]), ((0, 1), (1, 3), (3, 10), (10, 31)))):
        print(f" by {name}: " + "  ".join(
            f"{b}: {sum((k, id(s)) in found for k, s in us) / max(len(us), 1):.3f} ({len(us)})"
            for b in bins for us in [[(k, s) for k, s in units
                                      if (f(s) == b if name == "identity" else b[0] <= f(s) < b[1])]]))
    multi = [k for k, ss in segs.items() if sum(bool(s["ga"]) for s in ss) > 1]
    both = sum(all((k, id(s)) in found for s in segs[k] if s["ga"]) for k in multi)
    multi_any = [k for k, ss in segs.items() if sum(bool(s["ok"]) for s in ss) > 1]
    both_any = sum(sum((k, id(s)) in found for s in segs[k] if s["ok"]) > 1 for k in multi_any)
    print(f"reads with GA domains in two gene segments: {len(multi)}, both found {both}; "
          f"with any domain in two segments: {len(multi_any)}, two found {both_any}")
    for kind in ("ga", "weak"):
        ts = [a for a, r in targets.items() if r["confirmed"] == kind]
        ok_any = {a for ss in segs.values() for s in ss for a in s["ok"]}
        print(f"confirmed by {kind}: {len(ts)}, of which a true domain of the sample {sum(a in ok_any for a in ts) / max(len(ts), 1):.3f}")


if __name__ == "__main__":
    main()
