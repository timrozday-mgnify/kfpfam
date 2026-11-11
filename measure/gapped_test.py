"""Does a gapped score beat the ungapped diagonal at the same number of candidates?

Rescores what a loose filter passes (checks/gapped_export.mjs) with `scan gapped`, then compares
recall by match strength at equal candidates per read.
Usage: gapped_test.py WORK SCAN_DIR PREP_DIR
"""

import csv
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

LOWS = (30, 10, 0, -5, -10, -20, float("-inf"))
NAMES = ("GA+30", "+10-30", "GA-+10", "-5-GA", "-10--5", "-20--10", "<-20")
GAPS = ("ungapped", "gap 8/2", "gap 6/1", "gap 4/1", "gap 3/.5")
TARGETS = (1.3, 3.5, 11.6, 40)
TRUTH_READS, SCANNED = 50000, 300
BAR = 18.0  # the export's bar: background pairs below it were not kept, so real ones must not count


def rows(path):
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def main():
    work, scan, prep = (Path(a) for a in sys.argv[1:4])
    seqs = [line.split("\t")[0] for line in open(scan / "pairs.seqs.tsv")]
    pairs, bg, decoy = rows(work / "pairs.tsv"), rows(work / "bg.tsv"), rows(work / "decoy.tsv")
    orfs = {}
    for name in ("bg", "decoy"):
        orfs[name] = defaultdict(list)
        for line in open(scan / "real" / f"{name}.orfs.tsv"):
            read, seq = line.rstrip("\n").split("\t")
            orfs[name][read].append(seq)
    jobs = work / "jobs.tsv"
    with open(jobs, "w") as o:
        for p in pairs:
            if float(p["diag"]) > 0:
                o.write(f"P{p['g']}\t{p['model']}\t{seqs[int(p['g'])]}\n")
        for name, table in (("bg", bg), ("decoy", decoy)):
            for r in table:
                for s in orfs[name][r["read"]]:
                    o.write(f"{name}:{r['read']}:{r['model']}\t{r['model']}\t{s}\n")
    if not (work / "scores.tsv").exists():
        subprocess.run(["measure/scan/target/release/scan", "gapped", str(prep), str(jobs), str(work / "scores.tsv")], check=True)
    best = {}
    for line in open(work / "scores.tsv"):
        c = line.rstrip("\n").split("\t")
        v = list(map(float, c[2:]))
        best[c[0]] = [max(a, b) for a, b in zip(best[c[0]], v)] if c[0] in best else v

    def binned(score, ga):
        return next(i for i, low in enumerate(LOWS) if score - ga >= low)

    real = [(p, best.get(f"P{p['g']}"), binned(float(p["score"]), float(p["ga"]))) for p in pairs if p["set"] == "real"]
    fake = [(p, best.get(f"P{p['g']}"), binned(float(p["score"]), float(p["ga"]))) for p in pairs if p["set"] == "realdecoy"]
    unrated = {n: [(best[f"{n}:{r['read']}:{r['model']}"], float(r["diag"])) for r in t if float(r["score"]) <= -900] for n, t in (("bg", bg), ("decoy", decoy))}
    total = [sum(float(p["weight"]) for p, _, b in real if b == i) for i in range(7)]

    def at(metric, t, funnel=BAR):
        """Rates at a score threshold, among pairs whose diagonal reached the funnel's bar."""
        ok = [0.0] * 7
        for p, s, b in real:
            if s and s[metric] >= t and float(p["diag"]) >= funnel:
                ok[b] += float(p["weight"])
        cand = sum(ok) / TRUTH_READS + sum(s[metric] >= t and d >= funnel for s, d in unrated["bg"]) / SCANNED
        fakes = sum(float(p["weight"]) for p, s, _ in fake if s and s[metric] >= t and float(p["diag"]) >= funnel) / TRUTH_READS
        fakes += sum(s[metric] >= t and d >= funnel for s, d in unrated["decoy"]) / SCANNED
        above = sum(ok[:3]) / max(sum(total[:3]), 1)
        return cand, fakes, above, [o / max(x, 1) for o, x in zip(ok, total)]

    print("hits per bin", dict(zip(NAMES, map(round, total))))
    for want in TARGETS:
        print(f"\nnear {want} candidates per read: threshold cand decoy | >=GA | {' '.join(NAMES)}")
        for metric, name in enumerate(GAPS):
            t = 14.0
            while at(metric, t)[0] > want:
                t += 0.25
            c, f, above, by = at(metric, t)
            print(f"{name:9} t={t:5.2f} {c:6.2f} {f:6.2f} | {above:.3f} | " + " ".join(f"{x:.2f}" for x in by))
    print("\ngap 6/1 rescoring only pairs whose diagonal reached a bar: bar aligned-per-read t cand decoy | >=GA | bins")
    for funnel in (18, 20, 22, 24):
        aligned = at(0, 0.0, funnel)[0]
        for want in (1.3, 3.5):
            t = 14.0
            while at(2, t, funnel)[0] > want:
                t += 0.25
            c, f, above, by = at(2, t, funnel)
            print(f"{funnel} {aligned:6.1f} t={t:5.2f} {c:5.2f} {f:5.2f} | {above:.3f} | " + " ".join(f"{x:.2f}" for x in by))


if __name__ == "__main__":
    main()
