"""How much of Pfam, and of real metagenome output, the local benchmark actually tests.

The local recall figure (measure/index_score.py, 50,000 ERR10889147 reads) is a fraction of the
(read, model) truth pairs at or above GA that the filter emits as candidates. That denominator is
built from one human gut sample, so it only exercises the families that sample happens to contain.
This reports the size of that blind spot, and checks the invariant the headline recall cannot see:
that every model in the index can produce a candidate at all.

  representativeness.py coverage PREP TRUTH.domtbl [FAMILY_READS.tsv]
      TRUTH.domtbl      hmmsearch at E<=1000 over all models on the benchmark's reads
      FAMILY_READS.tsv  optional `accession<TAB>reads` from a real run, to weight the blind spot
                        by the reads it carries (e.g. summed *_pfam.txt.gz from a RawR run)

  representativeness.py models PREP SEARCH.tsv [SEARCH.tsv ...]
      SEARCH.tsv        `filter join` output (query, model index, diagonal, gapped). Reports the
                        models that never appear. A model absent from every search is a build
                        defect, not a sensitivity setting: it can never be seeded at any bar.
"""

import csv
import sys
from pathlib import Path


def models(prep: Path) -> list[dict]:
    with open(prep / "models.tsv") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def truth_pairs(path: Path, ga: dict[str, float]) -> tuple[set, set]:
    """(read, acc) pairs whose best domain reaches GA, and the families among them."""
    best: dict[tuple[str, str], float] = {}
    with open(path) as f:
        for line in f:
            if line[0] == "#":
                continue
            c = line.split()
            k = (c[0].split("_")[0], c[4].split(".")[0])
            best[k] = max(best.get(k, float("-inf")), float(c[13]))
    pairs = {k for k, s in best.items() if k[1] in ga and s >= ga[k[1]]}
    return pairs, {k[1] for k in pairs}


def coverage(prep: Path, truth: Path, family_reads: Path | None) -> None:
    rows = models(prep)
    ga = {r["acc"]: float(r["ga"]) for r in rows}
    pairs, fams = truth_pairs(truth, ga)
    print(f"local truth: {len(pairs):,} (read, model) pairs at or above GA, "
          f"over {len(fams):,} families")
    print(f"  that is {len(fams) / len(rows):.1%} of the {len(rows):,} models in the index; "
          f"{len(rows) - len(fams):,} models are never tested at GA")
    # one family's worth of pairs is the smallest loss the headline recall can resolve
    print(f"  one pair is {1 / len(pairs):.3%} of recall, so a defect costing fewer than about "
          f"{int(0.005 * len(pairs))} pairs reads as 0.995 or better")

    if family_reads is None:
        return
    reads: dict[str, float] = {}
    with open(family_reads) as f:
        for line in f:
            a, n = line.split("\t")[:2]
            reads[a.split(".")[0]] = reads.get(a.split(".")[0], 0) + float(n)
    total = sum(reads.values())
    inside = sum(n for a, n in reads.items() if a in fams)
    untested = {a: n for a, n in reads.items() if a not in fams}
    print(f"\nagainst the real run: {len(reads):,} families, {total:,.0f} annotated reads")
    print(f"  in families the local truth tests at GA: {inside / total:.1%} of reads")
    print(f"  in families it does not: {1 - inside / total:.1%} of reads, "
          f"over {len(untested):,} families")
    top = sorted(untested.items(), key=lambda kv: -kv[1])[:10]
    print("  largest untested families by real reads:")
    for a, n in top:
        print(f"    {a:10s} {n:12,.0f}")


def model_coverage(prep: Path, searches: list[Path]) -> None:
    rows = models(prep)
    seen: set[int] = set()
    for p in searches:
        with open(p) as f:
            for line in f:
                seen.add(int(line.split()[1]))
        print(f"{p}: {len(seen):,} distinct models seen so far")
    never = [i for i in range(len(rows)) if i not in seen]
    print(f"\n{len(rows) - len(never):,} of {len(rows):,} models produced a candidate; "
          f"{len(never):,} never did")
    if never:
        print("  these cannot be seeded at any bar -- check prep and `filter build`:")
        for i in never[:25]:
            print(f"    index {i:6d}  {rows[i]['acc']:10s} {rows[i]['name']:22s} "
                  f"length {rows[i]['length']:>5s} GA {rows[i]['ga']:>6s}")
        if len(never) > 25:
            print(f"    ... and {len(never) - 25} more")


def _selfcheck() -> None:
    assert truth_pairs.__doc__
    assert Path(__file__).exists()


def main() -> None:
    _selfcheck()
    if len(sys.argv) < 4:
        sys.exit(__doc__)
    mode, prep = sys.argv[1], Path(sys.argv[2])
    if mode == "coverage":
        coverage(prep, Path(sys.argv[3]),
                 Path(sys.argv[4]) if len(sys.argv) > 4 else None)
    elif mode == "models":
        model_coverage(prep, [Path(p) for p in sys.argv[3:]])
    else:
        sys.exit(__doc__)


if __name__ == "__main__":
    main()
