"""Compare amino-acid sets: this project's six-frame ORFs against the pipeline's.

Usage: compare_orfs.py PROJECT.faa PIPELINE_RAW.faa PIPELINE.faa FRAGMENTS.fasta
"""

import gzip
import json
import sys
from collections import Counter
from pathlib import Path

import numpy as np


def fasta(path):
    name, seq = None, []
    with open(path) as f:
        for line in f:
            if line[0] == ">":
                if name is not None:
                    yield name, "".join(seq)
                name, seq = line[1:].split()[0], []
            else:
                seq.append(line.strip())
    if name is not None:
        yield name, "".join(seq)


def stats(label, c):
    keys, counts = zip(*c.items())
    lengths = np.repeat([len(k) for k in keys], counts)
    print(f"{label}: {lengths.size:,} sequences, {lengths.sum():,} residues, mean {lengths.mean():.1f}, "
          f"p10/50/90/99/max {np.percentile(lengths, [10, 50, 90, 99]).astype(int).tolist()} {lengths.max()}, "
          f">=50 aa {np.mean(lengths >= 50):.3f}, residues in >=50 aa {lengths[lengths >= 50].sum() / lengths.sum():.3f}, "
          f"with X {sum(v for k, v in c.items() if 'X' in k):,}")


def diff(label, a, b):
    shared = sum((a & b).values())
    only_a, only_b = a - b, b - a
    mean = lambda c: sum(len(k) * v for k, v in c.items()) / max(sum(c.values()), 1)  # noqa: E731
    print(f"{label}: shared {shared:,}, only first {sum(only_a.values()):,} (mean {mean(only_a):.1f} aa), "
          f"only second {sum(only_b.values()):,} (mean {mean(only_b):.1f} aa)")
    return only_a, only_b


def main():
    project, raw, pipeline, fragments = sys.argv[1:5]
    p = Counter(s for _, s in fasta(project))
    r = Counter(s for _, s in fasta(raw))
    m = Counter(s for _, s in fasta(pipeline))
    stats("project", p)
    stats("translator on untouched reads", r)
    stats("pipeline", m)
    a, b = diff("translator vs project", r, p)
    print("  only translator:", list(a)[:3], "only project:", list(b)[:3])
    diff("pipeline vs project", m, p)
    # reformat.sh drops read descriptions, so fastp's "merged_" marks are read from its FASTQ
    with gzip.open(Path(fragments).parent / "merged.fastq.gz", "rt") as f:
        rows = [(" merged_" in h, len(q) - 1) for h, q, _, _ in zip(f, f, f, f)]
    merged = np.array([x for x, _ in rows])
    lengths = np.array([x for _, x in rows])
    print(f"fragments {lengths.size:,}, nt mean {lengths.mean():.1f}; merged {merged.sum():,} nt mean "
          f"{lengths[merged].mean():.1f} p10/50/90 {np.percentile(lengths[merged], [10, 50, 90]).astype(int).tolist()}; "
          f"unmerged {(~merged).sum():,} nt mean {lengths[~merged].mean():.1f}")
    print(f"pipeline sequences per fragment {sum(m.values()) / lengths.size:.2f}")
    for name in ("qc.json", "merge.json"):
        j = json.loads((Path(fragments).parent / name).read_text())
        s = j["summary"]
        print(name, s["before_filtering"]["total_reads"], s["before_filtering"]["total_bases"], s["after_filtering"]["total_reads"], s["after_filtering"]["total_bases"], j["filtering_result"], "adapter-trimmed", j.get("adapter_cutting", {}).get("adapter_trimmed_reads"), "insert peak", j.get("insert_size", {}).get("peak"), "dup", j.get("duplication", {}).get("rate"))


if __name__ == "__main__":
    main()
