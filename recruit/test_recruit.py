"""Checks of recruit.py's logic on hand-made pairs. Run: python3 recruit/test_recruit.py (or pytest)."""

import os
import random
import tempfile

import numpy as np

import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "measure"))
from real_orfs import CODON, RC, orfs_at  # noqa: E402

from recruit import Pair, annotate, choose, fdr_cut, fit_lr, loci, place, reverse, score_loci  # noqa: E402


def pair(orf, acc, score=20.0, e=1.0, hf=1, ht=30, af=1, at=30, L=30):
    return Pair(orf, orf.rsplit("_", 1)[0], acc, L, score, e, hf, ht, af, at)


def test_choose():
    assert choose({"a": 0.9, "b": 0.1}, 0.8) == ("one", ["a"])
    assert choose({"a": 0.5, "b": 0.4, "c": 0.1}, 0.8) == ("set", ["a", "b"])


def test_loci():
    # two domains apart on one ORF split; overlapping hits and other frames join the best locus
    ps = [pair("R1_0", "A", 30, af=1, at=20), pair("R1_0", "B", 25, af=25, at=45),
          pair("R1_0", "C", 10, af=3, at=18), pair("R1_2", "D", 12)]
    lc = loci(ps)
    assert sorted(sorted(p.acc for p in g) for g in lc.values()) == [["A", "C", "D"], ["B"]]


def test_position():
    models = {"A": {"M": 100}}
    cut, inner, whole = pair("R1_0", "A", hf=50, ht=79, af=1, at=30), pair("R2_0", "A", hf=40, ht=60, af=5, at=25), \
        pair("R3_0", "A", hf=1, ht=100, af=2, at=29)
    annotate([cut, inner, whole], models, {"A": np.arange(101.0)})
    assert (cut.pos, inner.pos, whole.pos) == ("cut", "inner", "whole")
    assert abs(cut.ratio - 20 / 30) < 1e-9


def test_fdr_cut():
    real = [0.99] * 90 + [0.5] * 10
    null = [0.99] * 1 + [0.5] * 9
    cut, q = fdr_cut(real, null, 1.0, 0.02)
    assert cut == 0.99 and q[0.99] <= 0.02 and q[0.5] == 0.1


def test_lr_and_sets():
    # strong pairs are real, weak ones half null: a strong pair must outrank a weak one at equal prior
    rng = random.Random(0)
    real = [pair(f"R{i}_0", "A", e=10 ** rng.uniform(-6, -1)) for i in range(300)] + \
           [pair(f"R{i}_0", "A", e=10 ** rng.uniform(1, 3)) for i in range(300, 600)]
    null = [pair(f"D{i}_0", "A", e=10 ** rng.uniform(1, 3)) for i in range(250)]
    for p in real + null:
        p.pos, p.ratio = "cut", 0.5
    fit_lr(real, null, 1.0, ["e"], real + null)
    assert real[0].lr > real[-1].lr
    S = score_loci(loci([real[0]]), 0.3, {"A": 1})
    assert S[("R0", 0)][0] > 0.9


def translate(dna):
    return "".join(CODON[dna[i:i + 3]] for i in range(0, len(dna) - 2, 3))


def test_place():
    # a read: gene A forward, a stop, then gene B on the reverse strand; each ORF's alignment maps back
    # to the bases it was translated from, and the two genes become two loci
    rng = random.Random(3)
    sense = [c for c, a in CODON.items() if a != "*"]
    a = "".join(rng.choice(sense) for _ in range(20)) + "TAA"
    b = "".join(rng.choice(sense) for _ in range(20))
    read = "G" + a + b[::-1].translate(RC)  # B is read on the minus strand
    d = tempfile.mkdtemp()
    found = orfs_at(read)
    with open(os.path.join(d, "c.tsv"), "w") as f:
        f.write("orf\tstrand\tstart\tread_len\n")
        for j, (prot, sign, start) in enumerate(found):
            f.write(f"R1_{j}\t{sign}\t{start}\t{len(read)}\n")
    ia = next(j for j, (prot, sign, _) in enumerate(found) if sign == "+" and translate(a[:-3]) in prot)
    ib = next(j for j, (prot, sign, _) in enumerate(found) if sign == "-" and translate(b) in prot)
    pa = pair(f"R1_{ia}", "A", 30, af=2, at=11, L=len(found[ia][0]))
    pb = pair(f"R1_{ib}", "B", 20, af=1, at=10, L=len(found[ib][0]))
    place([pa, pb], os.path.join(d, "c.tsv"))
    for p, (prot, sign, _) in ((pa, found[ia]), (pb, found[ib])):
        dna = read[p.rf:p.rt]
        got = translate(dna if sign == "+" else dna[::-1].translate(RC))
        assert got == prot[p.af - 1:p.at], (got, prot[p.af - 1:p.at])
    assert len(loci([pa, pb])) == 2
    # without coordinates the other ORF joins the best locus
    assert len(loci([pair(f"R1_{ia}", "A", 30), pair(f"R1_{ib}", "B", 20)])) == 1


def test_reverse():
    d = tempfile.mkdtemp()
    with open(os.path.join(d, "o.faa"), "w") as f:
        f.write(">R1_0\nABC\n>R2_0\nDEF\n")
    reverse(os.path.join(d, "o.faa"), os.path.join(d, "r.faa"))
    assert open(os.path.join(d, "r.faa")).read() == ">R1_0\nCBA\n>R2_0\nFED\n"


if __name__ == "__main__":
    for name, f in list(globals().items()):
        if name.startswith("test_"):
            f()
    print("ok")
