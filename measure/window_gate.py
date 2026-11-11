"""Speed plan step 1: how wide must the window be, and what gate? From `filter window` output.

For each width: the highest gate (quarter bits) that keeps >= 99.5% of (query, model) pairs at each bar,
and the share of all hits whose window reaches it. Also how quarter-bit integer diagonals move pairs
across the bar. Usage: window_gate.py PREFIX
"""

import sys

import numpy as np

pre = sys.argv[1]
h = np.loadtxt(pre + ".hist.tsv", dtype=np.int64)
p = np.genfromtxt(pre + ".pairs.tsv", dtype=None, encoding=None)
d = np.array([r[2] for r in p], dtype=float)
dq = np.array([r[3] for r in p], dtype=float)
wb = np.array([[r[4 + i] for i in range(5)] for r in p], dtype=float)
widths = sorted(set(h[:, 0]))
total = h[h[:, 0] == widths[0], 2].sum()
print(f"{total:,} hits; pairs at diagonal >= 22: {(d >= 22).sum():,}, >= 24: {(d >= 24).sum():,}")
for bar in (22, 24):
    at = d >= bar + 0.005  # printed to 2 decimals: 21.996 shows as 22.00
    print(f"\nbar {bar}: width, gate keeping >= 99.5% / 99.9% of pairs, hits surviving each; hits surviving at gates 12-20 bits")
    for w, width in enumerate(widths):
        tail = h[h[:, 0] == width, 2][::-1].cumsum()[::-1] / total  # tail[b]: share with window >= b / 4 bits
        row = []
        for keep in (0.995, 0.999):
            gates = np.arange(0, 256)
            kept = np.array([(wb[at, w] >= g / 4).mean() for g in gates])
            g = gates[kept >= keep].max()
            row.append(f"{g / 4:5.2f} bits {tail[g]:.4f}")
        fixed = " ".join(f"{tail[int(g * 4)]:.4f}" for g in (12, 14, 16, 18, 20))
        print(f"  {width:3d}: " + "   ".join(row) + f"   | {fixed}")
for bar in (22, 24):
    a, b = d >= bar, dq >= bar
    print(f"quarter-bit diagonal at {bar}: {(a & ~b).sum()} pairs lost, {(b & ~a).sum()} gained of {a.sum()}; "
          f"largest difference {np.abs(d - dq).max():.2f} bits")
