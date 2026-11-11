"""GOR secondary structure prediction: GOR IV scoring, with GOR V's profile averaging and decision
constants. Triplets, GOR V's other addition, are left out.

A residue's state S (H, E, C) comes from a 17-residue window (Garnier 1996):
    2/(W+1) * sum over offset pairs log f(S,a,b)/f(!S,a,b)
    - (W-1)/(W+1) * sum over offsets log f(S,a)/f(!S,a) + log f(S)/f(!S)
Sequences use the tables directly. Profiles use the same sum's expected value over each column's
emissions, assuming independent columns, which is GOR V's averaging over an alignment without one.

Usage: gor.py check
"""

import sys

import numpy as np

W = 17
HALF = W // 2
AA = "ACDEFGHIKLMNPQRSTVWY"
STATES = "HEC"
CODE = np.full(256, 20, dtype=np.int8)
for _i, _c in enumerate(AA):
    CODE[ord(_c)] = _i
SCODE = np.full(256, 2, dtype=np.int8)
for _i, _c in enumerate(STATES):
    SCODE[ord(_c)] = _i
PAIRS = [(a, b) for a in range(W) for b in range(a + 1, W)]


def encode(seq):
    return CODE[np.frombuffer(seq.encode(), dtype=np.uint8)]


def encode_ss(ss):
    return SCODE[np.frombuffer(ss.encode(), dtype=np.uint8)].astype(np.int64)


def windows(x):
    """N x W residue codes; 20 past either end, which has no information."""
    padded = np.concatenate([np.full(HALF, 20, x.dtype), x, np.full(HALF, 20, x.dtype)])
    return np.lib.stride_tricks.sliding_window_view(padded, W)


class Gor:
    """Log-odds tables, trained from (residue codes, state codes) pairs."""

    def __init__(self, chains, pseudo=1.0):
        win = np.concatenate([windows(s) for s, _ in chains]).astype(np.int64)
        st = np.concatenate([t for _, t in chains])
        keep = (win != 20).all(1)  # whole windows only, so every table counts the same residues
        win, st = win[keep], st[keep]
        self.prior = self._odds(np.bincount(st, minlength=3).astype(float))
        self.one = np.zeros((W, 21, 3))
        self.two = np.zeros((len(PAIRS), 21, 21, 3))
        for pos in range(W):
            c = np.bincount(win[:, pos] * 3 + st, minlength=60).reshape(20, 3).astype(float)
            self.one[pos, :20] = self._odds(c + pseudo)
        for i, (a, b) in enumerate(PAIRS):
            c = np.bincount((win[:, a] * 20 + win[:, b]) * 3 + st, minlength=1200)
            self.two[i, :20, :20] = self._odds(c.reshape(20, 20, 3).astype(float) + pseudo / 20)
        self.const = np.zeros(3)

    @staticmethod
    def _odds(c):
        return np.log(c / (c.sum(-1, keepdims=True) - c))

    def score_windows(self, win):
        """N x 3 scores of N x W residue codes."""
        win = win.astype(np.int64)
        two = np.zeros((len(win), 3))
        one = np.zeros((len(win), 3))
        for i, (a, b) in enumerate(PAIRS):
            two += self.two[i, win[:, a], win[:, b]]
        for pos in range(W):
            one += self.one[pos, win[:, pos]]
        return self._combine(two, one, (win != 20).sum(1))

    def _combine(self, two, one, n):
        """The formula with n, the residues actually in the window. A fixed W at a truncated window
        leaves the prior uncancelled and the singlets over-subtracted."""
        n = np.asarray(n, dtype=float)[:, None]
        return self.prior + 2 / (n + 1) * two - (n - 1) / (n + 1) * one

    def score_seq(self, x):
        """L x 3 scores of one sequence."""
        return self.score_windows(windows(x))

    def predict_seq(self, x):
        return (self.score_seq(x) + self.const).argmax(1)

    def score_profile(self, p):
        """L x 3 scores of an L x 20 emission profile, expected over the emissions."""
        n = len(p)
        q = np.zeros((n + W - 1, 21))
        q[HALF:HALF + n, :20] = p
        q[:HALF, 20] = q[HALF + n:, 20] = 1  # past the ends: the zero rows of the tables
        win = np.lib.stride_tricks.sliding_window_view(q, W, axis=0)  # L x 21 x W
        two = np.zeros((n, 3))
        one = np.zeros((n, 3))
        for i, (a, b) in enumerate(PAIRS):
            # L x 21 x 3 @ L x 21: contract one column, then the other
            part = np.einsum("la,abs->lbs", win[:, :, a], self.two[i], optimize=True)
            two += np.einsum("lbs,lb->ls", part, win[:, :, b], optimize=True)
        for pos in range(W):
            one += win[:, :, pos] @ self.one[pos]
        return self._combine(two, one, W - win[:, 20, :].sum(1))

    def predict_profile(self, p):
        return (self.score_profile(p) + self.const).argmax(1)


def fit_constants(scores, truth, grid=np.arange(-1.5, 1.51, 0.1)):
    """GOR V's decision constants: offsets for H and E that maximise Q3, relative to C."""
    best, out = -1.0, np.zeros(3)
    for h in grid:
        for e in grid:
            c = np.array([h, e, 0.0])
            q = ((scores + c).argmax(1) == truth).mean()
            if q > best:
                best, out = q, c
    return out


def check():
    """Planted signal: A favours H, V favours E, G favours C. Also profile against sequence."""
    rng = np.random.default_rng(1)
    chains = []
    for _ in range(300):
        st = np.repeat(rng.integers(0, 3, 30), rng.integers(4, 12, 30))
        favoured = np.array([0, 17, 5])[st]  # A, V, G
        x = np.where(rng.random(len(st)) < 0.5, favoured, rng.integers(0, 20, len(st)))
        chains.append((x.astype(np.int8), st.astype(np.int64)))
    g = Gor(chains[:250])
    truth = np.concatenate([t for _, t in chains[250:]])
    pred = np.concatenate([g.predict_seq(s) for s, _ in chains[250:]])
    q3 = (pred == truth).mean()
    assert q3 > 0.7, q3
    x = chains[260][0]
    p = np.eye(20)[x]  # a one-residue profile must score as the sequence does
    assert np.allclose(g.score_profile(p), g.score_seq(x)), "profile != sequence"
    assert len(g.predict_seq(x[:5])) == 5  # shorter than the window
    edge = centre = total = 0
    for s, t in chains[250:]:
        p = g.predict_seq(s)
        edge += (p[:3] == t[:3]).sum()
        centre += (p[8:-8] == t[8:-8]).sum()
        total += len(t[8:-8])
    assert edge / (3 * 50) > 0.8 * centre / total, "truncated windows must stay informative"
    print(f"gor check ok, planted Q3 {q3:.3f}")


if __name__ == "__main__":
    if sys.argv[1:] == ["check"]:
        check()
