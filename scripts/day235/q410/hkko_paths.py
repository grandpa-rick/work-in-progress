#!/usr/bin/env python3
"""(a) Cross-check the determinant engine (warnaar_calpha.py) against HKKO Thm 3.3 at t=1, in the
lattice-path model (HKKO Prop 2.7, as transcribed by Warnaar §6 / Clio §3): paths from 0 in
R = {x1>=...>=xh>=x1-w}, step a a 0/1 vector with alpha_a ones; signs c^-_{h,w}(lambda) as in
Warnaar §6 transposed to row-strict shape (checked by hand in review).
(b) Transport each path step to Korff's literal (lambda, mu, d) and evaluate Korff's (5.13) with the
code of korff_vs_macdonald.py; dictionary mu'_i = x_i - x_1 + w (our own derivation: Korff fixes
mu'_1 = k_K), d = theta_1.  Compare with Clio's formula (8).
(c) nu(P) = n at alpha=(1^n); sum_P c^- Psi_P vs c_alpha(t).
"""
import itertools
import sympy as sp
from warnaar_calpha import calpha, ord1, t
from korff_vs_macdonald import korff


def conj(p, length):
    return [sum(1 for x in p if x >= j) for j in range(1, length + 1)]


def cminus(lam, k, w):
    h = 2 * k
    if all(lam[2 * i] == lam[2 * i + 1] for i in range(k)):
        return 1
    if lam[0] - lam[h - 1] == w and all(lam[2 * i - 1] == lam[2 * i] for i in range(1, k)):
        return -1
    return 0


def paths(k, l, alpha):
    h, w = 2 * k, 2 * l + 2
    def rec(x, a, steps):
        if a == len(alpha):
            yield x, steps; return
        for S in itertools.combinations(range(h), alpha[a]):
            th = tuple(1 if i in S else 0 for i in range(h))
            y = tuple(x[i] + th[i] for i in range(h))
            if all(y[i] >= y[i + 1] for i in range(h - 1)) and y[h - 1] >= y[0] - w:
                yield from rec(y, a + 1, steps + [(x, th)])
    yield from rec(tuple([0] * h), 0, [])


def psi8(x, th, w):
    h = len(x); out = 1; cnt = 0
    for i in range(h):
        ip = (i + 1) % h
        if th[i] == 0 and th[ip] == 1:
            m = x[i] - x[i + 1] if i < h - 1 else w - (x[0] - x[h - 1])
            assert m >= 1
            out *= (1 - t ** m); cnt += 1
    return out, cnt


def korff_literal(x, th, w):
    h = len(x)
    mup = [x[i] - x[0] + w for i in range(h)]
    y = [x[i] + th[i] for i in range(h)]
    lamp = [y[i] - y[0] + w for i in range(h)]
    mu = tuple(conj(mup, w)); lam = tuple(conj(lamp, w))  # conj of conjugate = partition (k_K=w parts)
    r = korff(lam, mu, th[0], h, w)
    assert r is not None
    thK = tuple(r[0][i] for i in range(1, h + 1))
    assert thK == tuple(th), (x, th, thK)
    return r[3]


if __name__ == '__main__':
    print("(a) c_alpha(1) [determinant] vs HKKO signed path count, all alpha in {0..h}^n")
    for (k, l, n) in [(1, 0, 3), (1, 0, 4), (1, 1, 4), (1, 1, 5), (2, 0, 3), (2, 0, 4)]:
        h = 2 * k; bad = 0; tot = 0
        for alpha in itertools.product(range(h + 1), repeat=n):
            s = sum(cminus(xf, k, 2 * l + 2) for xf, _ in paths(k, l, list(alpha)))
            c1 = sp.expand(calpha(k, l, list(alpha))).subs(t, 1)
            tot += 1; bad += (s != c1)
        print(f"  (k,l,n)=({k},{l},{n}): {tot} contents, mismatches {bad}")
    print("(b),(c) k=1, alpha=(1^n), n=2M")
    for (k, l) in [(1, 0), (1, 1), (1, 2), (2, 0), (2, 1)]:
        M = k + l + 1; w = 2 * l + 2; n = 2 * M if k == 1 else 2 * M - 2
        alpha = [1] * n
        nus = set(); S = 0; kb = 0; cnt = 0
        for xf, steps in paths(k, l, alpha):
            P = 1; nu = 0; PK = 1
            for (x, th) in steps:
                f, c = psi8(x, th, w); P *= f; nu += c
                PK *= korff_literal(x, th, w)
            if sp.expand(P - PK) != 0:
                kb += 1
            nus.add(nu); cnt += 1
            S += cminus(xf, k, w) * P
        c = calpha(k, l, alpha)
        print(f"  (k,l,n)=({k},{l},{n}): paths {cnt}, nu values {sorted(nus)}, (8) vs Korff-literal mismatches {kb}; "
              f"sum c^- Psi = {sp.factor(S)}  vs  c = {sp.factor(c)} (ord {ord1(c)})")
