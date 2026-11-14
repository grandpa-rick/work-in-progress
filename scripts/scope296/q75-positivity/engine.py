#!/usr/bin/env python3
"""Exact Schur-basis engine for G = Sym/Ann(R), R = ribbon r_(2^k), Hall pairing, f acts by f^perp.
Image of f in G_i  <->  functional g |-> <f g, R> on Sym_{d-i}; we store it as the vector
(<f s_nu, R>)_{nu |- d-i}.  Ann = kernel of this map (definition of Q75 Corollary / q75_ring.py).
Integer LR coefficients via h-basis + Kostka.  Grade: computed.
"""
import itertools
from functools import lru_cache
from fractions import Fraction

@lru_cache(None)
def partitions(n, maxpart=None):
    if maxpart is None: maxpart = n
    if n == 0: return ((),)
    out = []
    for f in range(min(n, maxpart), 0, -1):
        for rest in partitions(n - f, f):
            out.append((f,) + rest)
    return tuple(out)

@lru_cache(None)
def kostka(lam, mu):
    """# SSYT shape lam content mu (mu any composition)."""
    if sum(lam) != sum(mu): return 0
    if not mu: return 1 if not lam else 0
    m = mu[-1]; rest = mu[:-1]
    # remove horizontal strip of size m from lam
    tot = 0
    L = len(lam)
    def rec(i, rem, new):
        nonlocal tot
        if i == L:
            if rem == 0:
                nl = tuple(x for x in new if x > 0)
                tot += kostka(nl, rest)
            return
        nxt = lam[i + 1] if i + 1 < L else 0
        for r in range(0, min(rem, lam[i] - nxt) + 1):
            rec(i + 1, rem - r, new + [lam[i] - r])
    rec(0, m, [])
    return tot

@lru_cache(None)
def Kmats(n):
    P = partitions(n)
    K = [[kostka(l, m) for m in P] for l in P]  # K[l][m]: h_m = sum_l K[l][m] s_l
    # invert (unitriangular w.r.t. dominance; use Fractions gaussian)
    N = len(P)
    A = [[Fraction(K[i][j]) for j in range(N)] + [Fraction(int(i == j)) for j in range(N)] for i in range(N)]
    for c in range(N):
        p = next(r for r in range(c, N) if A[r][c] != 0)
        A[c], A[p] = A[p], A[c]
        pv = A[c][c]; A[c] = [x / pv for x in A[c]]
        for r in range(N):
            if r != c and A[r][c] != 0:
                f = A[r][c]; A[r] = [a - f * b for a, b in zip(A[r], A[c])]
    Kinv = [[int(A[i][N + j]) for j in range(N)] for i in range(N)]
    return P, K, Kinv  # s_l = sum_m Kinv'?  see s_to_h

@lru_cache(None)
def s_to_h(lam):
    n = sum(lam); P, K, Kinv = Kmats(n)
    # h = K^T-ish: h_m = sum_l K[l][m] s_l  => vector of s in terms of h: s = (K)^{-1} applied: s_l = sum_m Kinv[m][l] h_m
    i = P.index(lam)
    return {P[m]: Kinv[m][i] for m in range(len(P)) if Kinv[m][i] != 0}

def h_to_s(hdict):
    out = {}
    for mu, c in hdict.items():
        n = sum(mu); P, K, _ = Kmats(n)
        for li, l in enumerate(P):
            v = K[li][P.index(mu)] if True else 0
            if v: out[l] = out.get(l, 0) + c * v
    return {k: v for k, v in out.items() if v}

@lru_cache(None)
def lr(lam, mu):
    """s_lam s_mu in Schur basis, dict."""
    a, b = s_to_h(lam), s_to_h(mu)
    prod = {}
    for m1, c1 in a.items():
        for m2, c2 in b.items():
            m = tuple(sorted(m1 + m2, reverse=True))
            prod[m] = prod.get(m, 0) + c1 * c2
    return h_to_s(prod)

def ribbon(alpha):
    L = len(alpha); hd = {}
    for cuts in itertools.product([0, 1], repeat=L - 1):
        beta, cur = [], alpha[0]
        for j, c in enumerate(cuts):
            if c: cur += alpha[j + 1]
            else: beta.append(cur); cur = alpha[j + 1]
        beta.append(cur)
        m = tuple(sorted(beta, reverse=True))
        hd[m] = hd.get(m, 0) + (-1) ** (L - len(beta))
    return h_to_s(hd)

# ---------- linear algebra over Q ----------
def rref(rows):
    M = [[Fraction(x) for x in r] for r in rows]
    piv = []; r = 0
    ncol = len(M[0]) if M else 0
    for c in range(ncol):
        p = next((i for i in range(r, len(M)) if M[i][c] != 0), None)
        if p is None: continue
        M[r], M[p] = M[p], M[r]
        pv = M[r][c]; M[r] = [x / pv for x in M[r]]
        for i in range(len(M)):
            if i != r and M[i][c] != 0:
                f = M[i][c]; M[i] = [a - f * b for a, b in zip(M[i], M[r])]
        piv.append(c); r += 1
    return M[:r], piv

def rank(rows):
    return len(rref(rows)[0]) if rows else 0

def solve(basis_rows, v):
    """coefficients x with sum x_i basis_rows[i] = v (basis independent); None if not in span."""
    n = len(basis_rows)
    # augmented: columns = coordinates; solve B^T x = v
    m = len(v)
    A = [[Fraction(basis_rows[i][j]) for i in range(n)] + [Fraction(v[j])] for j in range(m)]
    sol, piv = rref(A)
    if n in piv: return None
    x = [Fraction(0)] * n
    for row, c in zip(sol, piv): x[c] = row[n]
    return x

class G:
    def __init__(self, R, d):
        self.R, self.d = R, d
        self._img = {}
    def pair(self, f):  # <f,R> for f dict in Schur basis
        return sum(c * self.R.get(l, 0) for l, c in f.items())
    def img(self, lam):
        """vector of s_lam in G_i : (<s_lam s_nu, R>)_{nu |- d-i}"""
        if lam not in self._img:
            i = sum(lam)
            self._img[lam] = tuple(self.pair(lr(lam, nu)) for nu in partitions(self.d - i)) if i <= self.d else ()
        return self._img[lam]
    def img_of(self, f, i):
        v = [0] * len(partitions(self.d - i))
        for l, c in f.items():
            for j, x in enumerate(self.img(l)): v[j] += c * x
        return v
