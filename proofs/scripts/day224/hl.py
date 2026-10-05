"""Day 224 engine: symmetric functions at a numeric t (Fraction), HL P via Gram-Schmidt, parabolic symmetrizer
T_a(P_rho) = P_{rho+1^a} (Macdonald III (2.2), Day 223 Lemma 1.3). Everything stored in the p-basis.
Use: Sym(n, t) per degree. Grade of anything produced here: computed."""
from fractions import Fraction as Fr
from functools import lru_cache
from math import factorial
from collections import Counter
import itertools

def parts(n, m=None):
    if m is None: m = n
    if n == 0: yield (); return
    for k in range(min(n, m), 0, -1):
        for r in parts(n-k, k): yield (k,)+r

def zee(lam):
    c = Counter(lam); z = 1
    for k, m in c.items(): z *= k**m * factorial(m)
    return z

def dominates(a, b):  # a >= b
    sa = sb = 0
    for i in range(max(len(a), len(b))):
        sa += a[i] if i < len(a) else 0; sb += b[i] if i < len(b) else 0
        if sa < sb: return False
    return True

@lru_cache(None)
def p_to_m(lam):  # p_lam = sum_mu R[mu] m_mu ; R = #maps parts->rows of mu (exponent vector) with exact row sums
    n = sum(lam); out = {}
    for mu in parts(n):
        @lru_cache(None)
        def cnt(i, rem):
            if i == len(lam): return 1 if not any(rem) else 0
            tot = 0
            for j in range(len(rem)):
                if rem[j] >= lam[i]:
                    r = list(rem); r[j] -= lam[i]; tot += cnt(i+1, tuple(r))
            return tot
        c = cnt(0, mu)
        if c: out[mu] = Fr(c)
    return out

class Sym:
    def __init__(self, n, t):
        self.n, self.t = n, Fr(t)
        self.P = list(parts(n)); self.idx = {l: i for i, l in enumerate(self.P)}
        self._m = None; self._hl = None; self._e = None
    def ip(self, f, g):  # HL scalar product, f,g dicts in p basis
        t = self.t; s = Fr(0)
        for k, v in f.items():
            if k in g:
                w = Fr(zee(k))
                for r in k: w /= (1 - t**r)
                s += v*g[k]*w
        return s
    def m_basis(self):  # m_mu in p basis: invert p_to_m (triangular in dominance)
        if self._m: return self._m
        P = self.P; N = len(P)
        A = [[p_to_m(l).get(mu, Fr(0)) for mu in P] for l in P]  # row l: p_l = sum A[l][mu] m_mu
        # invert A
        inv = matinv(A)
        # m_mu = sum_l inv[mu][l] p_l
        self._m = {mu: {P[j]: inv[i][j] for j in range(N) if inv[i][j]} for i, mu in enumerate(P)}
        return self._m
    def hl(self):  # P_lam (HL) in p basis, Gram-Schmidt from bottom of dominance
        if self._hl: return self._hl
        m = self.m_basis(); P = self.P
        order = sorted(P, key=lambda l: l)  # lex increasing = linear extension of dominance (small first)
        hl = {}; norms = {}
        for lam in order:
            f = dict(m[lam])
            for mu in hl:
                if dominates(lam, mu) and mu != lam:
                    c = self.ip(m[lam], hl[mu]) / norms[mu]
                    if c:
                        for k, v in hl[mu].items(): f[k] = f.get(k, 0) - c*v
            f = {k: v for k, v in f.items() if v}
            hl[lam] = f; norms[lam] = self.ip(f, f)
        self._hl = hl; self._norms = norms
        return hl
    def to_hl(self, f):  # coefficients of f in P basis
        hl = self.hl(); return {l: self.ip(f, hl[l])/self._norms[l] for l in self.P if self.ip(f, hl[l])}
    def e_basis(self):
        if self._e: return self._e
        P = self.P; N = len(P)
        E = {}
        for lam in P:
            f = {(): Fr(1)}
            for k in lam: f = mul(f, e_p(k))
            E[lam] = f
        self._e = E
        M = [[E[l].get(mu, Fr(0)) for mu in P] for l in P]
        self._einv = matinv(M)  # p_mu = sum_l einv[mu][l] e_l
        return E
    def to_e(self, f):
        self.e_basis(); P = self.P; out = {}
        for i, mu in enumerate(P):
            v = f.get(mu)
            if v:
                for j, l in enumerate(P):
                    if self._einv[i][j]: out[l] = out.get(l, 0) + v*self._einv[i][j]
        return {k: v for k, v in out.items() if v}

def matinv(A):
    N = len(A); M = [list(r) + [Fr(int(i == j)) for j in range(N)] for i, r in enumerate(A)]
    for c in range(N):
        p = next(r for r in range(c, N) if M[r][c] != 0); M[c], M[p] = M[p], M[c]
        pv = M[c][c]; M[c] = [x/pv for x in M[c]]
        for r in range(N):
            if r != c and M[r][c] != 0:
                f = M[r][c]; M[r] = [x - f*y for x, y in zip(M[r], M[c])]
    return [r[N:] for r in M]

def mul(f, g):
    out = {}
    for a, x in f.items():
        for b, y in g.items():
            k = tuple(sorted(a+b, reverse=True)); out[k] = out.get(k, 0) + x*y
    return {k: v for k, v in out.items() if v}

def add(f, g, c=1):
    out = dict(f)
    for k, v in g.items(): out[k] = out.get(k, 0) + c*v
    return {k: v for k, v in out.items() if v}

@lru_cache(None)
def _e_p(k):
    out = {}
    for lam in parts(k):
        out[lam] = Fr((-1)**(k-len(lam)), zee(lam))
    return out
def e_p(k): return dict(_e_p(k)) if k > 0 else {(): Fr(1)}

_SYM = {}
def S(n, t):
    if (n, t) not in _SYM: _SYM[(n, t)] = Sym(n, t)
    return _SYM[(n, t)]

def T(a, f, t):
    """T_a(f) for f a homogeneous symmetric function (p-basis dict), restricted to a variables. Returns p-basis dict."""
    if not f: return {}
    d = sum(next(iter(f))); n = d + a
    if d == 0: c = f.get((), 0); return {k: c*v for k, v in S(a, t).hl()[(1,)*a].items()} if a else {(): c}
    co = S(d, t).to_hl(f); out = {}
    hl = S(n, t).hl()
    for rho, c in co.items():
        if len(rho) > a: continue
        lam = tuple(r+1 for r in rho) + (1,)*(a-len(rho))
        for k, v in hl[lam].items(): out[k] = out.get(k, 0) + c*v
    return {k: v for k, v in out.items() if v}
