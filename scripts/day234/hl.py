# Day 234 independent symmetric-function toolkit (fresh code): transitions among m, p, e, s, Hall-Littlewood P,
# Green polynomials X^lambda_mu(t) (p_mu = sum_lambda X^lambda_mu P_lambda), Kostka-Foulkes K(t) (s_lam = sum K_{lam,mu}(t) P_mu).
# HL P built by Gram-Schmidt from <p_l,p_m>_t = delta z_l prod (1-t^{l_i})^{-1}, unitriangular in m (Macdonald III (2.6)/(4.9)).
from functools import lru_cache
from math import factorial
from collections import Counter
from sympy.polys.fields import field
from sympy import QQ
from engine import partitions

F, tt = field("t", QQ)

def z(lam):
    c = Counter(lam); r = 1
    for k, m in c.items(): r *= k**m * factorial(m)
    return r

def nfun(lam): return sum(i * l for i, l in enumerate(lam))

def conj(lam):
    return tuple(sum(1 for l in lam if l > i) for i in range(lam[0])) if lam else ()

def dominates(mu, lam):  # mu >= lam
    if sum(mu) != sum(lam): return False
    a = b = 0
    for i in range(max(len(mu), len(lam))):
        a += mu[i] if i < len(mu) else 0
        b += lam[i] if i < len(lam) else 0
        if a < b: return False
    return True

def p_in_m(lam, mu):
    """coefficient of m_mu in p_lam: # ways to put parts of lam into rows with row sums mu"""
    mu = list(mu)
    def rec(i, rem):
        if i == len(lam): return 1 if all(r == 0 for r in rem) else 0
        tot = 0
        for j in range(len(rem)):
            if rem[j] >= lam[i]:
                rem[j] -= lam[i]; tot += rec(i + 1, rem); rem[j] += lam[i]
        return tot
    return rec(0, mu)

def e_in_m(lam, mu):
    """# 0-1 matrices with row sums lam, column sums mu"""
    mu = list(mu)
    import itertools
    def rec(i, rem):
        if i == len(lam): return 1 if all(r == 0 for r in rem) else 0
        tot = 0
        for S in itertools.combinations(range(len(rem)), lam[i]):
            if all(rem[j] > 0 for j in S):
                for j in S: rem[j] -= 1
                tot += rec(i + 1, rem)
                for j in S: rem[j] += 1
        return tot
    return rec(0, mu)

def kostka(lam, mu):
    """# SSYT of shape lam, content mu (remove horizontal strips of the largest letter)"""
    @lru_cache(None)
    def rec(shape, k):
        if k == 0: return 1 if sum(shape) == 0 else 0
        c = mu[k - 1]; tot = 0
        # remove horizontal strip of size c from shape
        def strips(i, sh, left):
            if i == len(sh):
                if left == 0: yield tuple(sh)
                return
            nxt = sh[i + 1] if i + 1 < len(sh) else 0
            for r in range(0, min(left, sh[i] - nxt) + 1):
                sh2 = list(sh); sh2[i] -= r
                yield from strips(i + 1, sh2, left - r)
        for sh in set(strips(0, list(shape), c)):
            tot += rec(sh, k - 1)
        return tot
    return rec(tuple(lam), len(mu))

class SymN:
    def __init__(self, n):
        self.n = n
        self.B = list(partitions(n))  # decreasing lex
        self.idx = {b: i for i, b in enumerate(self.B)}
        from sympy import Matrix
        self.Pm = Matrix([[p_in_m(l, m) for m in self.B] for l in self.B])  # p_l = sum Pm[l,m] m_m
        self.Em = Matrix([[e_in_m(l, m) for m in self.B] for l in self.B])
        self.Sm = Matrix([[kostka(l, m) for m in self.B] for l in self.B])
        self.Mp = self.Pm.inv()  # m_l = sum Mp[l,q] p_q
        self._hl = None

    def to_p(self, vec_m):
        """vec_m: dict partition->coef in m-basis; returns dict in p basis"""
        out = {}
        for l, c in vec_m.items():
            i = self.idx[l]
            for j, q in enumerate(self.B):
                if self.Mp[i, j] != 0:
                    out[q] = out.get(q, 0) + c * F(self.Mp[i, j])
        return out

    def e_to_m(self, vec_e):
        out = {}
        for l, c in vec_e.items():
            i = self.idx[l]
            for j, q in enumerate(self.B):
                if self.Em[i, j] != 0: out[q] = out.get(q, 0) + c * self.Em[i, j]
        return out

    def ip_t(self, a, b):  # a,b dicts in p basis (field elements)
        tot = F(0)
        for q, c in a.items():
            if q in b:
                w = F(z(q))
                for qi in q: w = w / (1 - tt**qi)
                tot += c * b[q] * w
        return tot

    def ip_hall(self, a, b):
        return sum((c * b[q] * z(q) for q, c in a.items() if q in b), F(0))

    def hl(self):
        """P_lambda in m basis: dict lam -> {mu: coef in F}"""
        if self._hl: return self._hl
        order = list(reversed(self.B))  # 1^n first (increasing revlex)
        P = {}; Pp = {}; norms = {}
        for lam in order:
            v = {lam: F(1)}
            vp = self.to_p({lam: F(1)})
            for mu in list(P.keys()):
                coef = self.ip_t(vp, Pp[mu]) / norms[mu]
                if coef != 0:
                    for k2, c2 in P[mu].items(): v[k2] = v.get(k2, F(0)) - coef * c2
                    for k2, c2 in Pp[mu].items(): vp[k2] = vp.get(k2, F(0)) - coef * c2
            P[lam] = {k: c for k, c in v.items() if c != 0}
            Pp[lam] = {k: c for k, c in vp.items() if c != 0}
            norms[lam] = self.ip_t(Pp[lam], Pp[lam])
        self._hl = P; self._hlp = Pp; self._norms = norms
        return P

    def expand_in_P(self, vec_m):
        """write vec (m-basis, field coefs) as sum c_lam P_lam using unitriangularity (revlex)"""
        P = self.hl()
        v = dict(vec_m); out = {}
        for lam in self.B:  # decreasing revlex: leading term first
            c = v.get(lam, F(0))
            if c != 0:
                out[lam] = c
                for k2, c2 in P[lam].items(): v[k2] = v.get(k2, F(0)) - c * c2
        assert all(c == 0 for c in v.values())
        return out

    def green(self):
        """X[lam][mu] with p_mu = sum_lam X^lam_mu P_lam"""
        X = {}
        for mu in self.B:
            vm = {q: F(self.Pm[self.idx[mu], j]) for j, q in enumerate(self.B) if self.Pm[self.idx[mu], j] != 0}
            for lam, c in self.expand_in_P(vm).items():
                X.setdefault(lam, {})[mu] = c
        return X

    def kostka_foulkes(self):
        """K[lam][mu](t): s_lam = sum_mu K_{lam mu}(t) P_mu"""
        K = {}
        for lam in self.B:
            vm = {q: F(self.Sm[self.idx[lam], j]) for j, q in enumerate(self.B) if self.Sm[self.idx[lam], j] != 0}
            K[lam] = self.expand_in_P(vm)
        return K

def b_lam(lam):
    r = F(1)
    for k, m in Counter(lam).items():
        for i in range(1, m + 1): r *= (1 - tt**i)
    return r
