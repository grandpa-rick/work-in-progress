# Wake 251 compute: order polynomials of shifted-shape posets vs H_kappa(j)=chi(X_kappa,O(j)) on LG(m,2m).
# Omega(P,t) = sum_k c_k binom(t,k), c_k = # strict chains emptyset=I0<...<Ik=P in J(P).
import sys, itertools
from fractions import Fraction as Fr
from functools import lru_cache
import sympy as sp
t = sp.Symbol('t')

def ideals(elts, less):
    # elts: list; less(a,b) strict order. Return list of frozensets (order ideals).
    n = len(elts); below = {a: frozenset(b for b in elts if less(b, a)) for a in elts}
    res = [frozenset()]; seen = {frozenset()}; frontier = [frozenset()]
    while frontier:
        new = []
        for I in frontier:
            for a in elts:
                if a not in I and below[a] <= I:
                    J = I | {a}
                    if J not in seen: seen.add(J); new.append(J)
        res += new; frontier = new
    return res

def omega_coeffs(elts, less):
    J = ideals(elts, less); full = frozenset(elts); n = len(elts)
    Jset = sorted(J, key=len)
    # f[I][k] = # strict chains from empty to I of length k
    f = {frozenset(): {0: 1}}
    for I in Jset[1:]:
        d = {}
        # predecessors: all ideals I' strictly contained in I
        for Ip in Jset:
            if len(Ip) >= len(I): break
            if Ip <= I:
                for k, v in f[Ip].items(): d[k + 1] = d.get(k + 1, 0) + v
        f[I] = d
    return f[full] if n else {0: 1}

def omega_poly(elts, less):
    cs = omega_coeffs(elts, less)
    return sp.expand(sum(v * sp.binomial(t, k) for k, v in cs.items()).rewrite(sp.factorial).expand(func=True)) if False else \
        sp.expand(sum(v * sp.ff(t, k) / sp.factorial(k) for k, v in cs.items()))

def shifted(lam):
    return [(r, r + c) for r, l in enumerate(lam) for c in range(l)]
def young(lam):
    return [(r, c) for r, l in enumerate(lam) for c in range(l)]
def boxless(a, b):
    return a != b and a[0] <= b[0] and a[1] <= b[1]

def longest_chain(elts, less):
    @lru_cache(None)
    def h(a): return 1 + max([h(b) for b in elts if less(b, a)], default=0)
    return max([h(a) for a in elts], default=0)

def maximal_chain_lengths(elts, less):
    mins = [a for a in elts if not any(less(b, a) for b in elts)]
    def covers(a): return [b for b in elts if less(a, b) and not any(less(a, c) and less(c, b) for c in elts)]
    lens = set()
    def go(a, L):
        cs = covers(a)
        if not cs: lens.add(L); return
        for b in cs: go(b, L + 1)
    for a in mins: go(a, 1)
    return lens
