# Wake 255 / Q34: is the half-step T=(I+C)^{1/2} (prove254/tmat.py) the shifted Murnaghan-Nakayama rule?
# Schur Q-functions built from marked shifted tableaux, one variable at a time:
#   Q_{lam/mu}(x) = 2^{#components of lam/mu} x^{|lam/mu|} if lam/mu is a (shifted) strip: lam_1>=mu_1>=lam_2>=mu_2>=..., else 0
# (Macdonald III (8.15)). MN coefficient c(kappa,rho;r) := [Q_kappa] p_r Q_rho, r odd, by triangular solve.
import sys, itertools, os
from fractions import Fraction as Fr
from functools import lru_cache
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '../prove254'))

def strict_parts(n, maxpart=None):
    out = []
    def rec(rem, mx, cur):
        if rem == 0: out.append(tuple(cur)); return
        for x in range(min(rem, mx), 0, -1): rec(rem - x, x - 1, cur + [x])
    rec(n, n if maxpart is None else maxpart, []); return out
def boxes(k): return frozenset((i, i + j) for i, li in enumerate(k) for j in range(li))
def ncomp(S):
    S = set(S); n = 0
    while S:
        n += 1; st = [S.pop()]
        while st:
            r, c = st.pop()
            for q in ((r+1, c), (r-1, c), (r, c+1), (r, c-1)):
                if q in S: S.remove(q); st.append(q)
    return n
def interlace(lam, mu):
    L = max(len(lam), len(mu)); a = list(lam) + [0]*(L+1-len(lam)); b = list(mu) + [0]*(L+1-len(mu))
    return all(a[i] >= b[i] >= a[i+1] for i in range(L))

NV = int(sys.argv[1]) if len(sys.argv) > 1 else 7   # number of variables = max degree used
@lru_cache(None)
def Q(lam, nv):  # dict exponent-tuple(len nv) -> int
    if nv == 0: return {(): 1} if not lam else {}
    out = {}
    n = sum(lam)
    for s in range(0, n + 1):
        for mu in strict_parts(s):
            if len(mu) > len(lam) or not boxes(mu) <= boxes(lam) or not interlace(lam, mu): continue
            d = n - s; w = 2 ** ncomp(boxes(lam) - boxes(mu))
            for e, c in Q(mu, nv - 1).items():
                out[e + (d,)] = out.get(e + (d,), 0) + c * w
    return {e: c for e, c in out.items() if c}
def mul_p(f, r, nv):
    out = {}
    for e, c in f.items():
        for i in range(nv):
            e2 = list(e); e2[i] += r; e2 = tuple(e2); out[e2] = out.get(e2, 0) + c
    return out
def expand_Q(f, nv):  # f symmetric, in Q-basis; requires nv >= degree
    f = {e: Fr(c) for e, c in f.items() if c}; res = {}
    while f:
        lead = max((e for e in f if list(e) == sorted(e, reverse=True)))
        mu = tuple(x for x in lead if x)
        assert len(set(mu)) == len(mu), ('non-strict leading term', lead)
        q = Q(mu, nv); c = f[lead] / q[lead]; res[mu] = c
        for e, v in q.items():
            f[e] = f.get(e, 0) - c * v
            if f[e] == 0: del f[e]
    return res
@lru_cache(None)
def MN(rho, r):
    n = sum(rho) + r
    return expand_Q(mul_p(Q(rho, n), r, n), n)
if __name__ == '__main__':
    # sanity: p_1 Q_rho
    print('p1 Q_(2,1):', MN((2, 1), 1)); print('p3 Q_():', MN((), 3)); print('p3 Q_(1):', MN((1,), 3))
