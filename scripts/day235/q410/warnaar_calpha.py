#!/usr/bin/env python3
"""Independent computation of c_alpha(t) = [x^alpha] (x1..xn)^k G_{k,l,n}(x;t) from Warnaar
arXiv:2511.17034 (6.1):
   G = sum_{y in Z^k} det_{1<=i,j<=k} ( t^{M y_i^2 - j y_i} ( e'_{n-i+j-N y_i} - e'_{n+i+j-N y_i} ) ),
   M = k+l+1, N = 2M, e'_r = e_r(x1, 1/x1, ..., xn, 1/xn).
Coefficient extraction: [x^beta] prod_s e'_{r_s} = [u^r] prod_v [x^{beta_v}] prod_s (1+u_s x)(1+u_s/x).
Usage: python3 warnaar_calpha.py
"""
import itertools, functools, sys
import sympy as sp
import sympy.combinatorics
from math import comb

t = sp.symbols('t')


def polymul(A, B, caps):
    C = {}
    for ea, ca in A.items():
        for eb, cb in B.items():
            e = tuple(x + y for x, y in zip(ea, eb))
            if any(x > c for x, c in zip(e, caps)):
                continue
            C[e] = C.get(e, 0) + ca * cb
    return {e: c for e, c in C.items() if c}


def per_var(k, b, caps):
    # [x^b] prod_s (1+u_s x)(1+u_s/x): each s contributes u_s^{a_s} x^{c_s}, (a,c) in {(0,0),(1,1),(1,-1),(2,0)}
    opts = [(0, 0, 1), (1, 1, 1), (1, -1, 1), (2, 0, 1)]
    out = {}
    for ch in itertools.product(opts, repeat=k):
        if sum(c[1] for c in ch) != b:
            continue
        e = tuple(c[0] for c in ch)
        out[e] = out.get(e, 0) + 1
    return out


def calpha(k, l, alpha):
    n = len(alpha); M = k + l + 1; N = 2 * M
    caps = tuple([2 * n] * k)
    P = {tuple([0] * k): 1}
    for a in alpha:
        P = polymul(P, per_var(k, a - k, caps), caps)
    E = lambda rs: 0 if any(r < 0 or r > 2 * n for r in rs) else P.get(tuple(rs), 0)
    total = 0
    # y_i range: need some index in [0,2n]: n - i + j - N y in [-(k), 2n+2k] roughly
    yr = range(-(2 * n // N) - 2, (2 * n // N) + 3)
    for y in itertools.product(yr, repeat=k):
        for sig in itertools.permutations(range(1, k + 1)):
            sgn = sp.combinatorics.Permutation([s - 1 for s in sig]).signature()
            texp = sum(M * y[i] ** 2 - sig[i] * y[i] for i in range(k))
            # expand prod_i (e_{a_i} - e_{b_i})
            acc = 0
            for choice in itertools.product((0, 1), repeat=k):
                rs = []
                for i in range(k):
                    ii = i + 1
                    a = n - ii + sig[i] - N * y[i]; b = n + ii + sig[i] - N * y[i]
                    rs.append(b if choice[i] else a)
                acc += (-1) ** sum(choice) * E(rs)
            if acc:
                total += sgn * acc * t ** texp
    return sp.expand(total)


def ord1(p):
    p = sp.Poly(sp.expand(p * t ** 50), t)  # clear negative powers
    if p.is_zero:
        return None
    o = 0
    while p.eval(1) == 0:
        p = sp.Poly(sp.quo(p.as_expr(), t - 1), t); o += 1
    return o


if __name__ == '__main__':
    import sympy.combinatorics  # noqa
    print("=== k=1, alpha=(1^n), n=N=2M ; claim c = C_M - (2M-1) t^{M-1} + t^{M+1}")
    for M in range(2, 7):
        l = M - 2; n = 2 * M
        c = calpha(1, l, [1] * n)
        CM = comb(2 * M, M) // (M + 1)
        claim = sp.expand(CM - (2 * M - 1) * t ** (M - 1) + t ** (M + 1))
        print(f"M={M} (l={l}, n={n}): c = {c}   matches claim: {sp.expand(c-claim)==0}   ord_(t=1) = {ord1(c)}  (Clio's bound would need >= {n})")
    print("=== k>=2 (Korff's own domain n_K=2k>2), alpha=(1^n): ord_(t=1) c vs required n")
    for (k, l) in [(2, 0), (2, 1), (2, 2), (3, 0), (3, 1)]:
        M = k + l + 1; N = 2 * M
        for n in range(max(2, N - 2), N + 3):
            c = calpha(k, l, [1] * n)
            print(f"k={k} l={l} n={n}: c = {c}   ord={ord1(c)}  required>={n}")
