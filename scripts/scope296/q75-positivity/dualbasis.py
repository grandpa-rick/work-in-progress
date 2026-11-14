#!/usr/bin/env python3
"""(e) Alternative bases.
 (1) Ray basis: per degree pick one Schur image on each extreme ray when the cone is simplicial;
     report whether ALL Schur images are >=0 in it.
 (2) Poincare-dual basis: b^lam in G_n with <b^lam s_mu, R> = delta (mu in Schur basis B_{d-n}).
     For Gr(2,N) this is again Schur.  Compute its structure constants and sign.
 (3) Triple numbers <s_lam s_mu s_nu, R> are >= 0 automatically (R Schur-positive, LR >= 0).
Grade: computed."""
import sys, itertools
from engine import *
from positivity import setup, coords

def greedy(g, nz, n, pref):
    S = []
    for l in sorted(nz[n], key=pref):
        if rank([g.img(x) for x in S + [l]]) > len(S): S.append(l)
    return S

for k in map(int, sys.argv[1:]):
    d, g, nz, r = setup(k)
    # Schur basis B (lex-desc greedy) in each degree
    B = {n: greedy(g, nz, n, lambda l: tuple(-x for x in l)) for n in range(d + 1)}
    # dual basis of G_n: functional on G_{d-n}.  Represent elements of G_n by their values on B_{d-n}
    # (Gorenstein: this is an iso G_n -> Q^{r}).  s_lam |-> (<s_lam s_mu,R>)_{mu in B_{d-n}}.
    def val(f, n):  # f Schur dict, homogeneous deg n
        return [sum(c * g.pair(lr(l, mu)) for l, c in f.items()) for mu in B[d - n]]
    # dual basis b^mu in these coords = unit vectors e_mu.  Product b^a b^b: need representatives.
    rep = {}
    for n in range(d + 1):
        M = [val({l: 1}, n) for l in B[n]]   # rows: s_l in dual coords
        # b^mu = sum_l x_l s_l with sum_l x_l M[l] = e_mu
        rep[n] = [dict(zip(B[n], solve(M, [int(i == j) for j in range(r[n])]))) for i in range(r[n])]
    neg = []; tot = 0
    for n in range(2, d + 1):
        for i in range(1, n // 2 + 1):
            for a, fa in enumerate(rep[i]):
                for b, fb in enumerate(rep[n - i]):
                    prod = {}
                    for l1, c1 in fa.items():
                        for l2, c2 in fb.items():
                            for nu, c in lr(l1, l2).items(): prod[nu] = prod.get(nu, 0) + c1 * c2 * c
                    cs = val(prod, n)  # coefficients in dual basis = values on B_{d-n}
                    tot += len(cs)
                    for j, x in enumerate(cs):
                        if x < 0: neg.append(((i, B[d - i][a]), (n - i, B[d - n + i][b]), B[d - n][j], x))
    print('k =', k, 'Poincare-dual-of-Schur basis: #constants', tot, '#negative', len(neg))
    for t in neg[:6]: print('   b^', t[0], '* b^', t[1], ' coeff on b^', t[2], '=', t[3])
