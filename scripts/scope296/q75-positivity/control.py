#!/usr/bin/env python3
"""(d) Negative/positive controls.
 C1: LR engine vs dimension identity sum_nu c^nu_{lam,mu} f^nu = C(n,|lam|) f^lam f^mu, all |lam|+|mu|<=8,
     plus hand values s21*s21, s2*s11.
 C2: same quotient pipeline with R = s_(k,k): Sym/Ann(s_kk) = H*Gr(2,k+2); expect nonzero s_lam exactly
     lam in the 2 x k box, Schur images independent, and structure constants = truncated LR (>=0).
Grade: computed."""
from math import comb, factorial
from engine import *
from positivity import coords

def flam(l):  # hook length
    n = sum(l); conj = [sum(1 for x in l if x > j) for j in range(l[0])] if l else []
    h = 1
    for i, a in enumerate(l):
        for j in range(a): h *= (a - j - 1) + (conj[j] - i - 1) + 1
    return factorial(n) // h

bad = 0; cnt = 0
for n in range(2, 9):
    for a in range(1, n):
        for lam in partitions(a):
            for mu in partitions(n - a):
                p = lr(lam, mu); cnt += 1
                if min(p.values()) < 0 or sum(c * flam(nu) for nu, c in p.items()) != comb(n, a) * flam(lam) * flam(mu): bad += 1
print('C1 dimension identity: products checked', cnt, 'failures', bad)
print('   s21*s21 =', lr((2, 1), (2, 1)))
print('   s2*s11  =', lr((2,), (1, 1)))

for k in range(2, 6):
    d = 2 * k; g = G({(k, k): 1}, d)
    ok = True; neg = 0; tot = 0; lrmatch = True
    B = {}
    for n in range(d + 1):
        nz = [l for l in partitions(n) if any(g.img(l))]
        box = [l for l in partitions(n) if len(l) <= 2 and (not l or l[0] <= k)]
        if sorted(nz) != sorted(box) or rank([g.img(l) for l in nz]) != len(nz): ok = False
        B[n] = nz
    for n in range(2, d + 1):
        for i in range(1, n):
            for lam in B[i]:
                for mu in B[n - i]:
                    c = coords(g, B[n], g.img_of(lr(lam, mu), n)); tot += len(c)
                    neg += sum(1 for x in c if x < 0)
                    full = lr(lam, mu)
                    if any(x != full.get(nu, 0) for nu, x in zip(B[n], c)): lrmatch = False
    print(f'C2 k={k}: nonzero Schur = 2x{k} box & independent: {ok};  constants {tot}, negative {neg}, equal to LR: {lrmatch}')
