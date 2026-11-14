#!/usr/bin/env python3
"""(b),(c),(e): structure constants of Schur-image bases of G_(2^k); exhaustive search for a positive Schur basis.
Coordinates: G_n is identified with the row space of the image vectors; we fix coordinates by an rref basis.
Grade: computed.
"""
import sys, itertools
from engine import *

def setup(k):
    d = 2 * k; R = ribbon((2,) * k); g = G(R, d)
    nz, r = {}, {}
    for n in range(d + 1):
        nz[n] = [l for l in partitions(n) if any(g.img(l))]
        r[n] = rank([g.img(l) for l in nz[n]])
    return d, g, nz, r

def coords(g, basis, v):
    return solve([g.img(l) for l in basis], v)

def struct(g, B, i, j):
    """dict (lam,mu) -> coefficient vector of s_lam s_mu in B[i+j]."""
    out = {}
    for lam in B[i]:
        for mu in B[j]:
            if i == j and mu < lam: continue
            v = g.img_of(lr(lam, mu), i + j)
            out[(lam, mu)] = coords(g, B[i + j], v)
    return out

def negatives(g, B, d, nmax=None):
    neg = []
    for n in range(2, d + 1):
        for i in range(1, n // 2 + 1):
            for (lam, mu), c in struct(g, B, i, n - i).items():
                for nu, x in zip(B[n], c):
                    if x < 0: neg.append((lam, mu, nu, x))
    return neg

def bases_in_degree(g, nz, r, n):
    return [S for S in itertools.combinations(nz[n], r[n]) if rank([g.img(l) for l in S]) == r[n]]

def search(k, verbose=True):
    """Backtrack over degrees; constraints for degree n only involve B_<=n."""
    d, g, nz, r = setup(k)
    choices = {n: bases_in_degree(g, nz, r, n) for n in range(d + 1)}
    # precompute: for each degree n, the set of Schur-product images landing there depends on B_i,B_j.
    count = [0]; sols = []
    surv = {n: 0 for n in range(d + 1)}
    def ok(B, n):
        for i in range(1, n // 2 + 1):
            for c in struct(g, B, i, n - i).values():
                if any(x < 0 for x in c): return False
        return True
    def rec(n, B):
        if n > d:
            sols.append(dict(B)); return
        for S in choices[n]:
            B[n] = S
            if ok(B, n):
                surv[n] += 1
                rec(n + 1, B)
            del B[n]
    rec(0, {})
    return d, g, nz, r, choices, sols, surv

if __name__ == '__main__':
    for k in map(int, sys.argv[1:]):
        d, g, nz, r, choices, sols, surv = search(k)
        print('=' * 60)
        print('k =', k, 'Hilb', [r[n] for n in range(d + 1)])
        print(' #Schur bases per degree:', [len(choices[n]) for n in range(d + 1)])
        print(' total combos:', eval('*'.join(str(len(choices[n])) for n in range(d + 1))))
        print(' # fully POSITIVE Schur bases:', len(sols))
        for s in sols[:6]:
            print('   ', {n: s[n] for n in range(d + 1) if r[n] > 1 or len(choices[n]) > 1})
        # natural greedy bases: lex-largest first, lex-smallest first
        for name, order in [('lex-desc (row-heavy first)', lambda L: L), ('lex-asc (column-heavy first)', lambda L: L[::-1])]:
            B = {}
            for n in range(d + 1):
                S = []
                for l in order(list(nz[n])):
                    if rank([g.img(x) for x in S + [l]]) > len(S): S.append(l)
                B[n] = tuple(S)
            neg = negatives(g, B, d)
            print(f' greedy {name}: B =', {n: B[n] for n in range(d + 1) if r[n] > 1})
            print(f'   negatives ({len(neg)}):', [(a, b, c, str(x)) for a, b, c, x in neg[:8]])
