#!/usr/bin/env python3
"""(a)/(b): nonzero Schur images in G_(2^k), ranks, number of Schur-subset bases per degree."""
import sys, itertools
from engine import *
for k in map(int, sys.argv[1:]):
    d = 2 * k; R = ribbon((2,) * k); g = G(R, d)
    print('=' * 60); print('k =', k, ' R =', R)
    hil = []
    for i in range(d + 1):
        P = partitions(i)
        nz = [l for l in P if any(g.img(l))]
        r = rank([g.img(l) for l in nz]) if nz else 0
        hil.append(r)
        nb = sum(1 for S in itertools.combinations(nz, r) if rank([g.img(l) for l in S]) == r) if len(nz) <= 14 else 'skip'
        zero = [l for l in P if l not in nz]
        print(f' i={i} r_i={r} #nonzero s={len(nz)}/{len(P)} #Schur bases={nb}')
        print('    nonzero:', nz)
        if len(zero) <= 12: print('    zero   :', zero)
    print(' Hilb =', hil)
