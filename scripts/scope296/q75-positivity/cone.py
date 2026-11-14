#!/usr/bin/env python3
"""Cone of Schur images in each G_n: extreme rays (exact, Caratheodory enumeration).
If the cone is simplicial (#rays == r_n) in every degree, then the ray basis has ALL Schur images
>= 0 in it, hence all products s_lam s_mu (Schur-positive in Sym) have >= 0 structure constants.
Grade: computed."""
import sys, itertools
from engine import *

def in_cone(v, gens, g):
    for m in range(1, len(gens) + 1):
        for S in itertools.combinations(gens, m):
            rows = [g.img(l) for l in S]
            if rank(rows) < m: continue
            x = solve(rows, v)
            if x is not None and all(t >= 0 for t in x): return True
    return False

for k in map(int, sys.argv[1:]):
    d = 2 * k; R = ribbon((2,) * k); g = G(R, d)
    print('=' * 60); print('k =', k)
    for n in range(d + 1):
        nz = [l for l in partitions(n) if any(g.img(l))]
        r = rank([g.img(l) for l in nz])
        # group parallel images
        classes = []
        for l in nz:
            for c in classes:
                if rank([g.img(c[0]), g.img(l)]) == 1: c.append(l); break
            else: classes.append([l])
        reps = [c[0] for c in classes]
        rays = [c for c in classes if not in_cone(g.img(c[0]), [x for x in reps if x != c[0]], g)]
        print(f' n={n} r={r} #directions={len(classes)} #extreme rays={len(rays)} simplicial={len(rays)==r}')
        if r >= 2: print('    rays:', rays)
