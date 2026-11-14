#!/usr/bin/env python3
"""k=6 chain-uniform certificate: products into degree 8 coming from bases common to ALL 128 positive
partial chains (B_2..B_6 parts shared by every chain).  Show: no Schur-subset B_8 makes all of them >= 0.
Grade: computed."""
import itertools
from engine import *
from positivity import setup, bases_in_degree
d, g, nz, r = setup(6)
# common elements across the 128 chains (read from obstruction_k6.txt structure: recompute quickly)
exec(open('obstruction_k6.py').read().split("print('#positive chains")[0].split("rec(0, {})")[0] + "rec(0, {})")
common = {n: set.intersection(*[set(c[n]) for c in chains]) for n in range(1, 8)}
print('elements common to all', len(chains), 'chains:', {n: sorted(common[n], reverse=True) for n in range(1, 8)})
prods = [(l, m) for i in range(1, 5) for l in common[i] for m in common[8 - i] if i < 8 - i or l <= m]
print('#common products into degree 8:', len(prods))
fail_all = True; hist = {}
for S in choices8 if False else bases_in_degree(g, nz, r, 8):
    w = next(((l, m) for l, m in prods if min(coef(S, 8, l, m)) < 0), None)
    if w is None: fail_all = False; print('   B_8 surviving common products:', S)
    else: hist[w] = hist.get(w, 0) + 1
print('every Schur B_8 killed by a common product:', fail_all)
print('witness histogram:', hist)
# smallest certificate set: greedy cover
left = set(bases_in_degree(g, nz, r, 8)); cover = []
while left:
    best = max(prods, key=lambda p: sum(1 for S in left if min(coef(S, 8, *p)) < 0))
    killed = {S for S in left if min(coef(S, 8, *best)) < 0}
    if not killed: break
    cover.append(best); left -= killed
print('greedy cover of all', len(bases_in_degree(g, nz, r, 8)), 'B_8 by products:', cover, 'uncovered:', len(left))
for p in cover:
    print('  ', p, ' = in Sym:', lr(*p))
