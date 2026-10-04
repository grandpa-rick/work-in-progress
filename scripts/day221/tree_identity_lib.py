# Day 221 PROVE: check (i) increasing-tree identity with GENERIC edge weights x_ij;
# (ii) history sum with Thm W weights == telescoped tree formula == Conj G, symbolic t.
import itertools, sympy as sp, random
from fractions import Fraction as Fr
def inc_trees(l):  # parent arrays: par[c] < c for c=1..l-1 (0-indexed, root 0)
    for par in itertools.product(*[range(c) for c in range(1,l)]):
        yield (None,)+par
def subtree(par,l):
    S=[{i} for i in range(l)]
    for c in range(l-1,0,-1): S[par[c]] |= S[c]
    return S
def conn_graphs(l):
    E=list(itertools.combinations(range(l),2))
    for m in range(1<<len(E)):
        H=[E[i] for i in range(len(E)) if m>>i&1]
        p=list(range(l))
        def f(a):
            while p[a]!=a: a=p[a]
            return a
        for a,b in H: p[f(a)]=f(b)
        if len({f(a) for a in range(l)})==1: yield H
