# Wake 251 reciprocity / period-collapse compute for F_kappa(N) (shifted diagram, diag entries even).
import sys, itertools
sys.path.insert(0, '/home/agent/projects/work-in-progress/scripts/wake251/orderpoly')
from orderpoly import ideals, shifted, boxless, longest_chain
import sympy as sp
x = sp.Symbol('x')

def setup(kap):
    P = shifted(sorted(kap, reverse=True)); J = ideals(P, boxless)
    diag = lambda b: b[0] == b[1]
    def antichain(S): return all(not boxless(a, b) for a in S for b in S)
    # transitions I -> K (I subset K); D = K-I; record (has_diag, is_antichain)
    tr = []
    for I in J:
        for K in J:
            if I <= K:
                D = K - I
                tr.append((I, K, any(diag(b) for b in D), antichain(D)))
    return P, J, tr

def count(P, J, tr, lo, hi, strict, diagpar=0):
    # #maps g:P->{lo..hi}, order-preserving (strict if strict), g(diag) = diagpar mod 2 (None = no cond).
    # level-set DP: after processing value v, ideal = {g<=v}.
    vec = {frozenset(): 1}
    for v in range(lo, hi + 1):
        nv = {}
        ok_d = diagpar is None or v % 2 == diagpar
        for I, K, hd, ac in tr:
            c = vec.get(I)
            if not c: continue
            if hd and not ok_d: continue
            if strict and not ac: continue
            nv[K] = nv.get(K, 0) + c
        vec = nv
    return vec.get(frozenset(P), 0)

def F(S, N, par=0): return count(*S, 0, N, False, par)
def Fint(S, M, par=0, lo=1): return count(*S, lo, M - 1, True, par)

def fitpoly(pts):
    return sp.expand(sp.interpolate(pts, x))
