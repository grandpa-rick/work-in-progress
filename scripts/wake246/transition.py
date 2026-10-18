"""T := (G * phi0(w^-1))^{-1} phiinf,  G(M) = prod_{m in M} w_m^{-4}.  Also checks deg_q Ã_{M,I} = dim - codim I."""
import sys, pickle, sympy as sp
sys.path.insert(0, '/home/agent/projects/work-in-progress/scripts/wake246')
from common import *
N = int(sys.argv[1]); S = subsets(N); w = sp.symbols('w1:%d' % (N+1)); dim = N*(N+1)//2
A = load(N); inv = {x: 1/x for x in w}
bad = [(M, I, qdeg(A[M, I])) for M in S for I in S if A[M, I] != 0 and qdeg(A[M, I]) != dim - sum(I)]
print('entries with deg_q != dim-codim I:', bad)
Pinf = {(M, I): (lead_inf(A[M, I], dim - sum(I)) if A[M, I] != 0 else 0) for M in S for I in S}
P0b = {(M, I): sp.cancel(A[M, I].subs(q, 0).subs(inv, simultaneous=True) * sp.prod([w[m-1]**-4 for m in M])) for M in S for I in S}
nonmono = 0
for J in S:
    c = tri_solve(P0b, {M: Pinf[M, J] for M in S}, S, simp=lambda e: sp.factor(sp.cancel(e)))
    nz = {I: c[I] for I in c if c[I] != 0}
    for v in nz.values():
        n, d = sp.fraction(v)
        if len(sp.Poly(n, *w).terms()) > 1 or len(sp.Poly(d, *w).terms()) > 1: nonmono += 1
    print('col J=%s:' % (J,), nz, flush=True)
print('non-monomial T entries:', nonmono)
pickle.dump((Pinf, P0b), open('phis_N%d.pkl' % N, 'wb'))
