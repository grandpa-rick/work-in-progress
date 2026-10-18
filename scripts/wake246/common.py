"""Shared helpers (Wake 246).  phi0_I(M) = Ã_{M,I}|_{q=0};  phiinf_I(M) = lim_{q->inf} q^{-(dim-codim I)} Ã_{M,I}.
tri_solve: coordinates in a basis Phi[(M,I)] triangular for Motegi's order (Phi[M,I]=0 unless M>=I)."""
import sys, pickle, sympy as sp
sys.path.insert(0, '/home/agent/projects/work-in-progress/scripts/day242/q484')
from engine import q, subsets, geq
p = sp.Symbol('p')
def load(N):
    return pickle.load(open('/home/agent/projects/work-in-progress/scripts/wake245/A_full_N%d.pkl' % N, 'rb'))
def order(S): return sorted(S, key=lambda T: (len(T), sum(T)))
def qdeg(a):
    n, d = sp.fraction(sp.cancel(a)); return sp.Poly(n, q).degree() - sp.Poly(d, q).degree()
def lead_inf(a, k):
    return sp.cancel(sp.cancel(a.subs(q, 1/p) * p**k).subs(p, 0))
def tri_solve(Phi, F, S, simp=sp.cancel):
    """F: dict M->value.  returns c: dict I->coef with F(M) = sum_I c_I Phi[M,I]."""
    c = {}
    for M in order(S):
        v = F[M] - sum(c[I] * Phi[M, I] for I in c if Phi[M, I] != 0 and c[I] != 0)
        c[M] = simp(v / Phi[M, M])
    return c
