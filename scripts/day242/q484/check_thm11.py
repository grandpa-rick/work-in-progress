"""Thm 11 (localisation formula via Motegi's all-cardinality inverse, M4) vs recursion; N=2, q symbolic.
Negative control: drop q^{|lambda(nu)|-|lambda(M)|} (must break)."""
import sys, sympy as sp
from engine import *
X = sys.argv[1] if len(sys.argv) > 1 else 'C'
N = 2
S = subsets(N)
A, C, w = structure_constants(X, N, q)
Ainv_q = {}  # A(w, q^{-1}) row-normalised
for M in S:
    for I in S:
        Ainv_q[M, I] = sp.cancel((A[M, I] / A[M, M]).subs(q, 1 / q))
comp = lambda T: tuple(j for j in range(1, N + 1) if j not in T)
ok = bad = tot = 0
# first: M4 itself, sum_I q^{..} A_{M,I}(q) A_{L^c,I^c}(q^{-1}) = delta
m4 = sum(sp.cancel(sum(q**(sum(I) - sum(L)) * A[M, I] / A[M, M] * Ainv_q[comp(L), comp(I)] for I in S)
                   - (1 if M == L else 0)) == 0 for M in S for L in S)
m4_plant = sum(sp.cancel(sum(A[M, I] / A[M, M] * Ainv_q[comp(L), comp(I)] for I in S)
                         - (1 if M == L else 0)) == 0 for M in S for L in S)
print(X, 'M4 biorthogonality', m4, '/', len(S)**2, '| planted (no q-power):', m4_plant, '/', len(S)**2, flush=True)
for l in S:
    for m in S:
        for nu in S:
            tot += 1
            f = sum(q**(sum(nu) - sum(M)) * Ainv_q[comp(M), comp(nu)] * A[M, l] * A[M, m] / A[M, M] for M in S)
            g = sum(Ainv_q[comp(M), comp(nu)] * A[M, l] * A[M, m] / A[M, M] for M in S)
            ok += sp.cancel(f - C[l, m][nu]) == 0
            bad += sp.cancel(g - C[l, m][nu]) == 0
print(X, 'Thm11 vs recursion', ok, '/', tot, '| planted (no q-power):', bad, '/', tot, flush=True)
