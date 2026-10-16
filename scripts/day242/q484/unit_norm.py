"""M-dependent renormalisation test: hat e_I = e_I / e_emptyset (makes the emptyset class the unit).
hat c^nu_{lam,mu} solves  sum_nu hat c e_nu(M) = e_lam(M) e_mu(M) / e_empty(M).
Also: denominators of the poles in the unnormalised (WC without E) census."""
import sys, sympy as sp
from engine import *
X, N, qarg, norm = sys.argv[1], int(sys.argv[2]), sys.argv[3], (sys.argv[4] == 'norm')
qq = q if qarg == 'q' else sp.Rational(qarg)
S = subsets(N)
A = {}
for M in S:
    for I in S:
        A[M, I], w = W(X, N, I, M, qq, normalise=norm)
assert all(A[M, ()] != 0 for M in S)
order = sorted(S, key=lambda T: (len(T), sum(T)))
lau = tot = unit_ok = 0; poles = set()
for l in S:
    for m in S:
        c = {}
        for nu in order:
            rhs = A[nu, l] * A[nu, m] / A[nu, ()] - sum(c[k] * A[nu, k] for k in c)
            c[nu] = sp.factor(sp.cancel(rhs / A[nu, nu]))
        if l == ():
            unit_ok += all(sp.cancel(c[nu] - (1 if nu == m else 0)) == 0 for nu in S)
        for nu in S:
            if geq(nu, l) and geq(nu, m):
                tot += 1
                L = is_laurent(c[nu], w); lau += L
                if not L:
                    poles |= {str(f) for f, e in sp.factor_list(sp.fraction(sp.cancel(c[nu]))[1])[1] if f.free_symbols & set(w)}
            else:
                assert c[nu] == 0
print(X, N, qarg, 'norm' if norm else 'unnorm', 'unit-normalised: Laurent', lau, '/', tot,
      '| hat e_empty is unit for', unit_ok, '/', len(S), 'mu | pole factors:', sorted(poles), flush=True)
