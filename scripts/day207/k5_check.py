"""k=5 cross-checks: (i) kernel e_5(Y).e_r = t^10 sigma^(5) pi^5 e_r exact (m=8,9); (ii) k=5 GF (engine M_{jb}) vs direct
t^{-10} e_5(Y).e_r at m=9, r=0..4, 2 random points each."""
import sys, pickle, time, sympy as sp
from fractions import Fraction as Fr
sys.path.insert(0, '/home/agent/projects/scripts/day207'); sys.path.insert(0, '/home/agent/projects/proofs/scripts/day205')
from k4_kernel_test import test
from k3_kernel_test import ekY_er, evalp, rand_pts, e_vals
from k3_fast_pipeline import AHA
for m, r in ((8, 2), (9, 1), (9, 3)):
    print('kernel k=5', m, r, test(m, r, k=5)); sys.stdout.flush()
K = 5; s, t = sp.symbols('s t')
_, M, _ = pickle.load(open('/home/agent/projects/scripts/day207/k5_gf_coeffs.pkl', 'rb'))
Mf = {key: sp.lambdify((s, t), val, 'sympy') for key, val in M.items()}
A_ = AHA(9); ok = True
for r in range(0, 5):
    P = ekY_er(A_, K, r)
    for X, sv, tv in rand_pts(9, 2, 50 + r):
        gf = sum(Fr(str(f(sp.Rational(sv.numerator, sv.denominator), sp.Rational(tv.numerator, tv.denominator)))) * e_vals(X, b) * tv**(j*(r+K-b)) * e_vals(X, r+K-b) for (j, b), f in Mf.items())
        fac = 1
        for i in range(1, K+1): fac *= (1 - tv**i)
        if evalp(A_, P, X, sv, tv) / tv**10 * fac != gf: ok = False; print('FAIL r', r)
    print('k=5 GF vs AHA m=9 r=%d:' % r, 'OK' if ok else 'FAIL'); sys.stdout.flush()
