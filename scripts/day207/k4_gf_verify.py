"""(a) Verify the telescoped k=4 GF (M_{jb} form, q->e substituted) against the direct AHA pipeline
t^{-6} e_4(Y).e_r at m=8 (r=0..4) and m=9 (r=0..5), 3 random rational points each; e_b, e_n from X directly.
Also check [z^r] for r=-1..-4 vanishes. Additionally checks the conjectural general-k M_{jb} formula (gen_M_conjecture)."""
import sys, pickle, time
import sympy as sp
from fractions import Fraction as Fr
sys.path.insert(0, '/home/agent/projects/scripts/day207'); sys.path.insert(0, '/home/agent/projects/proofs/scripts/day205')
from k3_kernel_test import ekY_er, evalp, rand_pts, e_vals
from k3_fast_pipeline import AHA
K = 4
s, t = sp.symbols('s t')
C, M, _ = pickle.load(open('/home/agent/projects/scripts/day207/k4_gf_coeffs.pkl', 'rb'))
Mf = {key: sp.lambdify((s, t), val, 'sympy') for key, val in M.items()}
def gf_coef(X, sv, tv, r):
    tot = Fr(0)
    for (j, b), f in Mf.items():
        n = r + K - b
        if n < 0: continue
        tot += Fr(str(f(sp.Rational(sv.numerator, sv.denominator), sp.Rational(tv.numerator, tv.denominator)))) * e_vals(X, b) * tv**(j*n) * e_vals(X, n)
    return tot
ok = True
for m, rs in ((8, range(0, 5)), (9, range(0, 6))):
    A_ = AHA(m)
    for r in rs:
        t0 = time.time()
        P = ekY_er(A_, K, r)
        for X, sv, tv in rand_pts(m, 3, 900 + 10*m + r):
            fac = (1-tv)*(1-tv**2)*(1-tv**3)*(1-tv**4)
            if evalp(A_, P, X, sv, tv) / tv**6 * fac != gf_coef(X, sv, tv, r): ok = False; print('FAIL', m, r)
        X, sv, tv = rand_pts(m, 1, 7)[0]
        for rr in range(-K, 0):
            if gf_coef(X, sv, tv, rr) != 0: ok = False; print('neg power nonzero', m, rr)
        print(f'm={m} r={r}: {"OK" if ok else "FAIL"} [{time.time()-t0:.1f}s]'); sys.stdout.flush()
print('ALL', 'OK' if ok else 'FAIL')
