"""(a) Verify the telescoped GF against the direct AHA pipeline e_3(Y).e_r (flint, exact poly -> random pts),
m=6,7, r=0..5; q_n computed directly from X (independent of the q->e conversion). Also check negative z-powers vanish."""
import sys, pickle, random, time
import sympy as sp
from fractions import Fraction as Fr
sys.path.insert(0, '/home/agent/projects/scripts/day207'); sys.path.insert(0, '/home/agent/projects/proofs/scripts/day205')
from k3_kernel_test import ekY_er, evalp, rand_pts, e_vals, qn_list
from k3_fast_pipeline import AHA
s, t, z = sp.symbols('s t z')
C, _ = pickle.load(open('/home/agent/projects/scripts/day207/k3_gf_coeffs.pkl', 'rb'))
Cz = {j: sp.expand(C[j] * z**3) for j in C}   # polynomial in z (deg <= 6)
qsyms = [sp.Symbol(f'q{n}') for n in (1, 2, 3)]
ok = True
for j in Cz: assert sp.Poly(Cz[j], z).degree() <= 6
# negative powers of z vanish identically: sum_j [z^{r}] for r<0  (use e's generic via direct eval below)
for m in (6, 7):
    A_ = AHA(m)
    for r in range(0, 6):
        t0 = time.time()
        P = ekY_er(A_, 3, r)
        pts = rand_pts(m, 3, 500 + 10*m + r)
        for X, sv, tv in pts:
            qn = qn_list(X, tv, 3)
            sub = {s: sv, t: tv, qsyms[0]: qn[1], qsyms[1]: qn[2], qsyms[2]: qn[3]}
            gf = Fr(0)
            for j, Cj in Cz.items():
                cj = sp.Poly(Cj.subs(sub), z)
                for (b,), cb in cj.terms():   # term cb z^{b-3} E(t^j z): [z^r] -> cb t^{j(r-b+3)} e_{r-b+3}
                    k = r - b + 3
                    if k < 0: continue
                    gf += Fr(str(cb)) * tv**(j*k) * e_vals(X, k)
            fact = (1 - tv)**3 * (1 + tv) * (1 + tv + tv**2)
            lhs = evalp(A_, P, X, sv, tv) / tv**3
            if lhs * fact != gf: ok = False; print('FAIL', m, r)
        # negative r: [z^r], r=-1,-2,-3 must vanish
        X, sv, tv = pts[0]; qn = qn_list(X, tv, 3)
        sub = {s: sv, t: tv, qsyms[0]: qn[1], qsyms[1]: qn[2], qsyms[2]: qn[3]}
        for rr in (-1, -2, -3):
            g = Fr(0)
            for j, Cj in Cz.items():
                for (b,), cb in sp.Poly(Cj.subs(sub), z).terms():
                    k = rr - b + 3
                    if k >= 0: g += Fr(str(cb)) * tv**(j*k) * e_vals(X, k)
            if g != 0: ok = False; print('neg power nonzero', rr)
        print(f'm={m} r={r}: {"OK" if ok else "FAIL"}  [{time.time()-t0:.1f}s]'); sys.stdout.flush()
print('ALL', 'OK' if ok else 'FAIL')
