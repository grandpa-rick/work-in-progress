"""Day 207 Test 3: k=3 generating-function derivation of t^{-3} e_3(Y).e_r by iterated Lemma 2.

Chain (argued as in Day 206b, k=3 kernel identity computed in k3_kernel_test.py):
  t^{-3} e_3(Y) e_r = sum_{|A|=3} Phi^(A) prod^x_A                       (K_T, computed m<=7)
  sum_{ordered (a,b,c)} a_ab a_ac a_bc = [3]_t!  (Poincare/Macdonald)      => divide by [3]!
  (1-t)^3 sum_ordered Phi(x_a,x_b,x_c) a_ab a_ac a_bc prod^x = Omega_x Omega_y Omega_w [Phi]
     Omega_w: ambient P(v)=Q(v)h(x,v)h(y,v); Omega_y: Q(v)h(x,v); Omega_x: Q(v); h(u,v)=(1-uv)/(1-tuv)
     Omega[u^n] = [v^n]P(v), n>=1  (Lemma 2 in the complement variables)
  Phi(x,y,w) = [z^r] E(z) phi(x)phi(y)phi(w),   phi(u) = u(1+s u z)/(1+u z).
Omega on rational f(u) with f(0)=0: partial fractions; A/(1-cu) -> A (P(c)-1); u^n -> P_n.
Q(c) kept as symbols Qc[c]; q_n symbols. At the end multiply by E(z) and try to telescope
  E(z) Q(-z) Q(-tz) ... Q(-t^{j-1} z) = E(t^j z).
"""
import sympy as sp, sys, itertools
s, t, z, v = sp.symbols('s t z v')
x, y, w = sp.symbols('x y w')
NQ = 10
qs = [sp.Integer(1)] + [sp.Symbol(f'q{n}') for n in range(1, NQ)]
Qsym = {}   # c -> symbol Q(c)
def Qof(c):
    c = sp.factor(c)
    if c not in Qsym:
        Qsym[c] = sp.Symbol(f'Q[{c}]')
    return Qsym[c]

def h(u, vv): return (1 - u * vv) / (1 - t * u * vv)

def consts():
    return list(Qsym.values()) + qs[1:]

def Omega(f, u, Aamb):
    """Aamb(v): explicit rational function in v (and other head vars); P = Q(v) Aamb(v)."""
    f = sp.expand(f)
    C = [c for c in consts() if f.has(c)]
    poly = sp.Poly(f, *C) if C else None
    terms = poly.terms() if C else [((), f)]
    out = 0
    for mon, coef in terms:
        cmon = sp.Mul(*[c**e for c, e in zip(C, mon)]) if C else 1
        g = sp.apart(sp.cancel(coef), u)
        res = 0; const0 = 0
        for term in sp.Add.make_args(g):
            num, den = sp.fraction(sp.factor(term))
            dp = sp.Poly(den, u)
            if dp.degree() == 0:
                pp = sp.Poly(sp.expand(term), u)
                for (n,), pn in pp.terms():
                    if n == 0: const0 += pn; continue
                    # P_n = sum_k A_k q_{n-k}
                    ser = sp.series(Aamb, v, 0, n + 1).removeO()
                    Pn = sum(ser.coeff(v, k) * qs[n - k] for k in range(n + 1))
                    res += pn * Pn
            else:
                assert dp.degree() == 1, ('higher-order pole', term)
                rho = sp.solve(den, u)[0]
                Acoef = sp.cancel(-(term * (u - rho)).subs(u, rho) / rho) if False else None
                B = sp.cancel(term * (u - rho))  # const in u
                B = sp.cancel(B)
                assert not B.has(u), term
                c = 1 / rho
                Ai = sp.cancel(-B * c)          # B/(u-rho) = Ai/(1-c u)
                const0 += Ai
                res += Ai * (Qof(c) * Aamb.subs(v, c) - 1)
        assert sp.simplify(const0) == 0, ('f(0) != 0', const0)
        out += cmon * res
    return sp.expand(out)

def phi(u): return u * (1 + s * u * z) / (1 + u * z)

if __name__ == '__main__':
    I3 = Omega(phi(w), w, qs and (sp.Symbol('dummy')*0 + 1) * 0 + h(x, v) * h(y, v))
    print('stage w done', len(sp.Add.make_args(I3)), 'terms'); sys.stdout.flush()
    I2 = Omega(sp.expand(phi(y) * I3), y, h(x, v))
    print('stage y done', len(sp.Add.make_args(I2))); sys.stdout.flush()
    I1 = Omega(sp.expand(phi(x) * I2), x, sp.Integer(1))
    print('stage x done', len(sp.Add.make_args(I1))); sys.stdout.flush()
    print('Q symbols:', Qsym)
    import pickle
    pickle.dump((I1, Qsym), open('/home/agent/projects/scripts/day207/k3_gf_I1.pkl', 'wb'))
    # group by monomials in Q symbols
    QC = list(Qsym.values())
    P = sp.Poly(I1, *QC)
    for mon, coef in P.terms():
        print(mon, sp.factor(coef))
