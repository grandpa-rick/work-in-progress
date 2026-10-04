"""Day 221: verify Clio's two defects (review 2026-10-03) directly from the E_k operator definition.
E_k F = sum_{|A|=k} prod_{i in A, j notin A}(x_i - t x_j)/(x_i - x_j) X_A F(X_{A^c}, s X_A).
Defect 1: E_2(e_1^2) has [e_{2,2}] = 0.  Defect 2: e*_{(1,1)} = E_1E_1(1) has e_2 coeff (1-s)(1+t)."""
import itertools, sys
from sympy import Integer, symbols, Poly, cancel, factor, together, expand, Matrix, Rational
N = int(sys.argv[1]) if len(sys.argv) > 1 else 4
s, t = symbols('s t'); x = symbols(f'x0:{N}')
def Ek(k, F):
    tot = 0
    for A in itertools.combinations(range(N), k):
        B = [j for j in range(N) if j not in A]
        c = 1
        for i in A:
            for j in B: c *= (x[i]-t*x[j])/(x[i]-x[j])
        XA = 1
        for i in A: XA *= x[i]
        tot += c*XA*F.subs({x[i]: s*x[i] for i in A}, simultaneous=True)
    return expand(cancel(together(tot)))
def esym(k): return sum((lambda p: p)(eval('*'.join(f'x[{i}]' for i in A))) for A in itertools.combinations(range(N), k)) if k else 1
def parts(n, m=None):
    if m is None: m = n
    if n == 0: yield (); return
    for k in range(min(n, m), 0, -1):
        for p in parts(n-k, k): yield (k,)+p
def e_expand(F, n):
    P = [p for p in parts(n) if len(p) <= N] ; lam = [p for p in parts(n) if p[0] <= N]
    # coefficient of x^mu (mu partition) in e_lam and F
    def mono(G, mu): return Poly(G, *x).coeff_monomial(tuple(list(mu)+[0]*(N-len(mu))))
    eL = {l: expand(eval('*'.join(f'esym({a})' for a in l))) for l in lam}
    M = Matrix([[mono(eL[l], mu) for l in lam] for mu in P])
    v = Matrix([mono(F, mu) for mu in P])
    sol = M.solve_least_squares(v) if M.shape[0] != M.shape[1] else M.LUsolve(v)
    check = expand(F - sum(sol[i]*eL[l] for i, l in enumerate(lam)))
    assert check == 0, 'e-expansion residual nonzero'
    return {l: factor(sol[i]) for i, l in enumerate(lam)}
e1 = esym(1)
print(f'N={N}')
r = e_expand(Ek(2, expand(e1**2)), 4)
print('E_2(e_1^2) =', r)
print('[e_{2,2}] =', r[(2, 2)])
r2 = e_expand(Ek(1, Ek(1, Integer(1))), 2)
print('e*_{(1,1)} = E_1E_1(1) =', r2)
print('at t=-1:', {k: factor(v.subs(t, -1)) for k, v in r2.items()})
r3 = e_expand(Ek(1, Ek(1, Ek(1, Integer(1)))), 3)
print('e*_{(1,1,1)} =', r3)
print('at t=-1:', {k: factor(v.subs(t, -1)) for k, v in r3.items()})
