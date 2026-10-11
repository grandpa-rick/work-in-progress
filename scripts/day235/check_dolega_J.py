# Day 235: printed FPSAC Sec. 4 separator claim, re-implemented from the printed sentence:
#  "Dolega's (+)-cumulant of columns (lowest (q-1)-coefficient of [g_n] Htilde_{lam'}, g_k = Htilde_{1^k}) gives,
#   by computation, J = 3/2 = J_star(0), whereas J_star(1) = 1", with J := L(1^4) L(2,2) / L(2,1,1)^2 and
#   J_star = (t^2+1)(t^3+3t^2+6t+6) / ((t+1)(t^2+t+2)^2) (printed, from Thm 4.1).
# Htilde_mu built from scratch by the Garsia-Haiman triangularity characterisation:
#   Htilde_mu[X(1-q)] in span{s_lam : lam >= mu}, Htilde_mu[X(1-t)] in span{s_lam : lam >= mu'}, <Htilde_mu, s_n> = 1.
# NEGATIVE CONTROL (changes the claim's input, not a symmetry of it): rows instead of columns, Htilde_lam in place
#   of Htilde_{lam'}; the t-independence / value 3/2 must then fail.
import sympy as sp, itertools
q, t = sp.symbols('q t')
NMAX = 4
p = sp.symbols('p1:%d' % (NMAX + 1))

def partitions(n, mx=None):
    if mx is None: mx = n
    if n == 0: yield (); return
    for k in range(min(n, mx), 0, -1):
        for r in partitions(n - k, k): yield (k,) + r
def conj(l): return tuple(sum(1 for x in l if x > i) for i in range(l[0])) if l else ()
def dom(mu, lam):
    a = b = 0
    for i in range(max(len(mu), len(lam))):
        a += mu[i] if i < len(mu) else 0; b += lam[i] if i < len(lam) else 0
        if a < b: return False
    return True
def z(rho):
    from collections import Counter
    r = 1
    for k, m in Counter(rho).items(): r *= k ** m * sp.factorial(m)
    return r
def h(n):
    if n < 0: return 0
    if n == 0: return 1
    return sum(sp.prod([p[i - 1] for i in rho]) / z(rho) for rho in partitions(n))
def schur(lam):
    l = len(lam)
    return sp.expand(sp.Matrix(l, l, lambda i, j: h(lam[i] - i + j)).det())
def pleth(F, a):  # F[X*a] where a is a 'monomial-like' factor: p_k -> p_k * a(k)
    return sp.expand(F.subs({p[k - 1]: p[k - 1] * a(k) for k in range(1, NMAX + 1)}, simultaneous=True))
def pmon(rho): return sp.prod([p[i - 1] for i in rho])

def schur_coeffs(F, n):
    parts = list(partitions(n)); S = [schur(l) for l in parts]
    cs = sp.symbols('c0:%d' % len(parts))
    eq = sp.Poly(sp.expand(F - sum(c * s for c, s in zip(cs, S))), *p).coeffs()
    sol = sp.solve(eq, cs, dict=True)[0]
    return {l: sp.simplify(sol[c]) for l, c in zip(parts, cs)}

def Htilde(mu):
    n = sum(mu); parts = list(partitions(n)); S = {l: schur(l) for l in parts}
    a = sp.symbols('a0:%d' % len(parts)); H = sum(ai * S[l] for ai, l in zip(a, parts))
    eqs = [a[parts.index((n,))] - 1]
    for plet, cond in ((lambda k: 1 - q ** k, mu), (lambda k: 1 - t ** k, conj(mu))):
        C = schur_coeffs(pleth(H, plet), n)
        eqs += [C[l] for l in parts if not dom(l, cond)]
    sol = sp.solve(eqs, a, dict=True); assert len(sol) == 1
    return sp.expand(sp.simplify(H.subs(sol[0])))

H = {}
for n in range(1, NMAX + 1):
    for mu in partitions(n): H[mu] = Htilde(mu)
# instrument checks: Htilde_2 = s2 + q s11, Htilde_11 = s2 + t s11; Htilde_mu(q,t) = Htilde_mu'(t,q)
print('instr Htilde_2 :', sp.simplify(H[(2,)] - schur((2,)) - q * schur((1, 1))) == 0)
print('instr Htilde_11:', sp.simplify(H[(1, 1)] - schur((2,)) - t * schur((1, 1))) == 0)
print('instr q<->t duality n<=4:', all(sp.expand(H[mu] - H[conj(mu)].subs({q: t, t: q}, simultaneous=True)) == 0 for mu in H))

def lead_cumulant(lam, rows=False):
    n = sum(lam); target = H[lam] if rows else H[conj(lam)]
    g = {k: H[(1,) * k] for k in range(1, n + 1)}
    parts = list(partitions(n)); cs = sp.symbols('d0:%d' % len(parts))
    G = sum(c * sp.prod([g[k] for k in nu]) for c, nu in zip(cs, parts))
    eq = sp.Poly(sp.expand(target - G), *p).coeffs()
    sol = sp.solve(eq, cs, dict=True); assert len(sol) == 1
    cn = sp.factor(sol[0][cs[parts.index((n,))]])           # [g_n]
    e = sp.symbols('e'); ser = sp.expand(sp.cancel(cn.subs(q, 1 + e)))
    num, den = sp.fraction(sp.cancel(cn.subs(q, 1 + e)))
    num = sp.Poly(num, e); den = sp.Poly(den, e)
    vn = min(m[0] for m in num.monoms()); vd = min(m[0] for m in den.monoms())
    lo = sp.factor(num.coeff_monomial(e ** vn) / den.coeff_monomial(e ** vd))
    return vn - vd, lo

Jstar = (t**2 + 1) * (t**3 + 3*t**2 + 6*t + 6) / ((t + 1) * (t**2 + t + 2)**2)
for label, rows in (('printed (columns, Htilde_{lam\'})', False), ('NEG CONTROL (rows, Htilde_lam)', True)):
    L = {lam: lead_cumulant(lam, rows) for lam in [(1, 1, 1, 1), (2, 2), (2, 1, 1)]}
    J = sp.factor(sp.cancel(L[(1, 1, 1, 1)][1] * L[(2, 2)][1] / L[(2, 1, 1)][1] ** 2))
    print(label, ': leads', {k: v for k, v in L.items()})
    print('   J =', J, '| t-independent:', sp.diff(J, t) == 0, '| J == 3/2:', sp.simplify(J - sp.Rational(3, 2)) == 0)
print('J_star(0) =', Jstar.subs(t, 0), ' J_star(1) =', Jstar.subs(t, 1), ' J_star at t=2:', Jstar.subs(t, 2))
