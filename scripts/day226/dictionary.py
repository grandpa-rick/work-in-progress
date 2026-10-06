# Dictionary test: Phi_a(P_rho; x,y) := <T_a P_rho, p_x p_y>  (Hall pairing)
#   claimed  = (1-t^x)(1-t^y) X^lam_{(x,y)}(t) / b_lam(t),  lam = rho + 1^a  (Day 223 Lemma 1.3: T_a P_rho = P_lam).
# LHS evaluated by Rick's Thm 2.5 (closed two-string formula) with g = P_rho in a variables.
# Mode S (symbolic t): G_A from P_rho = sum_nu X^rho_nu prod(1-t^nu_i)/(b_rho z_nu) p_nu.
# Mode N (numeric t, independent of X-engine on the g side): G_A(w) from P_rho by S_a-symmetrization at rational t,w, interpolated in w.
# grade: computed
import sys, itertools, sympy as sp
from fractions import Fraction as Fr
from green import t, parts, X, b, zee, mult, chi
w = sp.symbols('w')
def thm25(a, d, Gcoef, gpi, x, y, tt=t):
    """second form of Thm 2.5; Gcoef(A,k) = [w^k] G_A(w), gpi = g(pi_a)."""
    tot = gpi; corr = 0
    for A in range(1, a):
        B = a-A; s = 0
        for k in range(0, d+1):
            G = Gcoef(A, k)
            if G == 0: continue
            s += G*(tt**(A*(k-y+B)) if k >= y-B else tt**(B*(y-B-k)))
        tot += tt**(-A*B)*s
        corr += tt**(-A*B)*Gcoef(A, y-B)/((1-tt**A)*(1-tt**B))
    return (-1)**(a-1)*((1-tt**x)*(1-tt**y)/(1-tt**a)*tot-(1-tt**x)*(1-tt**y)*corr)
def P_in_p(rho):
    d = sum(rho)
    return {nu: X(rho, nu)*sp.prod([1-t**v for v in nu])/(b(rho)*zee(nu)) for nu in parts(d)} if d else {(): sp.Integer(1)}
def eval_p(coef, pts):
    return sum(c*sp.prod([sum(p**r for p in pts) for r in nu]) for nu, c in coef.items())
def phi_symbolic(rho, a, x, y):
    d = sum(rho); coef = P_in_p(rho)
    GA = {}
    for A in range(1, a):
        B = a-A
        GA[A] = sp.Poly(sp.expand(eval_p(coef, [t**i for i in range(A)]+[w*t**j for j in range(B)])), w)
    Gc = lambda A, k: GA[A].coeff_monomial(w**k) if 0 <= k <= d else 0
    gpi = eval_p(coef, [t**i for i in range(a)])
    return sp.cancel(thm25(a, d, Gc, gpi, x, y))
# ---- numeric, independent P_rho
def P_sym(rho, pts, tt):
    a = len(pts); lam = list(rho)+[0]*(a-len(rho))
    tot = Fr(0)
    for e in set(itertools.permutations(lam)):
        # e[v] = exponent carried by variable v (Macdonald III (2.2): sum over S_a/S_a^lam)
        term = Fr(1)
        for v in range(a): term *= pts[v]**e[v]
        for i in range(a):
            for j in range(a):
                if e[i] > e[j]: term *= (pts[i]-tt*pts[j])/(pts[i]-pts[j])
        tot += term
    return tot
def lagrange_coeffs(xs, ys):
    n = len(xs); M = sp.Matrix([[sp.Rational(xx)**k for k in range(n)] for xx in xs])
    return list(M.LUsolve(sp.Matrix([sp.Rational(yy) for yy in ys])))
def phi_numeric(rho, a, x, y, tt):
    d = sum(rho); GA = {}
    wv = [Fr(7, 3)+Fr(k, 5) for k in range(d+1)]
    for A in range(1, a):
        B = a-A
        ys = [P_sym(rho, [tt**i for i in range(A)]+[ww*tt**j for j in range(B)], tt) for ww in wv]
        GA[A] = [Fr(int(c.p), int(c.q)) for c in lagrange_coeffs(wv, ys)]
    Gc = lambda A, k: GA[A][k] if 0 <= k <= d else 0
    gpi = P_sym(rho, [tt**i for i in range(a)], tt)
    return thm25(a, d, Gc, gpi, x, y, tt)
if __name__ == '__main__':
    mode = sys.argv[1]; nmax = int(sys.argv[2])
    cnt = bad = 0
    for n in range(2, nmax+1):
        for lam in parts(n):
            a = len(lam); rho = tuple(p-1 for p in lam if p > 1)
            if a < 2:  # a=1: only the one-string term; Thm 2.5 needs A in 1..a-1 -> check separately
                pass
            for x in range(1, n):
                y = n-x
                target = (1-t**x)*(1-t**y)*X(lam, tuple(sorted((x, y), reverse=True)))/b(lam)
                if mode == 'S':
                    lhs = phi_symbolic(rho, a, x, y); ok = sp.cancel(lhs-target) == 0
                else:
                    ok = True
                    for tt in (Fr(1, 3), Fr(-2, 7), Fr(5, 2)):
                        lhs = phi_numeric(rho, a, x, y, tt)
                        if Fr(lhs) != Fr(str(sp.Rational(target.subs(t, sp.Rational(tt.numerator, tt.denominator))))):
                            ok = False
                cnt += 1
                if not ok: bad += 1; print('FAIL', lam, (x, y))
        print('n =', n, 'cumulative checked', cnt, 'bad', bad, flush=True)
