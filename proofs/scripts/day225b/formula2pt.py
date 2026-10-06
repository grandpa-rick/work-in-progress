# Conjectured closed two-point formula (Day 225). grade: computed
import sympy as sp
from sym2pt import Phi
t, w = sp.symbols('t w')
def gpoly(nu, pts):
    return sp.prod([sum(p**r for p in pts) for r in nu]) if nu else sp.Integer(1)
def Phi_formula(a, nu, x, y):
    d = sum(nu)
    tot = gpoly(nu, [t**i for i in range(a)])
    corr = 0
    for A in range(1, a):
        B = a - A
        G = sp.Poly(sp.expand(gpoly(nu, [t**i for i in range(A)] + [w*t**j for j in range(B)])), w)
        Gk = lambda k: G.coeff_monomial(w**k) if k >= 0 else 0
        s = 0
        for k in range(0, d+1):
            if k >= y-B: s += Gk(k)*t**(A*(k-y+B))
            else: s += Gk(k)*t**(B*(y-B-k))
        tot += t**(-A*B)*s
        corr += t**(-A*B)*Gk(y-B)/((1-t**A)*(1-t**B))
    val = (1-t**x)*(1-t**y)/(1-t**a)*tot - (1-t**x)*(1-t**y)*corr
    return sp.factor(sp.cancel((-1)**(a-1)*val))
if __name__ == '__main__':
    def parts(n, m=None):
        if m is None: m = n
        if n == 0: yield (); return
        for k in range(min(n, m), 0, -1):
            for r in parts(n-k, k): yield (k,)+r
    bad = 0; cnt = 0
    for a in (1, 2, 3, 4):
        for d in range(0, 5 if a < 4 else 4):
            for nu in parts(d):
                n = a + d
                for x in range(1, n):
                    y = n - x
                    f = Phi_formula(a, nu, x, y); e = Phi(a, nu, x, y)
                    cnt += 1
                    if sp.simplify(f - e) != 0: bad += 1; print('BAD', a, nu, x, y, f, e)
    print('checked', cnt, 'bad', bad)
