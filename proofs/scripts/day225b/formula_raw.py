# Raw form of the two-point theorem (as stated in the Day 225 proof file, Thm 2.5). grade: computed
import sympy as sp
from sym2pt import Phi
from formula2pt import gpoly, t, w
def Phi_raw(a, nu, x, y):
    n = a + sum(nu); tot = 0
    for A in range(1, a):
        B = a - A
        G = sp.Poly(sp.expand(gpoly(nu, [t**i for i in range(A)] + [w*t**j for j in range(B)])), w)
        Gk = lambda k: G.coeff_monomial(w**k) if k >= 0 else 0
        tot += t**(-A*B)*(Gk(y-B)/((1-t**A)*(1-t**B)) + sum((t**(-A*m)-t**(B*m))*Gk(y-B-m) for m in range(1, y-B+1))/(1-t**a))
    tot -= gpoly(nu, [t**i for i in range(a)])/(1-t**a)*sum(t**(-j*y) for j in range(a))
    return sp.factor(sp.cancel((-1)**a*(1-t**x)*(1-t**y)*tot))
if __name__ == '__main__':
    def parts(n, m=None):
        if m is None: m = n
        if n == 0: yield (); return
        for k in range(min(n, m), 0, -1):
            for r in parts(n-k, k): yield (k,)+r
    cnt = bad = 0
    for a in (1, 2, 3, 4, 5):
        for d in range(0, {1:6,2:5,3:5,4:4,5:3}[a]):
            for nu in parts(d):
                n = a + d
                for x in range(1, n):
                    cnt += 1
                    if sp.simplify(Phi_raw(a, nu, x, n-x) - Phi(a, nu, x, n-x)) != 0: bad += 1; print('BAD', a, nu, x)
    print('raw form checked', cnt, 'bad', bad)
    print('example a=2,p2,(2,2):', Phi_raw(2, (2,), 2, 2))
