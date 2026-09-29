"""Day 208 (Sep 29): structure fit of Gamma_k(z,w) = sum_{(lam,i,j)} R * e_lam E(t^i z)E(t^j w)  (gf_recursion.py output).
Tests: (S1) R^{(k)}_{(b),i,j} = s^{2b} t^{-b(i+j)} R^{(k-b)}_{(),i,j};
(S2) R^{(n)}_{(),i,j} = K_ij(z,w) * sum_{n1+n2=n} omega * N^{(n1)}_i z^{-n1} N^{(n2)}_j w^{-n2},  N^{(n)}_j = t^{-nj} c(n,j) (one-column),
K_ij = prod_{p<i}(t^p z - s w)/(t^p z - t^j w) * prod_{r<j}(s z - t^r w)/(t^i z - t^r w). Reports omega (should be monomials)."""
import pickle, sympy as sp
from gf_recursion import s, t, z, w
res = pickle.load(open('gf_Gamma_upto3.pkl', 'rb'))
poch = lambda x, N: sp.Mul(*[1 - x*t**i for i in range(N)])
al = lambda J: 0 if J < 0 else sp.Mul(*[s - t**i for i in range(1, J+1)]) / poch(t, J)
def c(n, j):
    if j < 0 or j > n: return 0
    return poch(s, n-j)/poch(t, n-j)*(al(j) - s*t**(n-j)*al(j-1))
N = lambda n, j: t**(-n*j)*c(n, j)
K = lambda i, j: sp.Mul(*[(t**p*z - s*w)/(t**p*z - t**j*w) for p in range(i)]) * sp.Mul(*[(s*z - t**r*w)/(t**i*z - t**r*w) for r in range(j)])
ok1 = True
for k in (2, 3):
    for (lam, i, j), R in res[k].items():
        if lam == (): continue
        b = lam[0]; assert len(lam) == 1
        R0 = res[k-b].get(((), i, j), 0)
        ok1 &= sp.cancel(R - s**(2*b)*t**(-b*(i+j))*R0) == 0
print('(S1) prefix rule R_(b),i,j = s^{2b} t^{-b(i+j)} R^{(k-b)}_(),i,j :', ok1)
Z, W = sp.symbols('Z W')  # Z=1/z, W=1/w
for n in (1, 2, 3):
    print(f'--- n={n}')
    for i in range(n+1):
        for j in range(n+1-i):
            R = res[n].get(((), i, j), 0)
            A = sp.cancel(R / K(i, j))
            A = sp.expand(sp.cancel(A.subs({z: 1/Z, w: 1/W})))
            num, den = sp.fraction(sp.together(A))
            P = sp.Poly(sp.expand(A), Z, W) if sp.Poly(den, Z, W).degree() == 0 else None
            if P is None: print('  (i,j)=', (i, j), 'NOT Laurent-poly after /K:', sp.factor(A)); continue
            out = []
            for (a1, a2), cf in P.terms():
                base = N(a1, i)*N(a2, j)
                om = sp.factor(sp.cancel(cf/base)) if base != 0 else ('base0', sp.factor(cf))
                out.append(((a1, a2), om))
            print('  (i,j)=', (i, j), ' omega[(n1,n2)] =', out)
