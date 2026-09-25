"""Expand P_{k-1}(u) (sign-normalised to positive leading coeff) in Newton basis N_j = prod_{i=1}^j (t^i u - 1)."""
import sympy as sp, pickle
q,t,u=sp.symbols('q t u')
def qi(n,x=t): return sum(x**i for i in range(n))
P1=q*(t*u-1)+qi(2)
P2=q**3*(t*u-1)*(t**2*u-1)+qi(3)*(q**2*(t*u-1)+q*(t-1)+1)
P3=-pickle.load(open('P3.pkl','rb'))
def newton(P,k):
    P=sp.Poly(sp.expand(P),u); coeffs={}
    rem=P.as_expr()
    for j in range(k-1,-1,-1):
        Nj=sp.prod([(t**i*u-1) for i in range(1,j+1)])
        lc=sp.Poly(rem,u).coeff_monomial(u**j)
        c=sp.cancel(lc/sp.Poly(Nj,u).LC())
        coeffs[j]=sp.factor(c); rem=sp.expand(rem-c*Nj)
    assert sp.simplify(rem)==0
    return coeffs
for k,P in ((2,P1),(3,P2),(4,P3)):
    print(f'k={k}:')
    for j,c in sorted(newton(P,k).items(),reverse=True): print(f'   N_{j}: {c}')
# alternative basis: powers of (u-1) with t-shifts? also try basis M_j = prod_{i=1}^{j}(t^{i} u - 1) but on u'=t u etc.
