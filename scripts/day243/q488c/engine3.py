"""Fast A0, A1 (type C), unpadded n=|M| (Motegi stability). A = pref*tot/((n-r)! E), E(0)=1."""
import itertools, sympy as sp
from engine import summand, denom_pairs, geq
q=sp.Symbol('q')
def A01(N,M,I):
    s,r=len(M),len(I)
    if not geq(M,I): return sp.Integer(0),sp.Integer(0)
    w=sp.symbols('w1:%d'%(N+1)); n=s; u=[w[m-1] for m in M]
    t0=0;t1=0
    for perm in itertools.permutations(range(n)):
        us=[u[p] for p in perm]
        t=sp.expand(summand(us,w,I,N,n,1,0,q))
        P=sp.Poly(t,q); c=P.all_coeffs()[::-1]
        if all(x==0 for x in c): continue
        d=sp.prod([us[a]-us[b] for (a,b) in denom_pairs(n,r)])
        t0+=c[0]/d; t1+=(c[1] if len(c)>1 else 0)/d
    E=(1-q)**(-(n-r))*sp.prod([1-q*a/b for a in u for b in u])*sp.prod([1-q*a*b for a in u for b in w])
    E/=sp.prod([1-q*u[j]*u[k] for j in range(n) for k in range(j+1,n)])
    E/=sp.prod([1-q*w[j]*w[k] for j in range(N) for k in range(j,N)])
    E1=sp.diff(E,q).subs(q,0)
    pref=sp.prod([x**(-(N+n)) for x in u])/sp.factorial(n-r)
    a0=sp.factor(sp.cancel(pref*t0)); a1=sp.factor(sp.cancel(pref*(t1-t0*E1)))
    return a0,a1
