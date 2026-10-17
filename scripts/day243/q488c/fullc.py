"""Full-q check for a single triple: c^nu_{lam,mu}(q) on the lower ideal of nu, numeric w, w_b=1/(w_a(1-eps)).
Engine = Motegi summand (engine.summand), unpadded n=|M|."""
import itertools, random, sys, sympy as sp
from engine import summand, denom_pairs, geq
q,eps=sp.symbols('q eps')
def A(N,M,I,wv):
    s,r=len(M),len(I)
    if not geq(M,I): return 0
    n=s; u=[wv[m-1] for m in M]; t=0
    for perm in itertools.permutations(range(n)):
        us=[u[p] for p in perm]
        t+=summand(us,wv,I,N,n,1,0,q)/sp.prod([us[a]-us[b] for (a,b) in denom_pairs(n,r)])
    E=(1-q)**(-(n-r))*sp.prod([1-q*a/b for a in u for b in u])*sp.prod([1-q*a*b for a in u for b in wv])
    E/=sp.prod([1-q*u[j]*u[k] for j in range(n) for k in range(j+1,n)])
    E/=sp.prod([1-q*wv[j]*wv[k] for j in range(N) for k in range(j,N)])
    return sp.cancel(t*sp.prod([x**(-(N+n)) for x in u])/sp.factorial(n-r)/E)
def lower(nu,N):
    return [tuple(c) for k in range(len(nu)+1) for c in itertools.combinations(range(1,N+1),k) if geq(nu,tuple(c))]
def cfull(lam,mu,nu,N,a,b,seed=5):
    random.seed(seed)
    wv=[sp.Rational(random.randint(2,97),random.randint(2,97))+random.randint(0,3) for c in range(N)]
    if a: wv[b-1]=1/(wv[a-1]*(1-eps))
    L=sorted(lower(nu,N),key=lambda T:(len(T),sum(T)))
    Am={(M,I):A(N,M,I,wv) for M in L for I in L}
    c={}
    for M in L:
        v=(Am[M,lam] if lam in L else 0)*(Am[M,mu] if mu in L else 0)
        c[M]=sp.cancel((v-sum(c[k]*Am[M,k] for k in c if Am[M,k]!=0 and c[k]!=0))/Am[M,M])
    return c[nu]
if __name__=='__main__':
    lam=eval(sys.argv[1]);mu=eval(sys.argv[2]);nu=eval(sys.argv[3]);N=int(sys.argv[4])
    for a,b in itertools.combinations(range(1,N+1),2):
        cv=cfull(lam,mu,nu,N,a,b)
        r=sp.factor(sp.cancel(cv*eps).subs(eps,0)) if True else None
        num,den=sp.fraction(sp.cancel(cv))
        dpoly=sp.Poly(den,eps); ordd=min(m[0] for m in dpoly.monoms())
        print((a,b),'ord_eps(den)=',ordd, ' residue(q)=', sp.factor(sp.cancel(cv*eps**ordd).subs(eps,0)) if ordd>0 else 0, flush=True)
