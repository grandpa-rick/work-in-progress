"""Ask 1, fast variant: pole order of FULL c^nu_{lam,mu} along every H_ab, N=5 (c)-failure triples.
Same construction as day243/q488c/fullc.cfull (Motegi summand, unpadded n=|M|, E^C normalisation) but with q
specialised to a random rational (two values) so only eps stays symbolic.  w generic rationals, w_b=1/(w_a(1-eps))."""
import sys, os, itertools, random, sympy as sp
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '../../day243/q488c'))
from engine import summand, denom_pairs, geq
eps=sp.Symbol('eps')
def A(N,M,I,wv,qq):
    s,r=len(M),len(I)
    if not geq(M,I): return 0
    n=s; u=[wv[m-1] for m in M]; t=0
    for perm in itertools.permutations(range(n)):
        us=[u[p] for p in perm]
        t+=summand(us,wv,I,N,n,1,0,qq)/sp.prod([us[a]-us[b] for (a,b) in denom_pairs(n,r)])
    E=(1-qq)**(-(n-r))*sp.prod([1-qq*a/b for a in u for b in u])*sp.prod([1-qq*a*b for a in u for b in wv])
    E/=sp.prod([1-qq*u[j]*u[k] for j in range(n) for k in range(j+1,n)])
    E/=sp.prod([1-qq*wv[j]*wv[k] for j in range(N) for k in range(j,N)])
    return sp.cancel(t*sp.prod([x**(-(N+n)) for x in u])/sp.factorial(n-r)/E)
def lower(nu,N):
    return [tuple(c) for k in range(len(nu)+1) for c in itertools.combinations(range(1,N+1),k) if geq(nu,tuple(c))]
def cnum(lam,mu,nu,N,a,b,qq,seed=5):
    random.seed(seed)
    wv=[sp.Rational(random.randint(2,97),random.randint(2,97))+random.randint(0,3) for c in range(N)]
    wv[b-1]=1/(wv[a-1]*(1-eps))
    L=sorted(lower(nu,N),key=lambda T:(len(T),sum(T)))
    Am={(M,I):A(N,M,I,wv,qq) for M in L for I in L}
    c={}
    for M in L:
        v=(Am[M,lam] if lam in L else 0)*(Am[M,mu] if mu in L else 0)
        c[M]=sp.cancel((v-sum(c[k]*Am[M,k] for k in c if Am[M,k]!=0 and c[k]!=0))/Am[M,M])
    return c[nu]
def order(cv):
    if cv==0: return 'zero'
    n_,d_=sp.fraction(sp.cancel(cv))
    on=min(m[0] for m in sp.Poly(n_,eps).monoms()); od=min(m[0] for m in sp.Poly(d_,eps).monoms())
    return max(od-on,0)
N=5
T=eval(open(os.path.join(os.path.dirname(os.path.abspath(__file__)),'../../day243/q488c/n5fails_partial.txt')).read())
for idx in map(int,sys.argv[1:]):
    lam,mu,nu=T[idx]
    for a,b in itertools.combinations(range(1,N+1),2):
        o=[order(cnum(lam,mu,nu,N,a,b,qq,seed=s)) for qq,s in ((sp.Rational(3,7),5),(sp.Rational(-5,11),9))]
        print(idx,lam,mu,nu,(a,b),'pole order (q=3/7,seed5; q=-5/11,seed9):',o,flush=True)
