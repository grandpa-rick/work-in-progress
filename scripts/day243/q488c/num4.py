"""Independent engine check of c1 poles at N=4 (and N=3 control).
A_{M,I} from Motegi's symmetrisation summand (engine.summand), unpadded n=|M| (stability), numeric random w,
except approach to H_ab:  w_b = 1/(w_a (1-eps)).  q symbolic -> A0=A|q=0, A1=dA/dq|q=0 as rational fns of eps.
c0,c1 by triangular solve; Res_H c1^nu := lim eps*c1^nu.  Not using any of the day-243 theory."""
import itertools, random, sys, sympy as sp
from engine import summand, denom_pairs, geq, subsets
q,eps=sp.symbols('q eps')
def A01(N,M,I,wv):
    s,r=len(M),len(I)
    if not geq(M,I): return 0,0
    n=s; u=[wv[m-1] for m in M]; t0=0; t1=0
    for perm in itertools.permutations(range(n)):
        us=[u[p] for p in perm]
        t=summand(us,wv,I,N,n,1,0,q)
        d=sp.prod([us[a]-us[b] for (a,b) in denom_pairs(n,r)])
        t0+=t.subs(q,0)/d; t1+=sp.diff(t,q).subs(q,0)/d
    E=(1-q)**(-(n-r))*sp.prod([1-q*a/b for a in u for b in u])*sp.prod([1-q*a*b for a in u for b in wv])
    E/=sp.prod([1-q*u[j]*u[k] for j in range(n) for k in range(j+1,n)])
    E/=sp.prod([1-q*wv[j]*wv[k] for j in range(N) for k in range(j,N)])
    E1=sp.diff(E,q).subs(q,0)
    pref=sp.prod([x**(-(N+n)) for x in u])/sp.factorial(n-r)
    return sp.cancel(pref*t0), sp.cancel(pref*(t1-t0*E1))
def run(N,a,b,seed=7):
    random.seed(seed)
    wv=[sp.Rational(random.randint(2,97),random.randint(2,97))+random.randint(0,3) for c in range(N)]
    wv[b-1]=1/(wv[a-1]*(1-eps))
    S=subsets(N); A0={};A1={}
    for M in S:
        for I in S: A0[M,I],A1[M,I]=A01(N,M,I,wv)
    order=sorted(S,key=lambda T:(len(T),sum(T)))
    def solve(v):
        c={}
        for nu in order:
            c[nu]=sp.cancel((v[nu]-sum(c[k]*A0[nu,k] for k in c if A0[nu,k]!=0 and c[k]!=0))/A0[nu,nu])
        return c
    out={}
    for l in S:
        for m in S:
            if order.index(m)<order.index(l): continue
            c0=solve({M:A0[M,l]*A0[M,m] for M in S})
            Phi={M:A0[M,l]*A1[M,m]+A1[M,l]*A0[M,m]-sum(c0[k]*A1[M,k] for k in S if c0[k]!=0) for M in S}
            c1=solve(Phi)
            for nu in S:
                if not(geq(nu,l) and geq(nu,m)): continue
                r1=sp.cancel(c1[nu]*eps)
                pole1 = sp.limit(r1,eps,0)!=0 if r1.has(eps) or r1!=0 else False
                r0=sp.cancel(c0[nu]*eps); pole0 = sp.limit(r0,eps,0)!=0 if r0!=0 else False
                out[l,m,nu]=(pole0,pole1, c0[nu].subs(eps,0)==0 if not sp.cancel(c0[nu]).has(eps) or True else None)
    return out
if __name__=='__main__':
    N=int(sys.argv[1]); a,b=int(sys.argv[2]),int(sys.argv[3])
    out=run(N,a,b)
    for k,v in out.items():
        l,m,nu=k; j=len(nu)-max(len(l),len(m))
        print(l,m,nu,'jump',j,'c0pole',v[0],'c1pole',v[1],flush=True)
