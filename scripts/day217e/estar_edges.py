"""Day 217e: e*_mu = E_{mu1}...E_{mul}(1) from Hikita operator matrices (mats_N5, b-basis), m-expanded;
test (A) lim_{s->oo} s^{-n(mu)} e*_mu = HL P_{mu'}(x;t); (B) e*_mu|_{t=0} = s^{n(mu)} P_{mu'}(x;1/s,0).
Macdonald P computed independently by Gram-Schmidt in the p-basis."""
import sys, pickle, itertools, sympy as sp
sys.path.insert(0,'/home/agent/projects/scripts/day214')
from ek_subset_engine import transpose, parts, dominates
s,t,q,T=sp.symbols('s t q T')
mats=pickle.load(open('scripts/day216/mats_N5.pkl','rb'))
N=int(sys.argv[1])
def nn(m): return sum(i*x for i,x in enumerate(m))
def apply(k,vec):
    res={}
    for mu,c in vec.items():
        for nu,d in mats[(k,sum(mu))][mu].items(): res[nu]=res.get(nu,0)+c*d
    return {a:sp.cancel(b) for a,b in res.items() if sp.cancel(b)!=0}
# symmetric function machinery in m-basis via power sums
from sympy.utilities.iterables import partitions
def plist(n): return [tuple(sorted(sum(([k]*v for k,v in p.items()),[]),reverse=True)) for p in partitions(n)]
def mono_coeffs_of_product(gens,n):
    # expand product of p_r or e_r in n variables, return m-coeffs
    pass
def m_of(expr,xs):
    P=sp.Poly(sp.expand(expr),*xs); d={}
    for mon,c in P.terms():
        if list(mon)==sorted(mon,reverse=True): d[tuple(a for a in mon if a>0)]=c
    return d
def zee(l):
    from collections import Counter
    r=1
    for k,v in Counter(l).items(): r*=k**v*sp.factorial(v)
    return r
cache={}
def macP(n):
    if n in cache: return cache[n]
    xs=sp.symbols(f'x1:{n+1}'); Ps=plist(n)
    pm={l:m_of(sp.prod([sum(x**r for x in xs) for r in l]),xs) for l in Ps}  # p_l in m-basis
    M=sp.Matrix([[pm[l].get(m,0) for m in Ps] for l in Ps])  # rows p, cols m
    Minv=M.inv()  # m_m = sum Minv[m_idx, l] p_l
    def ip(f,g):  # f,g dicts in m-basis
        fp=[sum(f.get(m,0)*Minv[j,i] for j,m in enumerate(Ps)) for i in range(len(Ps))]
        gp=[sum(g.get(m,0)*Minv[j,i] for j,m in enumerate(Ps)) for i in range(len(Ps))]
        return sum(fp[i]*gp[i]*zee(l)*sp.prod([(1-q**r)/(1-T**r) for r in l]) for i,l in enumerate(Ps))
    order=sorted(Ps)  # lex increasing = linear ext of dominance increasing
    P={}
    for lam in order:
        f={lam:sp.Integer(1)}
        for mu in order:
            if mu==lam: break
            if dominates(lam,mu):
                c=sp.cancel(ip({lam:1},P[mu])/ip(P[mu],P[mu]))
                for k,v in P[mu].items(): f[k]=sp.cancel(f.get(k,0)-c*v)
        # Gram-Schmidt against all lower in order (only dominated matter; others orthogonal anyway)
        P[lam]={k:sp.factor(v) for k,v in f.items() if v!=0}
    cache[n]=(P,xs); return cache[n]
ok=0;bad=0
for n in range(1,N+1):
    P,xs=macP(n)
    for mu in plist(n):
        v={():sp.Integer(1)}
        for k in reversed(mu): v=apply(k,v)
        # b_nu = s^{n(nu)} e_nu -> m-basis
        f=sum(c*s**nn(nu)*sp.prod([sum(sp.prod(cc) for cc in itertools.combinations(xs,r)) for r in nu]) for nu,c in v.items())
        md={k:sp.cancel(c) for k,c in m_of(sp.together(f)*1,xs).items()} if False else None
        num,den=sp.fraction(sp.cancel(sp.together(f)))
        md={k:sp.cancel(c/den) for k,c in m_of(num,xs).items()}
        lp=transpose(mu)
        HL={k:sp.cancel(v.subs(q,0)) for k,v in P[lp].items()}
        A={k:sp.cancel(sp.limit(c*s**(-nn(mu)),s,sp.oo)) for k,c in md.items()}
        TT=t if len(sys.argv)<3 else 1/t
        targA={k:HL.get(k,0).subs(T,TT) for k in set(A)|set(HL)}
        okA=all(sp.cancel(A.get(k,0)-targA[k])==0 for k in targA)
        B={k:sp.cancel(sp.limit(c,t,0)) for k,c in md.items()}
        QW={k:sp.cancel(v.subs(T,0).subs(q,1/s)*s**nn(mu)) for k,v in P[lp].items()}
        okB=all(sp.cancel(B.get(k,0)-QW.get(k,0))==0 for k in set(B)|set(QW))
        print(mu,'A',okA,'B',okB,flush=True)
        ok+=okA+okB; bad+=(not okA)+(not okB)
print('OK',ok,'BAD',bad)
