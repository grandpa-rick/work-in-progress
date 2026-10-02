"""check e*_mu|_{s=0} = [n;mu]_t e_n and t^{n(mu')-C(n,2)} e*_mu -> (1-s)^{l-1} e_n as t->oo"""
import sys,pickle,itertools,sympy as sp
sys.path.insert(0,'/home/agent/projects/scripts/day214')
from ek_subset_engine import transpose, parts
s,t=sp.symbols('s t'); mats=pickle.load(open('scripts/day216/mats_N5.pkl','rb'))
def nn(m): return sum(i*x for i,x in enumerate(m))
def apply(k,vec):
    res={}
    for mu,c in vec.items():
        for nu,d in mats[(k,sum(mu))][mu].items(): res[nu]=res.get(nu,0)+c*d
    return {a:sp.cancel(b) for a,b in res.items() if sp.cancel(b)!=0}
def tf(n): return sp.prod([sum(t**j for j in range(i)) for i in range(1,n+1)])
ok=bad=0
for n in range(1,6):
    for mu in parts(n):
        v={():sp.Integer(1)}
        for k in reversed(mu): v=apply(k,v)
        e={nu:sp.cancel(c*s**nn(nu)) for nu,c in v.items()}  # e-basis coefficients of e*_mu
        A={nu:sp.cancel(sp.limit(c,s,0)) for nu,c in e.items()}
        A={k:x for k,x in A.items() if x!=0}
        predA={(n,):sp.cancel(tf(n)/sp.prod([tf(r) for r in mu]))}
        B={nu:sp.cancel(sp.limit(c*t**(nn(transpose(mu))-sp.binomial(n,2)),t,sp.oo)) for nu,c in e.items()}
        B={k:x for k,x in B.items() if x!=0}
        predB={(n,):sp.expand((1-s)**(len(mu)-1))}
        a=all(sp.cancel(A.get(k,0)-predA.get(k,0))==0 for k in set(A)|set(predA))
        b=all(sp.cancel(B.get(k,0)-predB.get(k,0))==0 for k in set(B)|set(predB))
        print(mu,a,b,flush=True); ok+=a+b; bad+=(not a)+(not b)
print('OK',ok,'BAD',bad)
