import sys, re, sympy as sp
sys.path.insert(0,'../prove247')
import mp
c_=sp.Symbol('c')
N=int(sys.argv[1]); fn=sys.argv[2]; cs=eval(sys.argv[3])
inst=[]
for line in open(fn):
    m=re.match(r'\((\d+), (\d+)\) (\(.*?\)) (\(.*?\)) (\(.*?\)) d',line)
    if not m or ' P ' not in line: continue
    tl=lambda x: x if isinstance(x,tuple) else (x,)
    a,b=int(m.group(1)),int(m.group(2)); lam,mu,nu=[tl(eval(m.group(i))) for i in (3,4,5)]
    inst.append(((a,b,lam,mu,nu),sp.sympify(line.split(' P ',1)[1])))
systems={}
def S(a,b,seed):
    if (a,b,seed) not in systems: systems[(a,b,seed)]=mp.Sys(N,seed,a,b)
    return systems[(a,b,seed)]
for cc in cs:
    zX=set(); zP=set()
    for (a,b,lam,mu,nu),P in inst:
        if P.subs(c_,cc)==0: zP.add((a,b,lam,mu,nu))
        if all(mp.A(S(a,b,sd),a,b,lam,mu,nu,cc)==0 for sd in (11,12)): zX.add((a,b,lam,mu,nu))
    print('c',cc,'#X=0',len(zX),'#P(c)=0',len(zP),'X-zeros not P-zeros',len(zX-zP),'P-zeros not X-zeros',len(zP-zX),flush=True)
