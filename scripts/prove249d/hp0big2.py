# HP for lam=mu=0 at larger n via Theorem D / Demazure localisation at y = 1 + delta*a (exact), chi(X_kappa,O(m)) integers
import sys, random, itertools, sympy as sp
from fractions import Fraction as Fr
from dem import word
n=int(sys.argv[1]); random.seed(int(sys.argv[2]) if len(sys.argv)>2 else 1)
dl=Fr(1,10**40)
y0=tuple(1+dl*random.randint(1,97)*random.choice([1,-1]) for _ in range(n))
def act(i,y):
    y=list(y)
    if i<n: y[i-1],y[i]=y[i],y[i-1]
    else: y[n-1]=1/y[n-1]
    return tuple(y)
orbit={y0}; fr=[y0]
while fr:
    nf=[]
    for p in fr:
        for i in range(1,n+1):
            q=act(i,p)
            if q not in orbit: orbit.add(q); nf.append(q)
    fr=nf
orbit=list(orbit)
def ema(i,y): return y[i]/y[i-1] if i<n else y[n-1]**-2
def sv(i,v):
    v=list(v)
    if i<n: v[i-1],v[i]=v[i],v[i-1]
    else: v[n-1]=-v[n-1]
    return tuple(v)
def Dop(cls,i):
    keys=set(cls)|{sv(i,v) for v in cls}; out={}
    for vp in keys:
        a=cls.get(vp); b=cls.get(sv(i,vp)); d={}
        for y in orbit:
            e=ema(i,y); d[y]=((a[y] if a else 0)-e*(b[act(i,y)] if b else 0))/(1-e)
        out[vp]=d
    return out
c=sp.Symbol('c'); t=sp.Symbol('t'); U=2*n-2
subsl=[k for r in range(n+1) for k in itertools.combinations(range(1,n+1),r)]
if len(sys.argv)>3: subsl=[k for k in subsl if sum(k)>=int(sys.argv[3])]
bad=0
for kap in subsl:
    cls={tuple([1]*n):{y:Fr(1) for y in orbit}}
    for i in reversed(word(set(kap),n)): cls=Dop(cls,i)
    d=sum(kap); vals=[]
    for m in range(0,d+3):
        S=Fr(0)
        for v,co in cls.items():
            p=co[y0]
            for e,yy in zip(v,y0): p*=yy**(m*e)
            S+=p
        vals.append(round(S))
    h=sp.interpolate(list(zip(range(d+3),vals)),sp.Symbol('m'))
    P=sp.expand((-1)**d*h.subs(sp.Symbol('m'),-c/2-1))
    Q=sp.Poly(sp.expand(P.subs(c,U+1+t)),t).all_coeffs()
    nz=[q for q in Q if q!=0]; hp=(Q[-1]!=0) and (all(q>0 for q in nz) or all(q<0 for q in nz))
    try:
        roots=sp.Poly(P,c).nroots(n=30,maxsteps=500) if d>0 else []
    except Exception: roots=[]
    mx=max([complex(r).real for r in roots]) if roots else None
    over=(sp.degree(sp.interpolate(list(zip(range(d+3),vals)),sp.Symbol('m')))==d)
    from dem import chi as _chi
    if n<=4: assert vals==[_chi(set(kap),n,m)[0] for m in range(d+3)], (kap,vals)
    if not hp: bad+=1
    Q0=sp.Poly(sp.expand(P.subs(c,U+t)),t).all_coeffs(); nz0=[q for q in Q0 if q!=0]; hp0=all(q>0 for q in nz0) or all(q<0 for q in nz0)
    print(kap,'dim',d,'deg_ok',over,'HP+',hp,'HP0+',hp0,'P',sp.factor(P),flush=True)
print('n',n,'HP+ failures',bad)
