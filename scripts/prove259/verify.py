# PROVE 259: verify (V1) strict character formula, (V2) LS/G-side formula, (V3) P_comb = S C^{-1} S via n(j+1)=P n(j),
# (V4) negative control (perturbed beta).  Weights as integer tuples; characters as dicts.
import sys, itertools
from collections import Counter
from nchar import strict, n as nch
sys.path.insert(0,'../prove249d'); from dem import Di, word
m=int(sys.argv[1]); J=int(sys.argv[2]); pert=len(sys.argv)>3
def D(k): return [(a,b) for a in range(len(k)) for b in range(a,a+k[a])]
def sub(r,k): return len(r)<=len(k) and all(r[i]<=k[i] for i in range(len(r)))
def eps(c):  # c in 1..m
    v=[0]*m; v[c-1]=1; return v
def beta(x):
    d=x[1]-x[0]
    if d==0: return eps(m)
    v=[a-b for a,b in zip(eps(m-d),eps(m+1-d))]
    if pert and d==1: v=[a+b for a,b in zip(eps(m-d),eps(m+1-d))]  # negative control
    return v
def wom(k):
    v=[1]*m
    for p in k: v[m-p]=-1
    return v
def add(*vs): return tuple(sum(t) for t in zip(*vs))
def sc(c,v): return [c*t for t in v]
def maps(k,lo,hi,strictp):
    B=D(k); rel=[(x,y) for x in B for y in B if x!=y and x[0]<=y[0] and x[1]<=y[1]]
    B=sorted(B,key=lambda x:(x[0]+x[1],x))  # linear extension
    def rec(i,f):
        if i==len(B): yield dict(f); return
        x=B[i]; low=max([f[y]+(1 if strictp else 0) for y in f if y[0]<=x[0] and y[1]<=x[1]]+[lo])
        for v in range(low,hi+1):
            if x[0]==x[1] and v%2: continue
            f[x]=v; yield from rec(i+1,f); del f[x]
    yield from rec(0,{})
def ch(k,j):
    f={tuple([j]*m):1}
    for i in reversed(word(set(k),m)): f=Di(f,i,m)
    return Counter({e:c for e,c in f.items() if c})
ok1=ok2=True
for k in strict(m):
    for j in range(1,J+1):
        c=Counter()
        for f in maps(k,1,2*j-1,True):
            c[add(sc(-j,wom(k)),*[sc(-f[x],beta(x)) for x in D(k)])]+=1
        if Counter({e:v for e,v in nch(k,j,m).items() if v})!=+c: ok1=False; print('V1 FAIL',k,j)
    for j in range(0,J):
        c=Counter()
        for G in maps(k,0,2*j,False):
            c[add(sc(j,wom(k)),*[sc(G[x],beta(x)) for x in D(k)])]+=1
        if ch(k,j)!=+c: ok2=False; print('V2 FAIL',k,j)
print('m',m,'V1 strict-side formula:',ok1,' V2 LS-side formula:',ok2)
# V3: weight-preserving peel: n_k(j+1) = sum_rho P[k,rho] n_rho(j), P[k,rho] = sum_g x^{-w_k om - sum_{g=o} beta}
def Vmaps(th):
    th=list(th)
    for vals in itertools.product((0,1),repeat=len(th)):
        g=dict(zip(th,vals))
        if any(g[x] for x in th if x[0]==x[1]): continue
        if all(not(x!=y and x[0]<=y[0] and x[1]<=y[1]) or g[x]<g[y] for x in th for y in th): yield g
K=strict(m)
P={}
for k in K:
    for r in K:
        if not sub(r,k): continue
        th=set(D(k))-set(D(r)); c=Counter()
        for g in Vmaps(th): c[add(sc(-1,wom(k)),*[sc(-1,beta(x)) for x in th if g[x]==1])]+=1
        P[(k,r)]=c
def mul(a,b):
    c=Counter()
    for e1,v1 in a.items():
        for e2,v2 in b.items(): c[add(e1,e2)]+=v1*v2
    return c
ok3=True
for k in K:
    for j in range(1,J):
        lhs=Counter({e:v for e,v in nch(k,j+1,m).items() if v}); rhs=Counter()
        for r in K:
            if (k,r) in P: rhs+=mul(P[(k,r)],Counter(nch(r,j,m)))  # Counter + drops <=0; all positive here
        if lhs!=rhs: ok3=False; print('V3 FAIL',k,j)
print('V3 peel weight-preserving n(j+1)=P n(j):',ok3)
if m==3:
    for key in [((2,),(1,)),((3,1),(2,)),((3,2),(3,1))]:
        print(key,dict(P[key]))
