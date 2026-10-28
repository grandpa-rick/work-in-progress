# PROVE 261 PDF check: re-implement printed statements of §4.5 (ssec:equivariant) from scratch.
import sys, itertools, random
from fractions import Fraction
from collections import defaultdict
m=int(sys.argv[1]); PERT = len(sys.argv)>2
def add(a,b,c=1):
    r=defaultdict(int,a)
    for k,v in b.items(): r[k]+=c*v
    return {k:v for k,v in r.items() if v}
def mono(e,c=1): return {tuple(e):c}
def alpha(i):  # simple root i=1..m
    v=[0]*m
    if i<m: v[i-1]=1; v[i]=-1
    else: v[m-1]=2
    return v
def pair(l,i): return l[i-1]-l[i] if i<m else l[m-1]
def Dm(f,i):
    a=alpha(i); r={}
    for l,c in f.items():
        n=pair(l,i); t={}
        if n>=0:
            for k in range(n+1): t[tuple(l[q]-k*a[q] for q in range(m))]=c
        elif n<=-2:
            for k in range(1,-n): t[tuple(l[q]+k*a[q] for q in range(m))]=-c
        r=add(r,t)
    return r
def mul(a,b):
    r=defaultdict(int)
    for k,v in a.items():
        for k2,v2 in b.items(): r[tuple(x+y for x,y in zip(k,k2))]+=v*v2
    return {k:v for k,v in r.items() if v}
# strict partitions, diagrams
parts=[tuple(sorted(c,reverse=True)) for r in range(m+1) for c in itertools.combinations(range(1,m+1),r)]
def D(k): return [(a,a+d) for a in range(len(k)) for d in range(k[a])]
def v(k): return [(-1 if (m+1-(c+1)) in k else 1) for c in range(m)]  # position c+1 -> -1 iff m+1-(c+1) in k... check below
def vv(k):
    s=[1]*m
    for p in k: s[m+1-p-1]=-1
    return s
def beta(x):
    a,b=x; d=b-a; e=[0]*m
    if d==0: e[m-1]=1
    else:
        e[m-d-1]+=1; e[m-d]-=1
        if PERT and d==1: e=[0]*m; e[m-2]=1; e[m-1]=1
    return e
# ch via Demazure: build w_kappa by adding boxes along a linear extension: box offset d -> alpha_{m-d} (d>=1), alpha_m (d=0)
def word(k):
    # add boxes row by row? need each intermediate an order ideal: add by increasing a+b
    boxes=sorted(D(k),key=lambda x:(x[0]+x[1],x[0]))
    w=[]
    for (a,b) in boxes:
        d=b-a; w.append(m if d==0 else m-d)
    return w   # D_w = D_{w[-1]} o ... o D_{w[0]} applied: first-added acts first
def ch(k,j):
    f=mono([j]*m)
    for i in word(k): f=Dm(f,i)
    return f
def leq(x,y): return x[0]<=y[0] and x[1]<=y[1]
def maps(P,vals,strict,diag_even):
    P=list(P); res=[]
    def rec(i,cur):
        if i==len(P): res.append(dict(cur)); return
        x=P[i]
        for val in vals:
            if diag_even and x[0]==x[1] and val%2: continue
            ok=True
            for y,vy in cur.items():
                if leq(y,x) and y!=x and (vy>val if not strict else vy>=val): ok=False;break
                if leq(x,y) and y!=x and (val>vy if not strict else val>=vy): ok=False;break
            if ok:
                cur[x]=val; rec(i+1,cur); del cur[x]
    rec(0,{}); return res
def wsum(base,G):
    e=list(base)
    for x,g in G.items():
        b=beta(x); e=[e[q]+g*b[q] for q in range(m)]
    return e
ok=True; fails=defaultdict(int)
# printed v_kappa check vs w_kappa omega: lowest term of ch(k,1) (G=0) should be e^{v_kappa}
for k in parts:
    for j in range(0,4):
        lhs=ch(k,j); rhs={}
        for G in maps(D(k),range(0,2*j+1),False,True):
            rhs=add(rhs,mono(wsum([j*c for c in vv(k)],G)))
        if lhs!=rhs: fails['Tchar']+=1
    for j in range(1,4):
        lhs={kk:(-1)**sum(k)*c for kk,c in ch(k,-j).items()}; rhs={}
        for f in maps(D(k),range(1,2*j),True,True):
            rhs=add(rhs,mono([-c for c in wsum([j*c for c in vv(k)],f)]))
        if lhs!=rhs: fails['strict']+=1
# P matrix from printed definition
def Pent(k,r):
    th=[x for x in D(k) if x not in D(r)]
    if any(len(r)>len(k) or (i<len(r) and r[i]>k[i]) for i in range(len(r))) : return {}
    out={}
    for g in maps(th,[0,1],True,False):
        if any(g[x]==1 for x in th if x[0]==x[1]): continue
        e=[-c for c in vv(k)]
        for x in th:
            if g[x]==1: b=beta(x); e=[e[q]-b[q] for q in range(m)]
        out=add(out,mono(e))
    return out
def sub(r,k): return len(r)<=len(k) and all(r[i]<=k[i] for i in range(len(r)))
P={(k,r):Pent(k,r) for k in parts for r in parts if sub(r,k)}
for k in parts:
    for i in range(-2,4):
        rhs={}
        for r in parts:
            if sub(r,k): rhs=add(rhs,mul(P[(k,r)],ch(r,i)),(-1)**(sum(k)-sum(r)))
        if ch(k,i-1)!=rhs: fails['wpeel']+=1
# independence: numeric rank at random point
random.seed(1); pt=[Fraction(random.randint(2,30),random.randint(2,30)) for _ in range(m)]
def ev(f): 
    s=Fraction(0)
    for e,c in f.items():
        t=Fraction(c)
        for q in range(m): t*=pt[q]**e[q]
        s+=t
    return s
M=[[ev(ch(r,j)) for r in parts] for j in range(1,len(parts)+1)]
def rank(M):
    M=[row[:] for row in M]; rk=0; cols=len(M[0])
    for c in range(cols):
        piv=next((i for i in range(rk,len(M)) if M[i][c]!=0),None)
        if piv is None: continue
        M[rk],M[piv]=M[piv],M[rk]
        for i in range(len(M)):
            if i!=rk and M[i][c]!=0:
                fct=M[i][c]/M[rk][c]; M[i]=[a-fct*b for a,b in zip(M[i],M[rk])]
        rk+=1
    return rk
rk=rank(M); M1=[[sum(ch(r,j).values()) for r in parts] for j in range(1,len(parts)+1)]
print('m',m,'PERT',PERT,'fails',dict(fails),'equivariant rank',rk,'/',len(parts),'rank at e=1',rank([[Fraction(x) for x in row] for row in M1]))
# Theorem Tchev at e^mu=1: c_rho = (-1)^{|theta|-1} V(theta); print kappa=(3,2) at m=3 and P_(2),(1)
if m==3 and not PERT:
    k=(3,2); print('c for (3,2):',{r:-(-1)**(sum(k)-sum(r))*sum(P[(k,r)].values()) for r in parts if sub(r,k) and r!=k})
    print('P_(2),(1):',P[((2,),(1,))]); print('A_kk=e^{v_k} check', all(P[(k,k)]==mono([-c for c in vv(k)]) for k in parts))
    # Lemma hyperplane check: H_k(j)-H_k(j-1) = sum c H_r(j)
    bad=0
    for k in parts:
        for j in range(-3,4):
            H=lambda r,j: sum(ch(r,j).values())
            lhs=H(k,j)-H(k,j-1); rhs=sum(-(-1)**(sum(k)-sum(r))*sum(P[(k,r)].values())*H(r,j) for r in parts if sub(r,k) and r!=k)
            bad+= lhs!=rhs
    print('hyperplane identity failures',bad)
