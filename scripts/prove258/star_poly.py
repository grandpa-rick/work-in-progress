# PROVE 258: (star) as a polynomial identity, checked with Demazure Hilbert polynomials (no BR, no Q30 code).
# s_x(k) = (-1)^|k| H_k(-x). Check s_{x+1}(k) - s_x(k) - sum_{rho<k} V(D(k)\D(rho)) s_x(rho) == 0 identically.
# Negative control: perturb one V entry (+1) and find the first integer j>=1 where the identity breaks.
import sys, itertools
import sympy as sp
sys.path.insert(0,'/home/agent/projects/work-in-progress/scripts/prove249d')
from dem import chi as chiC
x=sp.Symbol('x')
def D(k): return {(i,i+j) for i,li in enumerate(k) for j in range(li)}
def V(th):
    th=list(th); n=len(th)
    tot=0
    for lab in itertools.product((0,1),repeat=n):   # 0=e<1=o
        ok=True
        for s in range(n):
            a,b=th[s]
            if a==b and lab[s]==1: ok=False;break
            for t in range(n):
                if s!=t and th[s][0]<=th[t][0] and th[s][1]<=th[t][1] and not lab[s]<lab[t]: ok=False;break
            if not ok: break
        if ok: tot+=1
    return tot
m=int(sys.argv[1])
parts=[tuple(sorted(c,reverse=True)) for L in range(m+1) for c in itertools.combinations(range(1,m+1),L)]
H={}
for k in parts:
    d=sum(k); pts=[(j,chiC(set(k),m,j)[0]) for j in range(d+1)]
    H[k]=sp.expand(sp.interpolate(pts,x)) if d>0 else sp.Integer(pts[0][1])
s=lambda k,arg: (-1)**sum(k)*H[k].subs(x,-arg)
sub=lambda r,k: all(i<len(k) and r[i]<=k[i] for i in range(len(r))) and r!=k
Vt={(k,r):V(D(k)-D(r)) for k in parts for r in parts if sub(r,k)}
bad=0
for k in parts:
    e=sp.expand(s(k,x+1)-s(k,x)-sum(Vt[(k,r)]*s(r,x) for r in parts if sub(r,k)))
    if e!=0: bad+=1; print('FAIL',k,e)
print('m',m,'shapes',len(parts),'polynomial (star) failures',bad)
# negative control: every perturbed entry with s_j(rho) not identically 0 should break at some integer j>=1
nc=0; caught=0; firstj=[]
for (k,r),v in Vt.items():
    nc+=1
    for j in range(1,40):
        lhs=s(k,j+1); rhs=s(k,j)+sum((Vt[(k,q)]+(1 if q==r else 0))*s(q,j) for q in parts if sub(q,k))
        if lhs!=rhs: caught+=1; firstj.append(j); break
print('negative control: perturbed entries',nc,'caught at an integer j>=1:',caught,'max first j',max(firstj) if firstj else None)
