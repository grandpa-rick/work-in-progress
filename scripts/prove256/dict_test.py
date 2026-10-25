# PROVE 256: PDF-is-the-claim test of the indexing dictionary (draft Lemma lem:dictionary)
# Inputs ONLY from printed statements: Table 1 polynomials (draft), Thm Q30 (draft), BR Thm 5.5 (printed recursion), Lemma S.
from fractions import Fraction as F
from itertools import combinations
import sympy as sp
j=sp.symbols('j')
m=3
H={(): sp.Integer(1),(1,):j+1,(2,):(j+1)**2,(3,):(j+1)*(j+2)*(4*j+3)/6,(2,1):(j+1)*(j+2)*(2*j+3)/6,
   (3,1):(j+1)**2*(j+2)*(2*j+3)/6,(3,2):(j+1)*(j+2)*(2*j+3)*(2*j**2+6*j+5)/30,
   (3,2,1):(j+1)*(j+2)**2*(j+3)*(2*j+3)*(2*j+5)/180}
def D(k): return {(a,b) for a in range(len(k)) for b in range(a,a+k[a])}
parts=list(H)
def dual(k): return tuple(sorted(set(range(1,m+1))-set(k),reverse=True))
def le(x,y): return x[0]<=y[0] and x[1]<=y[1]
# Thm Q30 as printed
def V(th):
    th=list(th)
    for x in th:
        for y in th:
            for z in th:
                if x!=y!=z and x!=z and le(x,y) and le(y,z): return 0
    for (c,d) in th:
        if c==d and (c-1,c) in th: return 0
    iso=[x for x in th if x[0]!=x[1] and not any(y!=x and (le(x,y) or le(y,x)) for y in th)]
    return 2**len(iso)
def cQ30(k,r):
    th=D(k)-D(r); return (-1)**(len(th)-1)*V(th)
# BR Thm 5.5 as printed (BR coordinates: same (row,col) shifted boxes)
def isrim(th): return not any(x!=y and x[0]<y[0] and x[1]<y[1] for x in th for y in th)
def isrow(th): return len({x[0] for x in th})==1
def iscol(th): return len({x[1] for x in th})==1
def C(th,p):
    th=frozenset(th)
    if not isrim(th): return 0
    if p<=0: return 1 if not th else 0
    if not th: return 0
    r=min(x[0] for x in th); c=max(x[1] for x in th if x[0]==r); n=len(th)
    best=None
    for s in range(1,2*m+2):
        sq={x for x in th if r<=x[0]<=r+s-1 and c-s+1<=x[1]<=c}
        if isrow(sq) or iscol(sq):
            if best is None or len(sq)>len(best): best=sq
    psi=best; hat=th-psi; a=len(psi)
    diag=any(x[0]==x[1] for x in th)
    if not hat:
        if diag: return (int(p==n)-int(p==n-1)) if iscol(th) else int(p==n)
        return 2*int(p==n)-int(p==n-1)
    conn=any(abs(x[0]-y[0])+abs(x[1]-y[1])==1 for x in psi for y in hat)
    if conn: return C(hat,p-a)-C(hat,p-a+1)
    if a==1: return 2*C(hat,p-1)-2*C(hat,p)
    return 2*C(hat,p-a)-3*C(hat,p-a+1)+C(hat,p-a+2)
def cBR(k,r,dic):
    lam,nu=dic(k),dic(r)
    if not D(lam)<=D(nu): return 0
    return C(D(nu)-D(lam),1)
ident=dual
swap=lambda k: dual({(3,):(2,1),(2,1):(3,)}.get(k,k))
ok=True
for k in parts:
    for r in parts:
        if r!=k and D(r)<=D(k):
            a,b=cQ30(k,r),cBR(k,r,ident)
            if a!=b: ok=False; print('MISMATCH',k,r,a,b)
print('Q30 vs BR Thm5.5 via kappa^vee, all pairs m=3:', ok)
for name,dic in [('kappa^vee',ident),('swapped (3)<->(2,1)',swap)]:
    bad=0
    for k in parts:
        lhs=sp.expand(H[k]-H[k].subs(j,j-1))
        rhs=sp.expand(sum(cBR(k,r,dic)*H[r] for r in parts if r!=k))
        if sp.simplify(lhs-rhs)!=0: bad+=1
    print(name,': hyperplane identity fails for',bad,'of',len(parts),'kappa')
print('kappa=(3,1):',{r:cBR((3,1),r,ident) for r in parts if r!=(3,1) and cBR((3,1),r,ident)})
