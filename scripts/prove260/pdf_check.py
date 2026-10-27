# PROVE 260: PDF-is-the-claim check of the FPSAC cut, re-implemented from the PRINTED statements.
# Independent of scripts/prove25x engines. LG(m,2m), type C Demazure operators on Laurent polys.
from fractions import Fraction as Fr
from itertools import product, combinations
from collections import defaultdict
import sympy as sp

def strict_parts(m):
    out=[]
    for k in range(m+1):
        for c in combinations(range(m,0,-1),k): out.append(tuple(c))
    return out

def D(i,f,m):  # i in 1..m ; f: dict exponent-tuple -> coeff ; e^lambda
    g=defaultdict(int)
    for lam,c in f.items():
        lam=list(lam)
        if i<m: n=lam[i-1]-lam[i]; al=[0]*m; al[i-1]=1; al[i]=-1
        else: n=lam[m-1]; al=[0]*m; al[m-1]=2
        if n>=0: ks=[(k,1) for k in range(0,n+1)]; sgn=-1
        elif n==-1: ks=[]
        else: ks=[(k,-1) for k in range(1,-n)]; sgn=1
        for k,s in ks:
            mu=tuple(lam[t]+sgn*k*al[t] for t in range(m)) if n>=0 else tuple(lam[t]+k*al[t] for t in range(m))
            g[mu]+=s*c
    return {k:v for k,v in g.items() if v}

def word(kappa,m):  # letters in order of application (first applied first)
    w=[]
    for k in kappa:  # c(kappa_1) acts first
        w+= list(range(m, m-k, -1))  # D_m first, then D_{m-1},...,D_{m-k+1}
    return w

def apply(kappa,f,m):
    for i in word(kappa,m): f=D(i,f,m)
    return f

def H_at(kappa,j,m):
    return sum(apply(kappa,{tuple([j]*m):1},m).values())

def Hpoly(kappa,m):
    J=sp.symbols('j'); d=sum(kappa)
    pts=[(j,H_at(kappa,j,m)) for j in range(d+1)]
    return sp.expand(sp.interpolate(pts,J)),J

def Fa(a,m):
    f={}
    for S in product([0,1],repeat=m):
        f[tuple(-a-S[i] for i in range(m))]=1
    return f

# --- combinatorics from printed definitions ---
def D_of(kappa): return [(a,b) for a in range(len(kappa)) for b in range(a,a+kappa[a])]
def leq(x,y): return x[0]<=y[0] and x[1]<=y[1]
def N_skew(kappa,nu,j):  # strict maps D(kappa)\D(nu) -> [2j-1], even on diagonal
    boxes=[x for x in D_of(kappa) if x not in set(D_of(nu))]
    boxes.sort(key=lambda x:(x[0]+x[1]))
    cnt=0
    def rec(i,f):
        nonlocal cnt
        if i==len(boxes): cnt+=1; return
        x=boxes[i]; lo=1
        for y in boxes[:i]:
            if leq(y,x): lo=max(lo,f[y]+1)
        for v in range(lo,2*j):
            if x[0]==x[1] and v%2: continue
            f[x]=v; rec(i+1,f)
        f.pop(x,None)
    rec(0,{})
    return cnt
def r_of(kappa): return 0 if not kappa else 1+max(a+b for a,b in D_of(kappa))
import math

m=3
table={  # printed Table 1 of the cut (H values)
(1,):[0,-1,-2,-3,-4, Fr(1,2),Fr(-1,2),Fr(-3,2),Fr(-5,2),Fr(-7,2)],
(2,):[0,1,4,9,16, Fr(1,4),Fr(1,4),Fr(9,4),Fr(25,4),Fr(49,4)],
(3,):[0,0,-3,-13,-34, Fr(1,8),Fr(1,8),Fr(-7,8),Fr(-55,8),Fr(-175,8)],
(2,1):[0,0,-1,-5,-14, Fr(1,4),0,Fr(-1,4),Fr(-5,2),Fr(-35,4)],
(3,1):[0,0,2,15,56, Fr(1,8),0,Fr(3,8),Fr(25,4),Fr(245,8)],
(3,2):[0,0,-1,-13,-70, Fr(1,8),0,Fr(-1,8),Fr(-17,4),Fr(-259,8)],
(3,2,1):[0,0,0,1,14, Fr(1,8),0,0,Fr(1,8),Fr(35,8)],
}
bad=0
for kappa,row in table.items():
    P,J=Hpoly(kappa,m); sz=sum(kappa); k1=kappa[0]; kp=kappa[1:]
    for t,j in enumerate([1,2,3,4,5]):
        h=Fr(int(P.subs(J,-j)))
        if h!=row[t]: print('TABLE int mismatch',kappa,j,h,row[t]); bad+=1
        if (-1)**sz*h!=N_skew(kappa,(),j): print('Thm1.1b fail',kappa,j); bad+=1
    for t,a in enumerate([0,1,2,3,4]):
        h=sp.Rational(P.subs(J,sp.Rational(-2*a-1,2)))
        hb=Fr(sum(apply(kappa,Fa(a,m),m).values()),2**m)   # Prop 6.1 bridge
        if Fr(h.p,h.q)!=row[5+t] or hb!=row[5+t]: print('TABLE half mismatch',kappa,a,h,hb,row[5+t]); bad+=1
        s=(-1)**sz*hb
        if a==0:
            if hb!=Fr(1,2**k1): print('ThmA fail',kappa); bad+=1
            continue
        # Cor 6.5 row formula
        rf=N_skew(kappa,(),a)+sum(Fr((-1)**(L+1),2**L)*N_skew(kappa,(L,),a) for L in range(1,k1+1))
        if rf!=s: print('row formula fail',kappa,a,rf,s); bad+=1
        # Thm 1.5
        Nk=N_skew(kp,(),a) if kp else 1
        zero = bool(kp) and a<=math.ceil(r_of(kp)/2)
        if (s==0)!=zero: print('HZ(i) fail',kappa,a); bad+=1
        if s!=0:
            sg = (-1)**(sz-k1) if 2*a<=k1 else (-1)**sz
            if (2*a<k1 or 2*a>k1) and (1 if hb>0 else -1)!=sg: print('HZ(ii) fail',kappa,a); bad+=1
        if abs(s)<Fr(Nk,2**k1): print('HZ(iii) bound fail',kappa,a); bad+=1
        if 2*a<=k1 and abs(s)!=Fr(Nk,2**k1): print('HZ(iii) eq fail',kappa,a); bad+=1
        if (s*2**k1).denominator!=1: print('HZ integrality fail'); bad+=1
    # Thm 1.4 window from table
    W=math.ceil(r_of(kappa)/2); sj=[(-1)**sz*row[t] for t in range(5)]
    for j in range(1,6):
        if (sj[j-1]==0)!=(j<=W): print('window fail',kappa,j); bad+=1
print('table/thm checks done, failures:',bad)
# Example 1.2
print('ex1.2 counts j=1..4:',[N_skew((2,1),(),j) for j in range(1,5)])
# Remark 4.4 V-values for kappa=(3,2), and Remark 5.1 recursion
def V(theta):
    th=list(theta); c=0
    for f in product('eo',repeat=len(th)):
        ok=all(not(th[i][0]==th[i][1] and f[i]=='o') for i in range(len(th)))
        for i in range(len(th)):
            for k in range(len(th)):
                if i!=k and leq(th[i],th[k]) and not (f[i]=='e' and f[k]=='o'): ok=False
        c+=ok
    return c
K=(3,2); Dk=set(D_of(K))
Vs={rho:V(Dk-set(D_of(rho))) for rho in strict_parts(3) if set(D_of(rho))<Dk}
print('V for (3,2):',Vs)
s=lambda kap,j:(-1)**sum(kap)*table[kap][j-1] if kap else 1
for j in (2,3):
    print('s_%d(3,2) via (*):'%(j+1), s(K,j)+sum(v*s(rho,j) for rho,v in Vs.items()), ' table:',s(K,j+1))
# check (*) as polynomial identity for all kappa, m=3, at half-integers too
for kappa in strict_parts(3):
    if not kappa: continue
    P,J=Hpoly(kappa,3); Dk=set(D_of(kappa))
    rhs=(-1)**sum(kappa)*P.subs(J,-J)
    tot=rhs
    for rho in strict_parts(3):
        if set(D_of(rho))<Dk:
            Q,_=Hpoly(rho,3) if rho else (sp.Integer(1),J)
            tot+=V(Dk-set(D_of(rho)))*(-1)**sum(rho)*Q.subs(J,-J)
    lhs=(-1)**sum(kappa)*P.subs(J,-J-1)
    if sp.expand(lhs-tot)!=0: print('(*) poly fail',kappa)
print('(*) polynomial identity checked m=3')
# Remark 4.4: c_rho=(-1)^{|theta|-1}V(theta) in the hyperplane identity (4.2), m=3,4, all kappa
for mm in (3,4):
    fails=0
    for kappa in strict_parts(mm):
        if not kappa: continue
        P,J=Hpoly(kappa,mm); Dk=set(D_of(kappa)); rhs=0
        for rho in strict_parts(mm):
            if set(D_of(rho))<Dk:
                th=Dk-set(D_of(rho)); Q=Hpoly(rho,mm)[0] if rho else sp.Integer(1)
                rhs+=(-1)**(len(th)-1)*V(th)*Q
        if sp.expand(P-P.subs(J,J-1)-rhs)!=0: fails+=1; print('hyperplane fail',kappa)
    print('hyperplane identity with printed c_rho, m=%d: failures %d'%(mm,fails))
