# (A) Jing CT: X^lam_mu = [z^lam] K'(z) p_mu(z), K' = prod_{i<j}(1-z_j/z_i)/(1-t z_j/z_i) expanded in z_j/z_i.
# (B) ell(lam)=3 collapse of JL sum via exponential formula.
# (C) slice-leakage: contribution of level-1 classes of length>=3 in (2.32) at two-part mu.
import sys; sys.path.insert(0,'../day226')
import sympy as sp, itertools
from green import X as Xg, parts, t
from jl import XJL, zt, subparts
def Kcoef_table(l, N):
    # coefficients of K' as dict gamma->coef, truncated: each factor sum_k c_k (z_j/z_i)^k, k<=N
    z=sp.symbols('z1:%d'%(l+1)); res={tuple([0]*l):sp.Integer(1)}
    for i in range(l):
        for j in range(i+1,l):
            new={}
            for g,c in res.items():
                for k in range(N+1):
                    ck = 1 if k==0 else (t**k - t**(k-1))
                    gg=list(g); gg[i]-=k; gg[j]+=k; gg=tuple(gg)
                    if any(abs(v)>N for v in gg): continue
                    new[gg]=new.get(gg,0)+c*ck
            res=new
    return res
def XCT(lam,mu):
    l=len(lam); n=sum(lam); K=Kcoef_table(l,n)
    tot=0
    for assign in itertools.product(range(l),repeat=len(mu)):
        g=list(lam)
        for k,i in enumerate(assign): g[i]-=mu[k]
        tot+=K.get(tuple(g),0)
    return sp.expand(tot)
bad=0;cnt=0
for n in range(1,7):
    for lam in parts(n):
        for mu in parts(n):
            cnt+=1
            if sp.simplify(XCT(lam,mu)-Xg(lam,mu))!=0: bad+=1; print('CT MISMATCH',lam,mu)
print('(A) Jing CT vs green:',cnt,'bad',bad)
# (B) ell(lam)=3: X = sum_{tau ⊳ mu,|tau|<=m} sum_i [s^i] prod_{tau}(1+s^tau_j) [u^{m-|tau|}] E(u,s) * w_i
u,s=sp.symbols('u s')
def X3(lam,mu):
    l1,l2,l3=lam; m=l2+l3
    E=(1-u)*(1-s*u)/((1-t*u)*(1-t*s*u))
    Es=sp.series(E,u,0,m+1).removeO()
    tot=0
    for tau in subparts(mu):
        if sum(tau)>m: continue
        k=m-sum(tau)
        f=sp.expand(sp.prod([1+s**x for x in tau])*Es.coeff(u,k))
        # X^{(l2,l3)}_nu = (t-1) sum_{i>=l2+1} D^(i) t^{i-l2-1} + D^(l3)
        for i in range(0,m+1):
            ci=f.coeff(s,i)
            w=(t-1)*t**(i-l2-1) if i>=l2+1 else 0
            if i==l3: w+=1
            tot+=ci*w
    return sp.factor(sp.simplify(tot))
bad=0;cnt=0
for n in range(3,10):
    for lam in parts(n):
        if len(lam)!=3: continue
        for x in range(1,n):
            y=n-x
            if x<y: continue
            cnt+=1
            if sp.simplify(X3(lam,(x,y))-Xg(lam,(x,y)))!=0: bad+=1; print('X3 MISMATCH',lam,(x,y))
print('(B) ell(lam)=3 collapse vs green, two-part classes n<=9:',cnt,'bad',bad)
# (C) leakage
for lam,mu in [((2,2,2),(3,3)),((3,2,1),(4,2)),((2,2,1,1),(3,3)),((3,3,2),(5,3))]:
    m=sum(lam[1:]); leak=0; ins=0; nterms=0
    for tau in subparts(mu):
        if sum(tau)>m: continue
        for rho in parts(m-sum(tau)):
            cls=tuple(sorted(tau+rho,reverse=True)); c=(-1)**len(rho)/zt(rho)*Xg(lam[1:],cls); nterms+=1
            if len(cls)>=3: leak+=c
            else: ins+=c
    print('(C)',lam,mu,'level-1 terms',nterms,' in-slice part',sp.factor(sp.simplify(ins)),' leak (len>=3)',sp.factor(sp.simplify(leak)))
