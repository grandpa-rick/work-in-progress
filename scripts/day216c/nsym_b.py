"""Day 216c: nonsymmetric test. Is gamma X_i gamma^{-1} = c * Y^bullet_j (Hikita/207b twisted Y)?
std rep: pi F = F(X2..Xm, s X1); bullet: pi F = X1 F(X2..Xm, s X1).  Y_i = t^{m-i} T_{i-1}..T_1 pi T_{m-1}^{-1}..T_i^{-1}."""
import sympy as sp, itertools, sys
from sympy import Rational as R
m=int(sys.argv[1]); d=int(sys.argv[2])
s=R(5,7); t=R(-3,11)
X=sp.symbols(f'x1:{m+1}')
def T(i,F):  # i 1-based
    a,b=X[i-1],X[i]; sF=F.subs({a:b,b:a},simultaneous=True)
    return sp.expand(t*sF+sp.cancel((t-1)*b*(sF-F)/(a-b)))
def Tinv(i,F): return sp.expand((T(i,F)-(t-1)*F)/t)
def pi(F,bul):
    G=F.subs({X[j]:(X[j+1] if j<m-1 else s*X[0]) for j in range(m)},simultaneous=True)
    return sp.expand(X[0]*G if bul else G)
def Y(i,F,bul):
    for j in range(i,m): F=Tinv(j,F)
    F=pi(F,bul)
    for j in range(1,i): F=T(j,F)
    return sp.expand(t**(m-i)*F)
def mons(deg): return [c for c in itertools.product(range(deg+1),repeat=m) if sum(c)==deg]
def vec(F,deg):
    P=sp.Poly(F,*X); D=dict(P.terms()); return sp.Matrix([D.get(c,0) for c in mons(deg)])
def mat(op,deg,deg2):
    B=mons(deg); return sp.Matrix.hstack(*[vec(op(sp.prod([x**e for x,e in zip(X,c)])),deg2) for c in B])
def fac(r):
    r=sp.Rational(r); a=0;b=0;n,dd=r.p,r.q
    for p,sign in ((2,1),(3,1)):
        pass
    def v(n,p):
        k=0
        while n%p==0: n//=p;k+=1
        return k
    a=v(abs(n),2)-v(dd,2); b=v(abs(n),3)-v(dd,3)
    assert sp.Rational(2)**a*sp.Rational(3)**b*sp.sign(n)==r,(r)
    return a,b
def eig(deg):
    Ys=[mat(lambda F,i=i:Y(i,F,False),deg,deg) for i in range(1,m+1)]
    lam=sp.Symbol('L'); vecs=[]
    cand=[sorted(sp.Poly(Yi.charpoly(lam).as_expr(),lam).ground_roots().keys()) for Yi in Ys]
    import itertools as it
    for ys in it.product(*cand):
        M=sp.Matrix.vstack(*[Ys[i]-ys[i]*sp.eye(Ys[0].rows) for i in range(m)])
        ns=M.nullspace()
        for v in ns: vecs.append((v,list(ys)))
    assert len(vecs)==Ys[0].rows,(len(vecs),Ys[0].rows)
    return vecs
res={}
for sg in (1,):
  for sb in (1,):
    ok=True
    for deg in range(d):
        E0=eig(deg); E1=eig(deg+1)
        def gam(ys):
            g=R(1)
            for y in ys:
                a,b=fac(y); g*= s**(sg*a*(a-1)//2)*t**(sb*a*b)
            return g
        P0=sp.Matrix.hstack(*[v for v,_ in E0]); P1=sp.Matrix.hstack(*[v for v,_ in E1])
        G0=sp.diag(*[gam(y) for _,y in E0]); G1=sp.diag(*[gam(y) for _,y in E1])
        Gam0=P0*G0*P0.inv(); Gam1=P1*G1*P1.inv()
        for i in range(1,m+1):
            Xi=mat(lambda F:sp.expand(X[i-1]*F),deg,deg+1)
            Gi=Gam1*Xi*Gam0.inv()
            found=[]
            for j in range(1,m+1):
                Yb=mat(lambda F,j=j:Y(j,F,True),deg,deg+1)
                # proportional?
                nz=[(a,b) for a in range(Yb.rows) for b in range(Yb.cols) if Yb[a,b]!=0]
                if not nz: continue
                c=Gi[nz[0]]/Yb[nz[0]]
                if Gi==c*Yb: found.append((j,c))
            print(sg,sb,deg,i,found,flush=True)
