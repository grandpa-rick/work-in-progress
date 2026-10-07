"""Day 228 Item 5: two-row specialisation of Thm 2.5 (two-point formula) vs green.py, |lam|<=10.
(R) raw: Thm 2.5 with a=2, A=B=1, G_1(w)=P_rho(1,w;t), Phi=(1-t^x)(1-t^y) X/b  =>  X = b*{...}.
(C) closed form (hand-derived): for y<=x, X = (t-1)t^{l2-1-y}(1+t^y) [y<l2]; m_xy-(1-t)t^{l2-1} [y=l2]; (t-1)t^{l2-1} [y>l2]."""
import sys; sys.path.insert(0,'/home/agent/projects/proofs/scripts/day226')
import sympy as sp
from green import X as Xgreen, b, t
w=sp.Symbol('w')
def G1(rho):
    r1,r2=rho; m=r1-r2
    if m==0: return w**r2
    return sp.expand(w**r2*(1+w**m+(1-t)*sum(w**k for k in range(1,m))))
def Xraw(lam,x,y):
    rho=(lam[0]-1,lam[1]-1); G=sp.Poly(G1(rho),w); Gk=lambda k: G.coeff_monomial(w**k) if k>=0 else 0
    a=2;A=1;B=1
    inner=Gk(y-B)/((1-t**A)*(1-t**B)) + sum((t**(-A*mm)-t**(B*mm))*Gk(y-B-mm) for mm in range(1,y+1))/(1-t**a)
    brace=t**(-A*B)*inner - G1(rho).subs(w,t)/(1-t**a)*sum(t**(-j*y) for j in range(a))
    Phi=(-1)**a*(1-t**x)*(1-t**y)*brace
    return sp.factor(sp.simplify(Phi*b(lam)/((1-t**x)*(1-t**y))))
def Xclosed(lam,x,y):
    l1,l2=lam
    if y>x: x,y=y,x
    if y<l2: return (t-1)*t**(l2-1-y)*(1+t**y)
    if y==l2: return (2 if x==y else 1)-(1-t)*t**(l2-1)
    return (t-1)*t**(l2-1)
bad=0;cnt=0
for n in range(2,11):
    for l2 in range(1,n//2+1):
        lam=(n-l2,l2)
        for y in range(1,n):
            x=n-y; cls=tuple(sorted((x,y),reverse=True))
            g=sp.expand(Xgreen(lam,cls)); r=sp.expand(Xraw(lam,x,y)); c=sp.expand(Xclosed(lam,x,y))
            cnt+=1
            if sp.simplify(g-r)!=0 or sp.simplify(g-c)!=0: bad+=1; print('FAIL',lam,(x,y),g,r,c)
print(f'checked {cnt} (lam,(x,y)) ordered pairs, n<=10; failures {bad}')
