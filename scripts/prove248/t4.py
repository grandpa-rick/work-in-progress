from psi import *
def Phi(kap,P,M):   # P: function of list x
    s=0
    for h in itertools.permutations(M,len(kap)):
        if any(h[j]<kap[j] for j in range(len(kap))): continue
        x=[w(m) for m in h]; s+=T0(kap,h,M)*P(x)
    return s
def g(u,Mm): return sp.Mul(*[(u-1/w(t))/(u-w(t)) for t in Mm])
# exact Lemma T*: Phi_kap[P](M) = sum_j s_j(u) Phi_{kap\kap_j}[P(..u..)prod_{i>=j}x_i](M-) + Tail
def check(kap,P,M):
    n=max(M); u=w(n); Mm=tuple(m for m in M if m!=n); r=len(kap)
    lhs=Phi(kap,P,M)
    rhs=0
    for j in range(r):
        J=j+1
        sj=(-1)**(r-J)*u**(-(r-J))*(1-u**-2)*sp.Mul(*[(1-w(k)/u) for k in range(1,kap[j])])*g(u,Mm)
        kk=kap[:j]+kap[j+1:]
        rhs+=sj*Phi(kk,lambda x: P(list(x[:j])+[u]+list(x[j:]))*sp.Mul(*x[j:]),Mm)
    rhs+=Phi(kap,lambda x:P(x)*sp.Mul(*[(u*xi-1)/(u*(xi-u)) for xi in x]),Mm)
    return sp.simplify(sp.cancel(sp.together(lhs-rhs)))
print(check((2,3),lambda x:x[0]**2*x[1]**-1,(2,3,4,5)))
print(check((1,3),lambda x:x[0]+3*x[1],(1,3,4,5)))
print(check((2,3,4),lambda x:1,(2,3,4,5)))
# D_u for B2: Phi_{(s,s+1)}[Q] heads = C(u) u Phi_{(s)}[-u x1(x1-u)]
s=2; M=(2,3,4,5); n=5; u=w(n); Mm=(2,3,4); a1,a2=w(s),w(s+1)
Q=lambda x: x[0]*x[1]-a2*x[0]+a1*a2
lhs=Phi((s,s+1),Q,M)-Phi((s,s+1),lambda x:Q(x)*sp.Mul(*[(u*xi-1)/(u*(xi-u)) for xi in x]),Mm)
C=(-1)**2*u**-2*(1-u**-2)*g(u,Mm)*sp.Mul(*[(1-w(k)/u) for k in range(1,s)])
rhs=C*u*Phi((s,),lambda x:-u*x[0]*(x[0]-u),Mm)
print(sp.simplify(sp.cancel(sp.together(lhs-rhs))))
