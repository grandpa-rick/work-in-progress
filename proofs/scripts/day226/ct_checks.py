# Day 226 independent checks: (i) CT[P_lam Q_mu(1/z) K] = phi_a delta (my orthogonality route to Thm 1.1)
# (ii) Phi_3(p_3^2;7,2) directly from the CT side of Thm 1.1 vs the Thm 2.5 value used in the kill test.
from sympy import *
from itertools import permutations
t=symbols('t')
def phi(m): return prod([1-t**i for i in range(1,m+1)])
def ct_K(H, z):
    a=len(z); H=expand(H)
    terms=Poly(H*prod(z)**60, *z).terms()   # shift to make polynomial
    # height h(beta)=sum_i i*beta_i (1-indexed); x_ij=z_i/z_j has h=i-j<0; need total h of H-monomial + K-monomial = 0 at beta=0
    res=0
    # build truncated K
    maxh=max(sum((i+1)*(e[i]-60) for i in range(a)) for e,_ in terms)
    Kt=1
    for i in range(a):
        for j in range(i+1,a):
            x=z[i]/z[j]; d=j-i; kmax=max(0,maxh//d)
            ser=1+sum((t**k-t**(k-1))*x**k for k in range(1,kmax+1))
            Kt=expand(Kt*ser)
    prodH=expand(H*Kt)
    # constant term
    s=0
    for term in Add.make_args(prodH):
        c,m = term.as_independent(*z)
        if m==1: s+=c
    return simplify(s)
def P(lam,z):
    a=len(z); lam=list(lam)+[0]*(a-len(lam))
    from collections import Counter
    vl=prod([prod([(1-t**i) for i in range(1,m+1)])/(1-t)**m for m in Counter(lam).values()])
    s=0
    for w in permutations(range(a)):
        zz=[z[w[i]] for i in range(a)]
        s+= prod([zz[i]**lam[i] for i in range(a)])*prod([(zz[i]-t*zz[j])/(zz[i]-zz[j]) for i in range(a) for j in range(i+1,a)])
    return factor(cancel(s/vl))
def b(lam):
    from collections import Counter
    return prod([phi(m) for p,m in Counter(lam).items() if p>0])
z=symbols('z1:4'); a=3
lams=[(2,1,1),(3,1,1),(2,2,1)]
mus=[(2,1,1),(4,),(3,1),(2,2),(2,1,1),(3,1,1),(2,2,1),(5,),(4,1)]
for lam in lams:
    Pl=expand(P(lam,z))
    for mu in mus:
        if sum(mu)!=sum(lam): continue
        Qm=b(mu)*P(mu,z); Qinv=Qm.subs({zi:1/zi for zi in z},simultaneous=True)
        v=ct_K(Pl*Qinv,z)
        print(lam,mu, factor(v), 'expect', factor(phi(a)) if lam==mu else 0)
# (ii)
g=sum(zi**3 for zi in z)**2
H=prod(z)*g*sum(zi**-7 for zi in z)*sum(zi**-2 for zi in z)
v=ct_K(H,z)
Phi=(1-t**7)*(1-t**2)/phi(3)*v
curly = t**-2/((1-t)*(1-t**2)) + t**-4*(1+t**3)**2 - t**-4*(1+t**3+t**6)**2*(1+t**2+t**4)/(1-t**3)
print('Phi direct-CT minus Thm2.5 value:', simplify(Phi - (-(1-t**7)*(1-t**2)*curly)))
