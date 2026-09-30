# Day 213: test log K_ij bilinear-affine in (X,Y)=(t^{-in},t^{jn}); exact exponents only.
import sympy as sp
s,t,u=sp.symbols('s t u')
N=7  # orders u^1..u^6
def K(i,j,z,w):
    e=sp.Integer(1)
    for p in range(i): e*= (t**p*z-s*w)/(t**p*z-t**j*w)
    for r in range(j): e*= (s*z-t**r*w)/(t**i*z-t**r*w)
    return e
def logcoeffs(i,j):
    # K(1,u) = s^j t^{-ij} * F(u), F(0)=1. Taylor F, then formal log. Independent of factor-wise logs.
    Ku=sp.cancel(K(i,j,1,u))
    c0=sp.simplify(Ku.subs(u,0))
    assert sp.simplify(c0 - s**j*t**(-i*j))==0, (i,j,c0)
    F=sp.series(Ku/c0,u,0,N).removeO()
    f=[sp.cancel(F.coeff(u,k)) for k in range(N)]; assert sp.simplify(f[0]-1)==0
    # formal log: L' = F'/F ; l_k = f_k - (1/k) sum_{m=1}^{k-1} m l_m f_{k-m}
    l=[0]*N
    for k in range(1,N):
        l[k]=sp.cancel(f[k]-sum(m*l[m]*f[k-m] for m in range(1,k))/k)
    return [sp.cancel(k*l[k]) for k in range(N)]  # coefficient of u^n/n
pairs=[(i,j) for i in range(4) for j in range(4)]
data={p:logcoeffs(*p) for p in pairs}
train=[(1,1),(1,2),(2,1),(2,2)]
held=[p for p in pairs if p not in train]
a,b,c,d=sp.symbols('a b c d')
ok=True
for n in range(1,N):
    eqs=[a+b*t**(-i*n)+c*t**(j*n)+d*t**(-i*n)*t**(j*n)-data[(i,j)][n] for (i,j) in train]
    sol=sp.solve(eqs,[a,b,c,d],dict=True)[0]
    sol={k:sp.factor(v) for k,v in sol.items()}
    T,S=t**n,s**n
    pred={a:(T*S-1/S)/(1-T),b:(1-T*S)/(1-T),c:(1/S-T)/(1-T),d:sp.Integer(-1)}
    dream=all(sp.simplify(sol[k]-pred[k])==0 for k in pred)
    bad=[p for p in held if sp.simplify(sol[a]+sol[b]*t**(-p[0]*n)+sol[c]*t**(p[1]*n)+sol[d]*t**(-p[0]*n+p[1]*n)-data[p][n])!=0]
    ok&= dream and not bad
    print(n,sol,'matches_dream_formula=',dream,'heldout_failures=',bad)
print('ALL OK' if ok else 'FAIL')
