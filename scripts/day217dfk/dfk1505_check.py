# Literal check of DFK arXiv:1505.01657 eq (5.15)/(5.25)/(5.27) vs Theorem B, N=4 variables, |mu|<=4.
import itertools, sympy as sp
N=4; z=sp.symbols('z1:%d'%(N+1)); q,t=sp.symbols('q t')
def M(k,f,n=1):  # DFK (5.15): sum_{|I|=k} z_I^n a_I(z) Gamma_I, a_I=prod_{i in I,j notin I} z_i/(z_i-z_j), Gamma_i: z_i->q z_i
    s=0
    for I in itertools.combinations(range(N),k):
        a=sp.Mul(*[z[i]/(z[i]-z[j]) for i in I for j in range(N) if j not in I])
        s+=sp.Mul(*[z[i]**n for i in I])*a*f.subs({z[i]:q*z[i] for i in I},simultaneous=True)
    return sp.expand(sp.cancel(sp.together(s)))
def E0(k,f,s=q):  # our E_k at t=0: sum_A prod x_i/(x_i-x_j) X_A F(X_{A^c}, sX_A)  (written independently)
    tot=0
    for A in itertools.combinations(range(N),k):
        Ac=[j for j in range(N) if j not in A]
        w=sp.Mul(*[z[i]/(z[i]-z[j]) for i in A for j in Ac])*sp.Mul(*[z[i] for i in A])
        tot+=w*f.xreplace({z[i]:s*z[i] for i in A})
    return sp.expand(sp.cancel(sp.together(tot)))
def parts(n,mx=None):
    mx=n if mx is None else mx
    if n==0: yield (); return
    for p in range(min(n,mx),0,-1):
        for r in parts(n-p,p): yield (p,)+r
def conj(l): return tuple(sum(1 for x in l if x>i) for i in range(l[0])) if l else ()
def nfun(l): return sum(i*x for i,x in enumerate(l))
def mono(l):
    l=tuple(l)+(0,)*(N-len(l))
    return sp.Add(*[sp.Mul(*[z[i]**e[i] for i in range(N)]) for e in set(itertools.permutations(l))])
def schur(l):
    l=list(l)+[0]*(N-len(l))
    if len(l)>N: return 0
    num=sp.Matrix(N,N,lambda i,j: z[i]**(l[j]+N-1-j)); den=sp.Matrix(N,N,lambda i,j: z[i]**(N-1-j))
    return sp.expand(sp.cancel(num.det()/den.det()))
def coeffs(f,d):
    P=sp.Poly(sp.expand(f),*z); return [P.coeff_monomial(sp.Mul(*[z[i]**(list(l)+[0]*N)[i] for i in range(N)])) for l in B[d]]
B={d:[l for l in parts(d) if len(l)<=N] for d in range(5)}
def D1(f):  # Macdonald D_1 (q,t) symbolic
    s=0
    for i in range(N):
        c=sp.Mul(*[(t*z[i]-z[j])/(z[i]-z[j]) for j in range(N) if j!=i])
        s+=c*f.subs(z[i],q*z[i])
    return sp.expand(sp.cancel(sp.together(s)))
def macP(lam,qq):  # P_lam(z;qq,0): eigenvector of D1 at generic (q,t), then t->0, q->qq
    d=sum(lam); idx=B[d].index(lam); low=B[d][idx:]  # B in reverse-lex => dominance-lower are later
    cs=sp.symbols('c1:%d'%len(low)); f=mono(lam)+sum(c*mono(l) for c,l in zip(cs,low[1:]))
    lp=list(lam)+[0]*(N-len(lam)); ev=sum(q**lp[i]*t**(N-1-i) for i in range(N))
    eqs=[sp.expand(sp.numer(sp.together(e))) for e in coeffs(sp.expand(D1(f)-ev*f),d)]
    sol=sp.solve([e for e in eqs if e!=0],cs,dict=True)[0] if cs else {}
    P=sp.expand(f.subs(sol))
    P=sum(sp.factor(sp.limit(c,t,0))*mono(l) for c,l in zip(coeffs(P,d),B[d]))
    return sp.expand(P.subs(q,qq))
KF={ (1,):{(1,):1}, (2,):{(2,):1},(1,1):{(2,):t,(1,1):1},
 (3,):{(3,):1},(2,1):{(3,):t,(2,1):1},(1,1,1):{(3,):t**3,(2,1):t+t**2,(1,1,1):1},
 (4,):{(4,):1},(3,1):{(4,):t,(3,1):1},(2,2):{(4,):t**2,(3,1):t,(2,2):1},
 (2,1,1):{(4,):t**3,(3,1):t+t**2,(2,2):t,(2,1,1):1},
 (1,1,1,1):{(4,):t**6,(3,1):t**3+t**4+t**5,(2,2):t**2+t**4,(2,1,1):t+t**2+t**3,(1,1,1,1):1}}
print('mu | DFK prod M_{mu_i,1}.1 == (a) q^n P_mu\'(1/q,0) | == (b) sum Ktilde s_lam\' | == (c) our E_k(t=0) | (a)==(b) | is P_mu\'(q,0)? ')
for d in range(1,5):
  for mu in parts(d):
    L=sp.Integer(1)
    for k in mu: L=M(k,L)
    a=sp.expand(q**nfun(mu)*macP(conj(mu),1/q))
    b=sp.expand(sum(sp.expand(q**nfun(mu)*sp.sympify(K).subs(t,1/q))*schur(conj(lam)) for lam,K in KF[mu].items()))
    C=sp.Integer(1)
    for k in reversed(mu): C=E0(k,C)
    wrong=sp.expand(macP(conj(mu),q))
    z0=lambda u: sp.simplify(u)==0
    print(mu, z0(L-a), z0(L-b), z0(L-C), z0(a-b), z0(L-wrong))
    if d==2 and mu==(1,1): print('   e1*e1 (DFK) =',sp.factor(L))
