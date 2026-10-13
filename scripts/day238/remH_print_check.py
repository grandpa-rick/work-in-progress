# Day 238: print check of the edited Rem H (longversion) — from printed statements only.
# (a) Kirillov factorisation M(e,P)_{lam,kap}(t) = sum_nu K_{nu,lam} K_{nu',kap}(t)
# (b) Cor dmatrix 1st formula == 2nd formula; d(1) = #0-1 matrices rows lam, cols mu'
# (c) corner mu=(n): d_{lam,(n)} = [n]_t!/prod [lam_i]_t!   (Hikita Lemma 6.3)
import itertools, sympy as sp
from sympy.utilities.iterables import partitions
t=sp.symbols('t')
def parts(n): return [tuple(sorted(sum(([k]*v for k,v in p.items()),[]),reverse=True)) for p in partitions(n)]
def conj(l): return tuple(sum(1 for x in l if x>i) for i in range(l[0])) if l else ()
def nfun(l): return sum(i*x for i,x in enumerate(l))
def HLP(lam,xs,u):
    m=len(xs); lam=tuple(lam)+(0,)*(m-len(lam)); tot=0
    seen=set()
    for w in itertools.permutations(range(m)):
        key=tuple(lam[w.index(i)] for i in range(m))
        if key in seen: continue
        seen.add(key)
        # sum over S_m/S_m^lam: w(x^lam prod_{i<j,lam_i>lam_j} (x_i-u x_j)/(x_i-x_j))
        term=sp.Mul(*[xs[w[i]]**lam[i] for i in range(m)])
        for i in range(m):
            for j in range(i+1,m):
                if lam[i]>lam[j]: term*= (xs[w[i]]-u*xs[w[j]])/(xs[w[i]]-xs[w[j]])
        tot+=term
    return sp.factor(sp.together(tot))
def ek(k,xs): return sum(sp.Mul(*c) for c in itertools.combinations(xs,k))
def monos(f,xs):
    return sp.Poly(sp.expand(f),*xs)
def coeff_in_basis(f,basis,xs,lam_list):
    # expand f in basis (dict lam->poly) by leading dominant monomial elimination via linear solve on monomials x^lam
    cs=sp.symbols('c0:%d'%len(lam_list))
    g=sp.expand(f-sum(c*basis[l] for c,l in zip(cs,lam_list)))
    P=sp.Poly(g,*xs)
    eqs=P.coeffs()
    sol=sp.solve(eqs,cs,dict=True)[0]
    return {l:sp.factor(sol[c]) for c,l in zip(cs,lam_list)}
def kostka(nu,lam):
    # number of SSYT shape nu weight lam via e/h? use s_nu = det h: compute by brute monomial coefficient
    return None
def tfact(n,u): return sp.Mul(*[sum(u**j for j in range(i)) for i in range(1,n+1)])
def zero_one(rows,cols):
    from functools import lru_cache
    rows=list(rows); cols=tuple(cols)
    def rec(i,cols):
        if i==len(rows): return 1 if all(c==0 for c in cols) else 0
        tot=0
        for S in itertools.combinations(range(len(cols)),rows[i]):
            if all(cols[j]>0 for j in S):
                c=list(cols)
                for j in S: c[j]-=1
                tot+=rec(i+1,tuple(c))
        return tot
    return rec(0,cols)
ok=True
for n in range(1,5):
    xs=sp.symbols('x1:%d'%(n+1)); L=parts(n)
    PT={l:HLP(l,xs,t) for l in L}; PTi={l:HLP(l,xs,1/t) for l in L}
    S={l:HLP(l,xs,0) for l in L}  # Schur = P at t=0
    for lam in L:
        el=sp.Mul(*[ek(k,xs) for k in lam])
        Me=coeff_in_basis(el,PT,xs,L)          # M(e,P)_{lam,.}(t)
        Mi=coeff_in_basis(el,PTi,xs,L)         # at t^{-1}
        Ks=coeff_in_basis(el,S,xs,L)           # K_{nu' lam}: e_lam = sum_nu K_{nu' lam} s_nu
        for kap in L:
            # (a) sum_nu K_{nu,lam} K_{nu',kap}(t) ; K_{nu lam} = Ks[nu'] ; K_{nu' kap}(t)=[P_kap] s_{nu'}
            rhs=0
            for nu in L:
                Knl=Ks[conj(nu)]
                if Knl==0: continue
                rhs+=Knl*coeff_in_basis(S[conj(nu)],PT,xs,L)[kap]
            if sp.simplify(Me[kap]-rhs)!=0: ok=False; print('(a) FAIL',lam,kap,Me[kap],rhs)
        for mu in L:
            mup=conj(mu)
            d1=sp.factor(t**(nfun(mup)-nfun(conj(lam)))*Mi[mup])
            # 2nd formula: t^{-n(lam')} sum_nu K_{nu' lam} Ktilde_{nu mu'}(t), Ktilde=t^{n(mu')}K(t^{-1})
            d2=0
            for nu in L:
                K=coeff_in_basis(S[nu],PTi,xs,L)[mup]   # K_{nu mu'}(t^{-1})
                d2+=Ks[nu]*t**nfun(mup)*K
            d2=sp.factor(t**(-nfun(conj(lam)))*d2)
            if sp.simplify(d1-d2)!=0: ok=False; print('(b) FAIL',lam,mu,d1,d2)
            if sp.simplify(d1.subs(t,1)-zero_one(lam,mup))!=0: ok=False; print('(b1) FAIL',lam,mu)
            if not sp.Poly(sp.expand(d1),t).is_polynomial if False else False: pass
            if mu==(n,):
                h=sp.factor(tfact(n,t)/sp.Mul(*[tfact(k,t) for k in lam]))
                if sp.simplify(d1-h)!=0: ok=False; print('(c) FAIL',lam,d1,h)
    print('n=',n,'done')
print('ALL OK' if ok else 'FAILURES')
