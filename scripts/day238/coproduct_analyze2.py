import pickle, sympy as sp
exec(open('coproduct_analyze.py').read().split('# (i) witness')[0].replace("print(","(lambda *a,**k:None)("))
def all_rowadd(nu):
    out=[]
    for d in range(1,sum(nu)):
        for a in parts(d):
            for b in parts(sum(nu)-d):
                if rowadd2(a,b,nu): out.append((a,b))
    return out
def limits(g):
    Lg={}
    for nu,dP in R.items():
        for (a,b),c in dP.items():
            if a==() or b==(): continue
            v=sp.factor(lim(c*w(g,nu,a,b)))
            if v!=0: Lg[(nu,a,b)]=v
    return Lg
LN=limits(gN)
miss=[(nu,a,b) for nu in R for (a,b) in all_rowadd(nu) if (nu,a,b) not in LN]
print('every row-additive split present:', not miss, miss)
# primitive-in-limit trivial terms: coefficient of nu⊗1 should be 1
print('nu⊗1 coefficients:', set(sp.simplify(R[nu][(nu,())]) for nu in R))
# coboundary solve: unknown g(lam), g(1)=1 normalisation; free generators = nu with no nontrivial split
G={(): sp.Integer(1)}
for d in range(1,5):
    for nu in parts(d):
        sp_=[(a,b) for (n2,a,b) in LN if n2==nu]
        if not sp_: G[nu]=sp.Integer(1); continue   # free: choose 1
        vals=set(sp.simplify(G[a]*G[b]/LN[(nu,a,b)]) for a,b in sp_)
        vals=set(sp.factor(v) for v in vals)
        assert len(vals)==1, (nu,vals)
        G[nu]=vals.pop()
print('(ii) coboundary g exists (all splits consistent):'); 
for k,v in G.items(): print('   g',k,'=',v)
# sanity: one-row g should be s-dependent 1/[n]_s! * s^{-C(n,2)} type
# controls
LC=limits(gC)
nonra=[k for k in LC if not rowadd2(k[1],k[2],k[0])]
print('CONTROL gauge t^|nu|: non-row-additive splits surviving =',len(nonra),'(must be >0)')
L0={}
for nu,dP in R.items():
    for (a,b),c in dP.items():
        if a and b:
            v=lim(c)
            if v!=0: L0[(nu,a,b)]=v
print('CONTROL identity gauge: non-row-additive surviving =',len([k for k in L0 if not rowadd2(k[1],k[2],k[0])]))
# swap control: does the N-limit support equal union-splits alpha ∪ beta = nu (Delta m_nu)?
def union(a,b,nu): return tuple(sorted(a+b,reverse=True))==nu
sw=[k for k in LN if not union(k[1],k[2],k[0])]
print('CONTROL swap nu<->nu\': N-limit splits NOT of form α∪β=ν:',len(sw),'(must be >0)')
