import pickle, sympy as sp
exec(open('coproduct_test.py').read().split('P={}; Pp={}')[0])
Pp,R=pickle.load(open('/tmp/day238_coprod.pkl','rb'))
# sanity: P_{11}=e_2, P_{111}=e_3 ; Schur at q=T
print('P11',Pp[(1,1)]); print('P2 at s=1/t ->', {k:sp.factor(v.subs(s,1/t)) for k,v in Pp[(2,)].items()})
def w(g,nu,a,b): return g(a)*g(b)/g(nu)
gN=lambda l: t**nf(l)*s**nf(conj(l)); gC=lambda l: t**sum(l)
def lim(expr): return sp.limit(sp.simplify(expr),t,sp.oo)
def rowadd(a,b,nu):
    m=max(len(a),len(b)); a2=list(a)+[0]*(m-len(a)); b2=list(b)+[0]*(m-len(b))
    return tuple(sorted([x+y for x,y in zip(a2,b2)],reverse=True)) == nu and all(x+y>0 for x,y in zip(a2,b2)) and tuple(x+y for x,y in zip(a2,b2))==tuple(sorted([x+y for x,y in zip(a2,b2)],reverse=True))
def rowadd2(a,b,nu):
    m=max(len(a),len(b),1); a2=list(a)+[0]*(m-len(a)); b2=list(b)+[0]*(m-len(b))
    return tuple(x for x in (u+v for u,v in zip(a2,b2)) if x>0)==nu
# (i) witness
c=R[(1,1)].get(((1,),(1,)))
print('(i) coef P1⊗P1 in ΔP11 =',sp.factor(c),' N-weight', sp.factor(w(gN,(1,1),(1,),(1,))))
bad=[]; L={}
for nu,dP in R.items():
    for (a,b),c in dP.items():
        if a==() or b==(): continue
        v=sp.factor(lim(c*w(gN,nu,a,b)))
        ra=rowadd2(a,b,nu)
        if (v!=0)!=ra: bad.append((nu,a,b,v,ra))
        if v!=0: L[(nu,a,b)]=v
print('(ii-a) support == row-additive splits:', 'YES' if not bad else bad)
for k,v in sorted(L.items()): print('  ',k,v)
# row-additive splits all present?
