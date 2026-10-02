import pickle, sympy as sp
s,t=sp.symbols('s t')
PSI=pickle.load(open('scripts/day216/psi_N5.pkl','rb'))
def nn(m): return sum(i*x for i,x in enumerate(m))
def lead(expr,var,inf=True):
    # leading exponent and coeff as var->inf (inf) or var->0
    num,den=sp.fraction(sp.cancel(expr))
    pn=sp.Poly(num,var); pd=sp.Poly(den,var)
    if inf:
        return pn.degree()-pd.degree(), sp.factor(pn.LC()/pd.LC())
    else:
        def low(p):
            d=min(m[0] for m in p.monoms()); return d, p.coeff_monomial(var**d)
        a,b=low(pn); c,d=low(pd); return a-c, sp.factor(b/d)
for mu,md in sorted(PSI.items(),key=lambda x:(sum(x[0]),x[0])):
    for name,var,inf in [('s->oo',s,True),('t->0',t,False)]:
        L={nu:lead(c,var,inf) for nu,c in md.items()}
        if inf: top=max(v[0] for v in L.values())
        else: top=min(v[0] for v in L.values())
        print(mu,name,'exp',top,{nu:v[1] for nu,v in L.items() if v[0]==top})
