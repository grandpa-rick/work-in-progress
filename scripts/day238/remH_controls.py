# negative controls for remH_print_check: these MUST fail somewhere at n<=4
exec(open('remH_print_check.py').read().split('ok=True')[0])
fails_a=fails_c=0
for n in range(2,5):
    xs=sp.symbols('x1:%d'%(n+1)); L=parts(n)
    PT={l:HLP(l,xs,t) for l in L}; PTi={l:HLP(l,xs,1/t) for l in L}; S={l:HLP(l,xs,0) for l in L}
    for lam in L:
        el=sp.Mul(*[ek(k,xs) for k in lam]); Me=coeff_in_basis(el,PT,xs,L); Ks=coeff_in_basis(el,S,xs,L); Mi=coeff_in_basis(el,PTi,xs,L)
        for kap in L:  # control (a'): drop the conjugate on the Kostka-Foulkes index
            rhs=sum(Ks[conj(nu)]*coeff_in_basis(S[nu],PT,xs,L)[kap] for nu in L)
            if sp.simplify(Me[kap]-rhs)!=0: fails_a+=1
        d1=sp.factor(t**(nfun((1,)*n)-nfun(conj(lam)))*Mi[(1,)*n])  # control (c'): t-factorials in t^{-1}
        h=tfact(n,1/t)/sp.Mul(*[tfact(k,1/t) for k in lam])
        if sp.simplify(d1-h)!=0: fails_c+=1
print('control a broke',fails_a,'times; control c broke',fails_c,'times')
