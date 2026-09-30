"""Day 212: Step B composite. For concrete I, exact s,t,z: with 207b's closed forms
 C_i(Y)/(D_i Y^i G(t^2Y)) = (1-stY)(1-sY/A)/((1-tY)(1-Y)), C_i(tY)/(...) = (A-stY)/(1-tY),
 C_{i-1}(tY)/(D_i Y^i G(t^2Y)) = (1-A)/(ts-A)(A/t)Y^{-1}(1-st^2Y/A)/(1-tY),
check: [kappa'_c x/(t^{i_c-1}z_c - ctx) * pref(I-e_c)/pref(I) * K_{I-e_c}/K_I * normalised C's] == S_c/prod(1-tX)
and L0 formula, as rational functions of x."""
import sympy as sp, itertools
from fractions import Fraction
x=sp.Symbol('x'); ok=True
for (s,t,Zall) in [(sp.Rational(3,7),sp.Rational(5,11),[sp.Rational(2,3),sp.Rational(9,8),sp.Rational(13,5),sp.Rational(17,19)])]:
  for ell in (2,3,4):
    Z=Zall[:ell]
    def Kp(i,j,z,w): return sp.prod([(t**p*z-s*w)/(t**p*z-t**j*w) for p in range(i)])*sp.prod([(s*z-t**r*w)/(t**i*z-t**r*w) for r in range(j)])
    def KI(I): return sp.prod([Kp(I[a],I[b],Z[a],Z[b]) for a in range(ell) for b in range(a+1,ell)])
    def kap(I,c):
        g=[t**I[d]*Z[d] for d in range(ell)]
        return sp.prod([g[c]-s*zd for zd in Z])/(g[c]*sp.prod([g[c]-g[d] for d in range(ell) if d!=c]))
    for I in itertools.product(range(3),repeat=ell):
        S=sum(I); A=[t**i for i in I]; cc=s**ell*t**(-S)
        X=[cc*x/(s*Z[c]) for c in range(ell)]; u=[s*t*Xc for Xc in X]
        pref=lambda J: s**(-(ell-1)*sum(J))*t**(sum(J)**2-sum(j*j for j in J))
        U=sp.prod([(1-X[c])/(1-s*X[c]/A[c]) for c in range(ell)])
        Ch=[(1-s*t*X[c])*(1-s*X[c]/A[c])/((1-t*X[c])*(1-X[c])) for c in range(ell)]
        Cht=[(A[c]-u[c])/(1-t*X[c]) for c in range(ell)]
        L0=sp.prod(Ch)*U-sp.prod(Cht)
        good=sp.cancel(L0-(sp.prod([1-uc for uc in u])-sp.prod([A[c]-u[c] for c in range(ell)]))/sp.prod([1-t*Xc for Xc in X]))==0
        tot=L0
        for c in range(ell):
            if I[c]==0: continue
            J=list(I); J[c]-=1; J=tuple(J)
            cols=sp.prod([Cht[d] for d in range(ell) if d!=c])*(1-A[c])/(t*s-A[c])*(A[c]/t)/X[c]*(1-s*t**2*X[c]/A[c])/(1-t*X[c])
            term=kap(J,c)*x/(t**(I[c]-1)*Z[c]-cc*t*x)*pref(J)/pref(I)*KI(J)/KI(I)*cols
            Psi=-sp.prod([(A[c]*X[d]-X[c])/(A[c]*X[d]-A[d]*X[c]) for d in range(ell) if d!=c])
            Sc=(1-A[c])*Psi*sp.prod([A[d]-u[d] for d in range(ell) if d!=c])/sp.prod([1-t*Xc for Xc in X])
            good&=sp.cancel(term-Sc)==0
            tot+=term
        good&=sp.cancel(tot)==0
        ok&=good
        if not good: print('FAIL',ell,I)
    print(f'ell={ell}: Step B composite + L0 + total=0 over I in [0,2]^ell:', 'OK' if ok else 'FAIL',flush=True)
print('ALL OK' if ok else 'FAIL')
