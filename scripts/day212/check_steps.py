"""Day 212: check each displayed step of the (★ℓ) proof.
(P1) partial fractions of F; (P2) U = xF(cx); (P3) U product; (P4) x/(t^{i-1}z-ctx);
(P5),(P6),(★ℓ coeff) moved to check_concrete.py (exact rational points);
(Z) closing identity, free lambda,mu, ell<=5;  (★ℓ-GF) normalised with free A_c, ell<=4."""
import sympy as sp, itertools, sys
from step_map_ell import Kp, cf, N, V, poch, al
s,t,x,w=sp.symbols('s t x w')
ok=True
def rep(name,val):
    global ok
    print(name, 'OK' if val else 'FAIL', flush=True); ok&=bool(val)
# ---- (Z)
for ell in range(1,6):
    lam=sp.symbols(f'l1:{ell+1}'); mu=sp.symbols(f'm1:{ell+1}')
    L=sp.prod([m-1 for m in mu])-sp.prod([l-1 for l in lam])
    R=sum((mu[c]-lam[c])*sp.prod([(lam[d]-1)*(lam[c]-mu[d])/(lam[c]-lam[d]) for d in range(ell) if d!=c]) for c in range(ell))
    rep(f'(Z) ell={ell}', sp.cancel(sp.together(L-R))==0)
# ---- normalised (★ℓ-GF) with free A: L0 + sum S_c = 0
for ell in range(1,5):
    A=sp.symbols(f'A1:{ell+1}'); X=sp.symbols(f'X1:{ell+1}')
    u=[s*t*Xc for Xc in X]
    L0=(sp.prod([1-uc for uc in u])-sp.prod([A[c]-u[c] for c in range(ell)]))
    Sc=sum((1-A[c])*(-sp.prod([(A[c]*X[d]-X[c])/(A[c]*X[d]-A[d]*X[c]) for d in range(ell) if d!=c]))*sp.prod([A[d]-u[d] for d in range(ell) if d!=c]) for c in range(ell))
    rep(f'(★ℓ-GF normalised, free A) ell={ell}', sp.cancel(sp.together(L0+Sc))==0)
# ---- (P1)-(P3) free gam, z, A
for ell in range(1,5):
    Z=sp.symbols(f'z1:{ell+1}'); A=sp.symbols(f'A1:{ell+1}')
    gam=[A[c]*Z[c] for c in range(ell)]
    cc=s**ell/sp.prod(A); y=sp.Symbol('y')
    F=sp.prod([y-s*zc for zc in Z])/(y*sp.prod([y-g for g in gam]))
    kap=[sp.prod([gam[c]-s*zd for zd in Z])/(gam[c]*sp.prod([gam[c]-gam[d] for d in range(ell) if d!=c])) for c in range(ell)]
    rep(f'(P1) ell={ell}', sp.cancel(F-cc/y-sum(kap[c]/(y-gam[c]) for c in range(ell)))==0)
    U=1-sum(kap[c]*x/(gam[c]-cc*x) for c in range(ell))
    rep(f'(P2) ell={ell}', sp.cancel(U-x*F.subs(y,cc*x))==0)
    Xc=[cc*x/(s*Z[c]) for c in range(ell)]
    rep(f'(P3) ell={ell}', sp.cancel(U-sp.prod([(1-Xc[c])/(1-s*Xc[c]/A[c]) for c in range(ell)]))==0)
    # (P4): x/(t^{i_c-1} z_c - c t x) = s t X_c/(c(A_c - s t^2 X_c))
    rep(f'(P4) ell={ell}', all(sp.cancel(x/(A[c]/t*Z[c]-cc*t*x)-s*t*Xc[c]/(cc*(A[c]-s*t**2*Xc[c])))==0 for c in range(ell)))
print('ALL OK' if ok else 'FAIL')
