"""Day 210 Task 1: sober re-check of connection-note §1: U := x F(cx) = prod_c (1-X_c)/(1-s X_c/A_c), for ell=1,2,3,4.
Definitions (Day 209 §2, ell-column extension):
  H(y;x) = y prod_c(1+s y z_c) / (prod_c(1+gam_c y)(1+x y)),  gam_c = A_c z_c, A_c = t^{i_c} (A_c free symbols here)
  c := lim_{y->oo} x*H(y;x)   (the infinity constant: H/y ~ c/(x y))
  F(y) = prod_c (y - s z_c)/(y prod_c (y - gam_c));  check Res_0 F == c, and partial fractions F = c/y + sum kappa_c/(y-gam_c)
  kappa_c := Res_{y=gam_c} F ;  U := 1 - sum_c kappa_c x/(gam_c - c x)  (the diagonal factor, as in (star2-GF))
  X_c := c x/(s z_c).  Claim: U == x F(cx) == prod_c (1-X_c)/(1 - s X_c/A_c).
"""
import sympy as sp
s, t, x, y = sp.symbols('s t x y')
ok = True
for ell in (1, 2, 3, 4):
    z = sp.symbols(f'z1:{ell+1}'); A = sp.symbols(f'A1:{ell+1}')
    gam = [A[c]*z[c] for c in range(ell)]
    H = y*sp.Mul(*[1+s*y*zc for zc in z])/(sp.Mul(*[1+g*y for g in gam])*(1+x*y))
    cinf = sp.simplify(sp.limit(x*H, y, sp.oo))
    F = sp.Mul(*[y - s*zc for zc in z])/(y*sp.Mul(*[y-g for g in gam]))
    res0 = sp.simplify(sp.residue(F, y, 0)) if ell < 4 else sp.cancel((F*y).subs(y, 0))
    kap = [sp.cancel((F*(y-g)).subs(y, g)) for g in gam]
    pf = sp.cancel(F - cinf/y - sum(kap[c]/(y-gam[c]) for c in range(ell)))
    U = 1 - sum(kap[c]*x/(gam[c]-cinf*x) for c in range(ell))
    Xc = [cinf*x/(s*zc) for zc in z]
    prod = sp.Mul(*[(1-Xc[c])/(1-s*Xc[c]/A[c]) for c in range(ell)])
    d1 = sp.cancel(U - x*F.subs(y, cinf*x)); d2 = sp.cancel(U - prod)
    good = sp.cancel(cinf - res0) == 0 and pf == 0 and d1 == 0 and d2 == 0
    ok &= good
    print(f'ell={ell}: c = {cinf};  Res0F==c: {sp.cancel(cinf-res0)==0};  F=c/y+sum kap/(y-gam): {pf==0};  U==xF(cx): {d1==0};  U==prod: {d2==0}')
print('ALL OK' if ok else 'FAIL')
