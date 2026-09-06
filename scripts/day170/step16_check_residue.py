"""Day 170 Step 16 — Check if 'residue' vanishes on the variety Y = T*phi(Y), q = 1-sT-2pTY."""
import sympy as sp

T, s, p, Y, q = sp.symbols('T s p Y q')

# Residue found in step 15 (using dE2_q = -2T^2/q)
resid = (48*T**6*Y*p**4 - 100*T**6*Y*p**3*s**2 + 74*T**6*Y*p**2*s**4 - 21*T**6*Y*p*s**6 + 2*T**6*Y*s**8 - 40*T**6*p**3*s + 54*T**6*p**2*s**3 - 19*T**6*p*s**5 + 2*T**6*s**7 + 32*T**5*Y*p**3*q*s + 136*T**5*Y*p**3*s - 44*T**5*Y*p**2*q*s**3 - 188*T**5*Y*p**2*s**3 + 17*T**5*Y*p*q*s**5 + 87*T**5*Y*p*s**5 - 2*T**5*Y*q*s**7 - 12*T**5*Y*s**7 + 24*T**5*p**3*q + 56*T**5*p**3 - 34*T**5*p**2*q*s**2 - 90*T**5*p**2*s**2 + 15*T**5*p*q*s**4 + 60*T**5*p*s**4 - 2*T**5*q*s**6 - 10*T**5*s**6 - 40*T**4*Y*p**3*q - 168*T**4*Y*p**3 + 78*T**4*Y*p**2*q*s**2 + 230*T**4*Y*p**2*s**2 - 54*T**4*Y*p*q*s**4 - 149*T**4*Y*p*s**4 + 10*T**4*Y*q*s**6 + 30*T**4*Y*s**6 + 32*T**4*p**2*q*s + 80*T**4*p**2*s - 33*T**4*p*q*s**3 - 77*T**4*p*s**3 + 8*T**4*q*s**5 + 20*T**4*s**5 - 88*T**3*Y*p**2*q*s - 192*T**3*Y*p**2*s + 71*T**3*Y*p*q*s**3 + 146*T**3*Y*p*s**3 - 20*T**3*Y*q*s**5 - 40*T**3*Y*s**5 - 64*T**3*p**2*q - 44*T**3*p**2 + 32*T**3*p*q*s**2 + 61*T**3*p*s**2 - 12*T**3*q*s**4 - 20*T**3*s**4 + 54*T**2*Y*p**2*q + 76*T**2*Y*p**2 - 59*T**2*Y*p*q*s**2 - 99*T**2*Y*p*s**2 + 20*T**2*Y*q*s**4 + 30*T**2*Y*s**4 - 25*T**2*p*q*s - 36*T**2*p*s + 8*T**2*q*s**3 + 10*T**2*s**3 + 36*T*Y*p*q*s + 47*T*Y*p*s - 10*T*Y*q*s**3 - 12*T*Y*s**3 + 11*T*p*q + 11*T*p - 2*T*q*s**2 - 2*T*s**2 - 11*Y*p*q - 11*Y*p + 2*Y*q*s**2 + 2*Y*s**2)

# Substitute q = 1 - sT - 2pTY
resid_q = resid.subs(q, 1 - s*T - 2*p*T*Y)
resid_q = sp.expand(resid_q)

# Now poly in Y with T, s, p params. Reduce Y^k for k >= 2.
Y_rel_val = ((1-s*T)*Y - T) / (p*T)   # Y^2 = ...

def reduce_Y(expr, max_iter=60):
    e = sp.expand(expr)
    for _ in range(max_iter):
        e_new = sp.expand(e.subs(Y**2, Y_rel_val))
        if e_new == e: break
        e = e_new
    return e

resid_reduced = reduce_Y(resid_q)
print(f"Residue (after q substituted, Y^2 reduced):")
if resid_reduced == 0:
    print("  = 0  --- IDENTITY HOLDS!")
else:
    print(f"  = {resid_reduced}")
    print(f"  factored: {sp.factor(resid_reduced)}")
    poly = sp.Poly(resid_reduced, Y)
    print(f"  Y-degree: {poly.degree()}")
    print(f"  As poly in Y:")
    for (dY,), c in sorted(poly.terms()):
        print(f"    Y^{dY}: {sp.factor(c)}")
