"""Day 170 Step 6 — Cleaner L_{-1} normal form via Groebner reduction.

Strategy: compute L_{-1} as rational function in T, s, p, Y (with q substituted),
then reduce numerator and denominator mod Y^2 relation, then simplify.
"""
import sympy as sp
import time

T, s, p, Y = sp.symbols('T s p Y')

q = 1 - s*T - 2*p*T*Y

phi_Y = 1 + s*Y + p*Y**2
Yp = phi_Y / q

def d_dT(expr):
    """Total T-derivative, treating Y as function of T with Y' = phi/q."""
    return sp.diff(expr, T) + sp.diff(expr, Y) * Yp

# H = pY/T; K_{-1} = -pY/q^2
H = p*Y/T
Hp = d_dT(H)
Hpp = d_dT(Hp)

K = -p*Y/q**2
Kp = d_dT(K)

R3 = -T**2 * q**2
R2 = q**2 * (1 - s*T)
coef_Hp = -11*T + 14*s*T**2 + (12*p - 3*s**2) * T**3
coef_H  = 1 + 12*s*T + (5*p - s**2) * T**2
coef_H2 = 23*T**2 + s*T**3
coef_K  = -s + (2*s**2 + 10*p) * T + (4*p*s - s**3) * T**2

SOURCE = (
    R3 * Hpp
    + coef_Hp * Hp
    + coef_H * H
    + 3 * R3 * (H * Kp + K * Hp)
    + coef_H2 * H**2
    + 18 * T**3 * H * Hp
    + T**4 * H**3
    + 3 * R3 * H * K**2
    + R2 * Kp
    + coef_Hp * 2 * H * K
    + coef_K * K
    + R2 * K**2
)

Lm1 = -SOURCE / (q**3 * H)

# Step 1: expand fully as a rational function
print("Step 1: get rational form...")
t0 = time.time()
Lm1_rat = sp.cancel(Lm1)  # Cancel common factors, get single fraction
print(f"  done in {time.time()-t0:.1f}s")

num = sp.expand(sp.numer(Lm1_rat))
den = sp.expand(sp.denom(Lm1_rat))
print(f"  numerator degree in Y: {sp.Poly(num, Y).degree()}")
print(f"  denominator degree in Y: {sp.Poly(den, Y).degree()}")

# Step 2: reduce mod Y^2 = ((1-sT)Y - T)/(pT), i.e. pTY^2 - (1-sT)Y + T = 0
# Represent as polynomial in Y with coefficients in Q(T, s, p).
# Use polynomial division: divide by Y^2's replacement.

print("\nStep 2: reduce mod Y^2 relation...")

Y2_val = ((1-s*T)*Y - T) / (p*T)

def reduce_Y_poly(expr, max_iter=100):
    e = sp.expand(expr)
    for _ in range(max_iter):
        e_new = sp.expand(e.subs(Y**2, Y2_val))
        if e_new == e:
            break
        e = e_new
    return e

t0 = time.time()
num_red = reduce_Y_poly(num)
den_red = reduce_Y_poly(den)
print(f"  reduced in {time.time()-t0:.1f}s")
print(f"  num after reduction: degree in Y = {sp.Poly(num_red, Y).degree()}")
print(f"  den after reduction: degree in Y = {sp.Poly(den_red, Y).degree() if sp.Poly(den_red, Y).degree() > 0 else 0}")

# Step 3: if denominator still has Y, rationalize by multiplying by conjugate
print("\nStep 3: rationalize denominator (multiply by conjugate if needed)...")
den_poly = sp.Poly(den_red, Y)
if den_poly.degree() >= 1:
    a_coef = den_poly.nth(0)
    b_coef = den_poly.nth(1)
    # conjugate: replace Y with Y_conj = (sum) - Y, where sum = (1-sT)/(pT).
    # (a + bY)(a + b(sum - Y)) = (a + b·sum/2)^2 - (b·sum/2 - a)(...) hmm simpler:
    # (a + bY)(a + b*Y_conj) = a^2 + ab(Y + Y_conj) + b^2 Y Y_conj
    #                        = a^2 + ab·(1-sT)/(pT) + b^2·(1/p)   [product of roots = 1/p from Y^2 - (1-sT)/(pT) Y + 1/p = 0]
    # This gives a polynomial in T, s, p only.
    conj_multiplier = a_coef + b_coef * ((1-s*T)/(p*T) - Y)
    num_r2 = reduce_Y_poly(sp.expand(num_red * conj_multiplier))
    den_r2 = reduce_Y_poly(sp.expand(den_red * conj_multiplier))
    # Check den_r2 is Y-free
    assert sp.Poly(den_r2, Y).degree() == 0, f"Still Y: {den_r2}"
    num_red = num_r2
    den_red = den_r2

# Now den_red is polynomial in T, s, p only, and num_red is A + B*Y form.
den_final = sp.factor(den_red)
num_pol = sp.Poly(num_red, Y)
A = sp.together(num_pol.nth(0))
B = sp.together(num_pol.nth(1) if num_pol.degree() >= 1 else 0)

print("\n=== L_{-1} = (A + B*Y) / D ===")
print(f"D = {den_final}")
print(f"A = {A}")
print(f"\nB = {B}")

# Now try to factor A/D and B/D
A_over_D = sp.cancel(A / den_red)
B_over_D = sp.cancel(B / den_red)
print(f"\nA/D simplified:")
print(f"  numerator = {sp.factor(sp.expand(sp.numer(A_over_D)))}")
print(f"  denominator = {sp.factor(sp.expand(sp.denom(A_over_D)))}")
print(f"\nB/D simplified:")
print(f"  numerator = {sp.factor(sp.expand(sp.numer(B_over_D)))}")
print(f"  denominator = {sp.factor(sp.expand(sp.denom(B_over_D)))}")
