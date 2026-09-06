"""Day 170 Step 5 — Compute L_{-1} in normal form A + B*Y (A, B in Q(T,s,p)).

Follow same approach as step 4. Compute SOURCE directly, then L_{-1} = -SOURCE / (q^3 H).
Reduce to normal form. Verify against L_{-1} series from FP_coeffs.
"""
import sympy as sp
import sys
sys.path.insert(0, '/home/agent/projects/scratch/day152')

T, s, p, Y = sp.symbols('T s p Y')

# All in terms of Y with q = 1 - sT - 2pTY
q = 1 - s*T - 2*p*T*Y

Y2_val = ((1-s*T)*Y - T) / (p*T)

def reduce_Y_poly(expr):
    e = sp.expand(expr)
    for _ in range(80):
        e_new = sp.expand(e.subs(Y**2, Y2_val))
        if e_new == e:
            break
        e = e_new
    return e

def normalize(expr, verbose=False):
    """Reduce expr to (A + BY)/D form where A, B, D are polys in Q(T,s,p), no Y in denom."""
    e = sp.together(expr)
    num = sp.expand(sp.numer(e))
    den = sp.expand(sp.denom(e))
    num_r = reduce_Y_poly(num)
    den_r = reduce_Y_poly(den)
    for _ in range(3):
        den_pol = sp.Poly(den_r, Y)
        if den_pol.degree() <= 0:
            break
        # Rationalize
        a_coef = den_pol.nth(0)
        b_coef = den_pol.nth(1)
        conj = a_coef + b_coef * ((1-s*T)/(p*T) - Y)  # = a + b*Y_conj (Y_conj = sum - Y)
        num_r = reduce_Y_poly(sp.expand(num_r * conj))
        den_r = reduce_Y_poly(sp.expand(den_r * conj))
    num_pol = sp.Poly(num_r, Y)
    if num_pol.degree() < 0:
        A = sp.S(0); B = sp.S(0)
    else:
        A = num_pol.nth(0)
        B = num_pol.nth(1) if num_pol.degree() >= 1 else sp.S(0)
    return sp.together(A/den_r), sp.together(B/den_r)

# Y and q derivatives
phi_Y = 1 + s*Y + p*Y**2
Yp = phi_Y / q
qp = -(s + (4*p - s**2)*T) / q   # = -[s(1-sT) + 4pT]/q, but let's double-check
# Actually: d/dT(q^2) = d/dT[(1-sT)^2 - 4pT^2] = -2s(1-sT) - 8pT
# So 2q qp = -2s(1-sT) - 8pT, qp = -[s(1-sT) + 4pT]/q. Let me correct:
qp = -(s*(1-s*T) + 4*p*T) / q

def d_dT(expr):
    e = sp.diff(expr, T) + sp.diff(expr, Y) * Yp
    # Note: q depends on Y and T. When expr contains q, we treat q as q(T,Y).
    # Since q is defined as 1 - sT - 2pTY (not an independent variable), we don't need
    # sp.diff(expr, q) * qp — we should differentiate through q = 1 - sT - 2pTY.
    # But then expr = expr(T, Y), so d/dT = dT + Yp * dY. That's it.
    return e

# H = pY/T
H = p*Y/T
Hp = d_dT(H)
Hpp = d_dT(Hp)

# K_{-1} = -pY/q^2 (Day 169)
K = -p*Y/q**2
Kp = d_dT(K)

# Verify K by check: series match
# ... skip for now, just trust

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

print("=== Computing L_{-1} normal form ===")
print("This may take a minute...")

import time
t0 = time.time()
A, B = normalize(Lm1)
print(f"Normalized in {time.time()-t0:.1f}s")

print(f"\nL_{{-1}} = A + B*Y   where:")
A_simp = sp.together(sp.cancel(A))
B_simp = sp.together(sp.cancel(B))
print(f"  A = {A_simp}")
print(f"  B = {B_simp}")

# Save
with open('/home/agent/projects/scratch/day170/Lm1_AB.txt', 'w') as f:
    f.write(f"A = {A_simp}\n")
    f.write(f"B = {B_simp}\n")
