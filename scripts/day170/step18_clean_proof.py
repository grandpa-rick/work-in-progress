"""Day 170 Step 18 — Clean proof by explicit q substitution + Y-reduction.

The RING is Q(T, s, p)[Y] / (pTY^2 + (sT-1)Y + T), with q = 1 - sT - 2pTY.

We prove partial_T Route_A - partial_T R^{(-1)} - (L_{-1} - L_0) = 0 in this ring.
"""
import sympy as sp
import time

T, s, p, Y, q = sp.symbols('T s p Y q')
E1, E2 = s, p

# Derivative rules (with q as independent symbol, tracked via qp = ..., dE1_q, dE2_q)
phi_Y = 1 + E1*Y + E2*Y**2
Yp = phi_Y / q
qp = -(E1*(1-E1*T) + 4*E2*T) / q

dE1_Y = T*Y/q
dE2_Y = T*Y**2/q
dE1_q = -T*(1-E1*T)/q
dE2_q = -2*T**2/q

def dE1(expr):
    return sp.diff(expr, E1) + sp.diff(expr, Y) * dE1_Y + sp.diff(expr, q) * dE1_q
def dE2(expr):
    return sp.diff(expr, E2) + sp.diff(expr, Y) * dE2_Y + sp.diff(expr, q) * dE2_q
def d_dT(expr):
    return sp.diff(expr, T) + sp.diff(expr, Y) * Yp + sp.diff(expr, q) * qp

# ---- Assemble 2 * partial_T Route_A ----
xi_0_T = E2 * Y / T
R1R2 = 1 - T**2 * (E1**2 - 4*E2)
term1 = 2 * dE1(dE1(xi_0_T))
term2 = 3 * E1 * dE1(dE2(xi_0_T))
term3 = E1**2 * dE2(dE2(xi_0_T))
term4 = 2 * dE1(qp / q)
term5 = E1 * dE2(qp / q)
term6 = dE2(xi_0_T)
term7 = -(1/q - T * qp / q**2)
term8 = d_dT(T*(q + R1R2)/q**3)
term9 = -E1 * d_dT(T*Y/q)
two_pT_RouteA = term1 + term2 + term3 + term4 + term5 + term6 + term7 + term8 + term9

# ---- Assemble 2 * (partial_T R^{(-1)} + (L_{-1} - L_0)) ----
K0 = (p*Y*(2*q + 1) + s*q) / q**2
theta_K0 = T * d_dT(K0)
L0 = (1 + 3*T*K0 + T**2 * K0**2 + T * theta_K0) / q

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
    R3 * Hpp + coef_Hp * Hp + coef_H * H
    + 3 * R3 * (H * Kp + K * Hp)
    + coef_H2 * H**2 + 18 * T**3 * H * Hp + T**4 * H**3
    + 3 * R3 * H * K**2
    + R2 * Kp + coef_Hp * 2 * H * K + coef_K * K + R2 * K**2
    + 18 * T**3 * H**2 * K
)
Lm1 = -SOURCE / (q**3 * H)
Rmn = T * (E2 * Y**2 * ((q+1)**2 - E1*T) + (q + R1R2)/2) / q**3

two_LHS = 2 * (d_dT(Rmn) + Lm1 - L0)

diff = two_pT_RouteA - two_LHS

# ---- Reduce in ring: substitute q = 1 - sT - 2pTY, reduce Y^2 ----
print("Cancel diff...")
t0 = time.time()
diff_can = sp.cancel(diff)
print(f"  cancel: {time.time()-t0:.1f}s")

print("\nSubstitute q = 1 - sT - 2pTY...")
t0 = time.time()
q_of_Y = 1 - s*T - 2*p*T*Y
diff_q = sp.expand(diff_can.subs(q, q_of_Y))
print(f"  subs: {time.time()-t0:.1f}s")

# Extract numerator and denominator
diff_q_together = sp.together(diff_q)
num = sp.expand(sp.numer(diff_q_together))
den = sp.expand(sp.denom(diff_q_together))

# Now reduce Y^k for k >= 2 using pTY^2 = (1-sT)Y - T
Y_rel = p*T*Y**2 - (1-s*T)*Y + T

# Polynomial division: express num, den as polynomials in Y and take remainder mod Y_rel
def reduce_by_Y_rel(expr):
    pol = sp.Poly(expr, Y)
    rel = sp.Poly(Y_rel, Y)
    _, rem = sp.div(pol, rel, Y)
    return sp.expand(rem.as_expr())

print("\nReduce num by Y-relation...")
t0 = time.time()
num_red = reduce_by_Y_rel(num)
print(f"  reduce num: {time.time()-t0:.1f}s")

print("Reduce den by Y-relation...")
t0 = time.time()
den_red = reduce_by_Y_rel(den)
print(f"  reduce den: {time.time()-t0:.1f}s")

print(f"\nnum reduced (poly in Y, deg <= 1): {num_red}")
print(f"den reduced (poly in Y, deg <= 1): {den_red}")

# The ring is Q(T,s,p)[Y]/(pTY^2 + (sT-1)Y + T), which is rank 2 over Q(T,s,p).
# Element is 0 iff its normal form is 0.
# So num_red should be 0 (if identity holds), regardless of den_red.

if sp.expand(num_red) == 0:
    print("\n" + "=" * 60)
    print("*** IDENTITY PROVED ALGEBRAICALLY ***")
    print("=" * 60)
    print("The difference partial_T Route_A - [partial_T R^{(-1)} + (L_{-1} - L_0)]")
    print("reduces to 0 in the ring Q(T, s, p)[Y] / (pTY^2 + (sT-1)Y + T)")
    print("under q = 1 - sT - 2pTY.")
else:
    print("\nUnexpectedly nonzero. Investigate further.")
    print(f"factored num_red: {sp.factor(num_red)}")
