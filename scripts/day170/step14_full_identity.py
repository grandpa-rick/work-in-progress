"""Day 170 Step 14 — Full symbolic verification of Prop 3 identity via T-derivative.

Identity: R^{(-1)}(T) = Route_A(T) - I_Delta(T)
  where I_Delta = int_0^T (L_{-1} - L_0) dT'.

Differentiate: partial_T R^{(-1)} = partial_T Route_A - (L_{-1} - L_0).

At T=0 both sides are 0, so if the derivative identity holds, done.

STRATEGY: Show equality of RATIONAL FUNCTIONS in Q(s,p,T)[Y,q]/(Y-rel, q^2-rel).
- L_{-1} - L_0 computed directly (rational in Y, q, T, s, p).
- partial_T R^{(-1)}: derivative of Day 162 closed form R^{(-1)} = T[E_2 Y^2((q+1)^2 - E_1T) + (q + R_1R_2)/2]/q^3.
- partial_T Route_A: closed-form via chain rule + Day 158/161 machinery.

But Route_A involves E-partials of xi_0, xi_1, xi_2. Rather than derive it directly,
use the alternative reformulation:

Prop 3 says: R^{(-1)}_n = (1/2) partial_{u_3}^2 Xi_n |_0 - [deg=n-1][T^n] log(F_{-1}/F_0).

And Day 165 proves: given Sigma_0 closed form, R^{(-1)} = Day 162 closed form (via corrected L3).
So instead of verifying Route_A - I_Delta = R^{(-1)}, we can verify Sigma_0 closed form.

Sigma_0 = ℓ^{top}_0(L_A F_1 / F_0), Day 165 closed form:
   -Sigma_0 = (q + 1 - u)(q^2 - 6q + 6 - 6u) / (2 q^4)   where u = E_1 T

This is a rational function in q, T, E_1 (no Y!). Only 3 variables.

Actually looking again at the target: verify EITHER Sigma_0 OR the full identity.

Let me first attempt the FULL identity check symbolically in the ring.
"""
import sympy as sp
import time

T, s, p, Y, q = sp.symbols('T s p Y q')

# ============================================================
# Common setup: derivative rules
# ============================================================
phi_Y = 1 + s*Y + p*Y**2
Yp = phi_Y / q
qp = -(s*(1-s*T) + 4*p*T) / q

def d_dT(expr):
    return sp.diff(expr, T) + sp.diff(expr, Y) * Yp + sp.diff(expr, q) * qp

# ============================================================
# L_0 (Day 168) and L_{-1} (Day 169, corrected)
# ============================================================
K0 = (p*Y*(2*q + 1) + s*q) / q**2      # Day 158 K_0
theta_K0 = T * d_dT(K0)
L0 = (1 + 3*T*K0 + T**2 * K0**2 + T * theta_K0) / q  # Day 168

H = p*Y/T
Hp = d_dT(H)
Hpp = d_dT(Hp)
K = -p*Y/q**2                             # Day 169 K_{-1}
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
    + 18 * T**3 * H**2 * K       # CORRECTED (P_3 G^3 e=1)
)
Lm1 = -SOURCE / (q**3 * H)                # Day 169

Delta = Lm1 - L0                          # L_{-1} - L_0

# ============================================================
# R^{(-1)} closed form (Day 162): T*[E2*Y^2*((q+1)^2 - E1*T) + (q + R1R2)/2] / q^3
# R1R2 = 1 - T^2*(E1^2 - 4E2)
# ============================================================
E1, E2 = s, p                             # aliases
R1R2 = 1 - T**2 * (E1**2 - 4*E2)
Rmn = T * (E2 * Y**2 * ((q+1)**2 - E1*T) + (q + R1R2)/2) / q**3
Rmn_prime = d_dT(Rmn)

# ============================================================
# Route A T-derivative: partial_T Route_A = (1/2)*partial_{u_3}^2 (partial_T Xi)|_0
# partial_T Xi = ell^top_2(G) = H (the top layer of G).
# But H is defined for general u_3. We need its 2nd u_3-derivative at u_3=0.
#
# Alternative approach (RECOMMENDED): use direct series verification for now, and
# check whether ring-based algebraic verification is tractable.
# ============================================================

# For today, verify the KEY EQUATION at series level to high N and note the algebra:
# LHS (derivative form): partial_T R^{(-1)} + Delta = partial_T Route_A
# If we can show LHS = partial_T Route_A as rational functions, done.

# Compute LHS = partial_T R^{(-1)} + Delta as rational in Y, q, T, s, p
LHS_deriv = Rmn_prime + Delta

# Simplify LHS_deriv
print("Computing LHS_deriv = partial_T R^{(-1)} + (L_{-1} - L_0)...")
t0 = time.time()
LHS_deriv_cancel = sp.cancel(LHS_deriv)
print(f"  cancel: {time.time()-t0:.1f}s")
num_lhs = sp.expand(sp.numer(LHS_deriv_cancel))
den_lhs = sp.expand(sp.denom(LHS_deriv_cancel))
print(f"  num degree Y = {sp.Poly(num_lhs, Y).degree() if Y in num_lhs.free_symbols else 0}, degree q = {sp.Poly(num_lhs, q).degree() if q in num_lhs.free_symbols else 0}")
print(f"  den = {sp.factor(den_lhs)}")

# Reduce mod ring relations
Y_rel = p*T*Y**2 - (1-s*T)*Y + T
q_rel = q**2 - ((1-s*T)**2 - 4*p*T**2)
gb = sp.groebner([Y_rel, q_rel], Y, q, order='lex', domain=sp.QQ.frac_field(T, s, p))

num_lhs_red = gb.reduce(sp.Poly(num_lhs, Y, q, domain=sp.QQ.frac_field(T, s, p)))[1].as_expr()
den_lhs_red = gb.reduce(sp.Poly(den_lhs, Y, q, domain=sp.QQ.frac_field(T, s, p)))[1].as_expr()

LHS_reduced = sp.cancel(num_lhs_red / den_lhs_red)
print(f"\nLHS_deriv reduced in ring:")
print(f"  = {LHS_reduced}")

# Print poly form
LHS_num = sp.numer(LHS_reduced)
LHS_den = sp.denom(LHS_reduced)
LHS_num_pol = sp.Poly(sp.expand(LHS_num), Y, q)
print(f"\nAs poly in Y, q:")
for (dY, dq), c in sorted(LHS_num_pol.terms()):
    print(f"  Y^{dY} q^{dq}: {sp.factor(c)}")
print(f"Denom: {sp.factor(sp.expand(LHS_den))}")

# Save
with open('/home/agent/projects/scratch/day170/LHS_deriv_reduced.txt', 'w') as f:
    f.write(f"LHS_deriv (reduced) = {LHS_reduced}\n")
    f.write(f"Numerator poly in (Y, q):\n")
    for (dY, dq), c in sorted(LHS_num_pol.terms()):
        f.write(f"  Y^{dY} q^{dq}: {sp.factor(c)}\n")
    f.write(f"Denom: {sp.factor(sp.expand(LHS_den))}\n")
