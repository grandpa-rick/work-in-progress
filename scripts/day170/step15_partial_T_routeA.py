"""Day 170 Step 15 — Compute partial_T Route_A as a rational function.

Derivation (using algebraic simplifications):

2 * Route_A = 2*d^2_E1 xi_0 + 3*E1*d_E1_E2 xi_0 + E1^2*d^2_E2 xi_0
              + 2*d_E1 log q + E1*d_E2 log q + d_E2 xi_0
              - T/q + T(q + R_1R_2)/q^3 - E_1*T*Y/q

Take partial_T:
2 * (partial_T Route_A) = 2*d^2_E1(d_T xi_0) + 3*E1*d_E1_E2(d_T xi_0) + E1^2*d^2_E2(d_T xi_0)
                          + 2*d_E1(q'/q) + E1*d_E2(q'/q) + d_E2(d_T xi_0)
                          - d_T(T/q) + d_T[T(q + R_1R_2)/q^3] - E_1*d_T(TY/q)

where d_T xi_0 = E_2 * Y/T = E_2*phi (rational in Y, s, p).

All E-partials computed via chain rule with:
  d_E1 Y = TY/q,  d_E2 Y = T Y^2/q,
  d_E1 q = -T(1-E1 T)/q,  d_E2 q = -TY(q+1-E1 T)/q.
"""
import sympy as sp
import time

T, s, p, Y, q = sp.symbols('T s p Y q')
E1, E2 = s, p

# Basic quantities
phi_Y = 1 + E1*Y + E2*Y**2  # = Y/T
Yp = phi_Y / q               # dY/dT (chain rule)
qp = -(E1*(1-E1*T) + 4*E2*T) / q  # dq/dT

# E-partials of Y and q
dE1_Y = T*Y/q
dE2_Y = T*Y**2/q
dE1_q = -T*(1-E1*T)/q
dE2_q = -2*T**2/q  # cleanest form from q^2 = (1-E1 T)^2 - 4 E2 T^2

# Verify: check dE1_q via q^2 relation: q^2 = (1-E1 T)^2 - 4 E2 T^2, so 2q * dE1_q = -2T(1-E1 T)
# -> dE1_q = -T(1-E1 T)/q ✓
# Verify dE2_q: 2q * dE2_q = -4 T^2, so dE2_q = -2T^2/q. But Day 162 says dE2_q = -TY(q+1-E1 T)/q.
# These should agree given the relation. Let me verify: q = 1 - E1 T - 2 E2 T Y, so
# dE2_q = -2 T Y - 2 E2 T dE2_Y = -2 T Y - 2 E2 T · T Y^2/q = -2TY - 2 E2 T^2 Y^2/q
#        = -2TY(1 + E2 T Y/q) = -2 T Y (q + E2 T Y)/q.
# And Day 162 says: dE2_q = -T Y (q + 1 - E1 T)/q = -T Y (q + q + 2 E2 T Y)/q [using 1-E1T=q+2E2TY]
#                        = -T Y · 2(q + E2 T Y)/q = -2 T Y (q + E2 T Y)/q. ✓
# So dE2_q = -T Y (q + 1 - E1 T)/q = -2 T Y (q + E2 T Y)/q. Both correct.

# General E-partial operator: for f(T, E1, E2, Y, q), take partial with
# implicit dependence on Y, q via E-partials.

def dE1(expr):
    return sp.diff(expr, E1) + sp.diff(expr, Y) * dE1_Y + sp.diff(expr, q) * dE1_q

def dE2(expr):
    return sp.diff(expr, E2) + sp.diff(expr, Y) * dE2_Y + sp.diff(expr, q) * dE2_q

# Since E1 = s and E2 = p, we need to be careful — for a raw expression in s, p, T, Y, q,
# we differentiate w.r.t. s or p and add Y, q chain rule terms.
# But sp.diff(expr, s) treats other symbols (including Y, q) as constants — correct.

# d_T xi_0 = E_2 * Y/T = p*Y/T = phi_Y * p
xi_0_T = E2 * Y / T  # = E_2 * phi (implicit in Y)

# Verify at series level
sv, pv = sp.Rational(2), sp.Rational(3)
N = 14
def Y_series(N):
    Ys = [sp.S(0)] * (N+1)
    for n in range(1, N+1):
        acc = sp.S(0)
        if n == 1: acc += 1
        acc += sv * Ys[n-1]
        for k in range(n):
            acc += pv * Ys[k] * Ys[n-1-k]
        Ys[n] = acc
    return Ys

Y_list = Y_series(N)
Y_ser_expr = sum(Y_list[n]*T**n for n in range(N+1))
q_ser_expr = sp.expand(1 - sv*T - 2*pv*T*Y_ser_expr)

# Compute 2*partial_T Route_A
print("Building partial_T Route_A...")
t0 = time.time()

# Terms:
# 2*d^2_E1(xi_0_T)
term1 = 2 * dE1(dE1(xi_0_T))

# 3*E1*d_E1_E2(xi_0_T)
term2 = 3 * E1 * dE1(dE2(xi_0_T))

# E1^2 * d^2_E2(xi_0_T)
term3 = E1**2 * dE2(dE2(xi_0_T))

# 2*d_E1(q'/q)
term4 = 2 * dE1(qp / q)  # q'/q = qp/q where qp = dq/dT

# E1*d_E2(q'/q)
term5 = E1 * dE2(qp / q)

# d_E2(xi_0_T)
term6 = dE2(xi_0_T)

# -d_T(T/q) = -(1/q - T*qp/q^2) = -1/q + T*qp/q^2
term7 = -(1/q - T * qp / q**2)  # = -1/q + T*qp/q^2

# d_T[T(q + R1R2)/q^3]
R1R2 = 1 - T**2 * (E1**2 - 4*E2)
tmp = T*(q + R1R2)/q**3
def d_dT(expr):
    return sp.diff(expr, T) + sp.diff(expr, Y) * Yp + sp.diff(expr, q) * qp
term8 = d_dT(tmp)

# -E_1 * d_T(TY/q)
term9 = -E1 * d_dT(T*Y/q)

two_pT_RouteA = term1 + term2 + term3 + term4 + term5 + term6 + term7 + term8 + term9

print(f"  built in {time.time()-t0:.1f}s")

# Simplify
print("Simplifying...")
t0 = time.time()
two_pT_RouteA_can = sp.cancel(two_pT_RouteA)
print(f"  cancel: {time.time()-t0:.1f}s")

# Now compare 2*(partial_T Route_A) with 2*(LHS_deriv) = 2*(partial_T R^{(-1)} + Delta)
# Load LHS_deriv from step 14... actually recompute here.

# L_0
K0 = (p*Y*(2*q + 1) + s*q) / q**2
theta_K0 = T * d_dT(K0)
L0 = (1 + 3*T*K0 + T**2 * K0**2 + T * theta_K0) / q

# L_{-1}
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
Delta = Lm1 - L0

# R^{(-1)} Day 162
Rmn = T * (E2 * Y**2 * ((q+1)**2 - E1*T) + (q + R1R2)/2) / q**3
Rmn_prime = d_dT(Rmn)

two_LHS = 2 * (Rmn_prime + Delta)

# Difference
diff = two_pT_RouteA - two_LHS
print("Cancel diff...")
t0 = time.time()
diff_can = sp.cancel(diff)
print(f"  cancel: {time.time()-t0:.1f}s")

# Reduce in ring
Y_rel = p*T*Y**2 - (1-s*T)*Y + T
q_rel = q**2 - ((1-s*T)**2 - 4*p*T**2)
gb = sp.groebner([Y_rel, q_rel], Y, q, order='lex', domain=sp.QQ.frac_field(T, s, p))

num_d = sp.expand(sp.numer(diff_can))
den_d = sp.expand(sp.denom(diff_can))

print(f"num degree: Y={sp.Poly(num_d, Y).degree() if Y in num_d.free_symbols else 0}, q={sp.Poly(num_d, q).degree() if q in num_d.free_symbols else 0}")

num_red = gb.reduce(sp.Poly(num_d, Y, q, domain=sp.QQ.frac_field(T, s, p)))[1].as_expr()
den_red = gb.reduce(sp.Poly(den_d, Y, q, domain=sp.QQ.frac_field(T, s, p)))[1].as_expr()

diff_reduced = sp.cancel(num_red / den_red)
print(f"\ndiff reduced = {diff_reduced}")
if diff_reduced == 0:
    print("\n*** IDENTITY PROVED! ***")
else:
    print("\nDiff nonzero — investigating residue...")
    diff_num = sp.numer(diff_reduced)
    diff_den = sp.denom(diff_reduced)
    diff_num_pol = sp.Poly(sp.expand(diff_num), Y, q)
    print("residue numerator as poly in Y, q:")
    for (dY, dq), c in sorted(diff_num_pol.terms()):
        print(f"  Y^{dY} q^{dq}: {sp.factor(c)}")
    print(f"denom: {sp.factor(sp.expand(diff_den))}")

# Also do series check
print("\n=== Series check: 2*d_T Route_A vs 2*LHS ===")
sub = {s: sv, p: pv, Y: Y_ser_expr, q: q_ser_expr}
two_pT_RA_ser = sp.series(two_pT_RouteA.subs(sub), T, 0, N+1).removeO()
two_LHS_ser = sp.series(two_LHS.subs(sub), T, 0, N+1).removeO()
for n in range(N+1):
    a = sp.expand(two_pT_RA_ser.coeff(T, n))
    b = sp.expand(two_LHS_ser.coeff(T, n))
    print(f"  [T^{n}]: 2 d_T RA = {a}, 2 LHS = {b}, diff = {a-b}")
