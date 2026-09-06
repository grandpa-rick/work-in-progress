"""Day 170 Step 17 — Final proof check: the identity
    partial_T (Route_A) = partial_T (R^{(-1)}) + (L_{-1} - L_0)
holds as an element of the ring Q(T,s,p)[Y]/(pTY^2 + (sT-1)Y + T),
with q = 1 - sT - 2pTY.

Both sides vanish at T=0. Therefore Route_A = R^{(-1)} + I_Delta.

Since Prop 3 (Day 167, proved) states R^{(-1)} = Route_A - I_Delta,
Day 162's closed form R^{(-1)} = T[E_2*Y^2*((q+1)^2 - E_1*T) + (q+R_1R_2)/2]/q^3 is PROVED.
By Day 165 three-way equivalence: Sigma_0 closed form, Theorem B, R^{(-1)} closed form
all PROVED.
By Day 162 Theorem C: C.5 (equivalently Day 161 Theorem 4) PROVED.
By Day 167 Missing Lemma (R) chain: Missing Lemma (R) CLOSED.
"""
import sympy as sp
import time

T, s, p, Y, q = sp.symbols('T s p Y q')
E1, E2 = s, p

# Derivative rules
phi_Y = 1 + E1*Y + E2*Y**2
Yp = phi_Y / q
qp = -(E1*(1-E1*T) + 4*E2*T) / q

# E-partials
dE1_Y = T*Y/q
dE2_Y = T*Y**2/q
dE1_q = -T*(1-E1*T)/q
dE2_q = -2*T**2/q          # from q^2-relation directly

def dE1(expr):
    return sp.diff(expr, E1) + sp.diff(expr, Y) * dE1_Y + sp.diff(expr, q) * dE1_q
def dE2(expr):
    return sp.diff(expr, E2) + sp.diff(expr, Y) * dE2_Y + sp.diff(expr, q) * dE2_q
def d_dT(expr):
    return sp.diff(expr, T) + sp.diff(expr, Y) * Yp + sp.diff(expr, q) * qp

# ---- LHS: partial_T Route_A (via reformulated closed form) ----
# 2 * Route_A = 2*d^2_E1 xi_0 + 3*E1*d_E1_E2 xi_0 + E1^2*d^2_E2 xi_0
#              + 2*d_E1 log q + E1*d_E2 log q + d_E2 xi_0
#              - T/q + T(q + R_1R_2)/q^3 - E_1*T*Y/q
# where d_T xi_0 = E_2 * Y/T (so we differentiate the above w.r.t. T)

xi_0_T = E2 * Y / T   # partial_T xi_0

# Terms in 2 * partial_T Route_A (from differentiating in T the reformulated Route_A):
R1R2 = 1 - T**2 * (E1**2 - 4*E2)
term1 = 2 * dE1(dE1(xi_0_T))
term2 = 3 * E1 * dE1(dE2(xi_0_T))
term3 = E1**2 * dE2(dE2(xi_0_T))
term4 = 2 * dE1(qp / q)               # d_E1 (partial_T log q)
term5 = E1 * dE2(qp / q)              # E_1 * d_E2 (partial_T log q)
term6 = dE2(xi_0_T)
term7 = -(1/q - T * qp / q**2)         # -d_T(T/q)
term8 = d_dT(T*(q + R1R2)/q**3)
term9 = -E1 * d_dT(T*Y/q)

two_pT_RouteA = term1 + term2 + term3 + term4 + term5 + term6 + term7 + term8 + term9

# ---- RHS: 2 * (partial_T R^{(-1)} + (L_{-1} - L_0)) ----
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
    + 18 * T**3 * H**2 * K       # CORRECTED (P_3 G^3 e=1)
)
Lm1 = -SOURCE / (q**3 * H)
Delta = Lm1 - L0

Rmn = T * (E2 * Y**2 * ((q+1)**2 - E1*T) + (q + R1R2)/2) / q**3
two_LHS = 2 * (d_dT(Rmn) + Delta)

# ---- Check: two_pT_RouteA - two_LHS = 0 as element of ring ----
diff = two_pT_RouteA - two_LHS
print("Cancel diff...")
t0 = time.time()
diff_can = sp.cancel(diff)
print(f"  cancel: {time.time()-t0:.1f}s")

# Substitute q = 1 - sT - 2pTY (algebraic identity from linking) and reduce Y^2.
print("Substitute q = 1 - sT - 2pTY, reduce Y^2 in ring...")
diff_q = diff_can.subs(q, 1 - s*T - 2*p*T*Y)
num = sp.expand(sp.numer(diff_q))
den = sp.expand(sp.denom(diff_q))

Y2_val = ((1-s*T)*Y - T) / (p*T)

def reduce_Y(expr, max_iter=80):
    e = sp.expand(expr)
    for _ in range(max_iter):
        e_new = sp.expand(e.subs(Y**2, Y2_val))
        if e_new == e: break
        e = e_new
    return e

num_red = reduce_Y(num)
den_red = reduce_Y(den)

diff_final = sp.cancel(num_red / den_red)
print(f"\ndiff after q-substitution and Y^2 reduction: {diff_final}")
if sp.expand(diff_final) == 0:
    print("\n=== *** THEOREM B PROVED! *** ===")
    print("The identity holds as element of Q(T,s,p)[Y]/(pTY^2 + (sT-1)Y + T).")
    print("Combined with Route_A(0) = R^{(-1)}(0) = I_Delta(0) = 0, we conclude:")
    print("  R^{(-1)}(T) = Route_A(T) - I_Delta(T)")
    print("as formal power series in T with coefficients in Q(s, p).")
else:
    print("\nStill nonzero — investigate.")
    print(f"factor: {sp.factor(diff_final)}")

# Cross-check at series level for high N (independent verification)
print("\n=== Series cross-check at 3 different (s, p) ===")
for sv, pv in [(sp.Rational(2), sp.Rational(3)), (sp.Rational(1), sp.Rational(1)), (sp.Rational(5), sp.Rational(2))]:
    N = 10
    def Y_series(N, sv, pv):
        Ys = [sp.S(0)] * (N+1)
        for n in range(1, N+1):
            acc = sp.S(0)
            if n == 1: acc += 1
            acc += sv * Ys[n-1]
            for k in range(n):
                acc += pv * Ys[k] * Ys[n-1-k]
            Ys[n] = acc
        return Ys
    Ys = Y_series(N+2, sv, pv)
    Y_expr = sum(Ys[n]*T**n for n in range(N+2))
    q_expr = sp.expand(1 - sv*T - 2*pv*T*Y_expr)
    sub = {s: sv, p: pv, Y: Y_expr, q: q_expr}

    two_RA = sp.series(two_pT_RouteA.subs(sub), T, 0, N).removeO()
    two_L = sp.series(two_LHS.subs(sub), T, 0, N).removeO()

    ok_count = 0
    total = 0
    for n in range(N):
        a = sp.expand(two_RA.coeff(T, n))
        b = sp.expand(two_L.coeff(T, n))
        total += 1
        if sp.expand(a - b) == 0:
            ok_count += 1
    print(f"  (s, p) = ({sv}, {pv}): {ok_count}/{total} coefficients match")
