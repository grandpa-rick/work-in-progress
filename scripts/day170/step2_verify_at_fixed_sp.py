"""Day 170 Step 2 — Verify identity at fixed (s, p) numerically to high N.

Substitute E_1 = 2, E_2 = 3 (or similar) and check the identity as pure T-series.
This is much faster than keeping s, p symbolic.
"""
import sympy as sp

N = 12

sv, pv = sp.Rational(2), sp.Rational(3)
T = sp.symbols('T')

# ----- Y as series in T with s=sv, p=pv -----
def make_Y_series(N):
    Y_series = [sp.S(0)] * (N + 2)
    for n in range(1, N + 2):
        s = sp.S(0)
        if n == 1:
            s += 1
        s += sv * Y_series[n-1]
        for k in range(n):
            s += pv * Y_series[k] * Y_series[n-1-k]
        Y_series[n] = s
    return Y_series

Y_ser = make_Y_series(N + 2)
Y_expr = sum(Y_ser[n] * T**n for n in range(N + 2))
Y_expr = sp.series(Y_expr, T, 0, N + 2).removeO()

q_sp = sp.expand(1 - sv*T - 2*T*pv*Y_expr)
q_sp = sp.series(q_sp, T, 0, N + 2).removeO()

phi_sp = 1 + sv*Y_expr + pv*Y_expr**2
log_q = sp.series(sp.log(q_sp), T, 0, N + 2).removeO()
log_YoverTq = sp.series(sp.log(phi_sp / q_sp), T, 0, N + 2).removeO()

# ----- Now we need xi_0, xi_1, xi_2 as series in T,  BUT keeping E1, E2 formal
# so we can take partial derivatives. But E1, E2 are s=2, p=3 numerically.
# So partial derivatives with respect to E1 or E2 don't make sense numerically.
#
# INSTEAD: compute Route A directly as [u_3^2] Xi (from Fact I).
# Route A = coefficient of u_3^2 in Xi, extracted from FP_coeffs.

# Load FP_coeffs
import sys
sys.path.insert(0, '/home/agent/projects/scratch/day152')
from lib import FP_coeffs, slog, pmul, padd, psub, pscal, pconst, ppow, sdiv, pzero
from fractions import Fraction as Fr

def substitute_u3(pcoef, u3_val, u1_val, u2_val):
    """Substitute u_1 = u1_val, u_2 = u2_val, u_3 = u3_val into a polynomial in u_1, u_2, u_3."""
    val = Fr(0)
    for (i, j, k), c in pcoef.items():
        val += c * (u1_val ** i) * (u2_val ** j) * (u3_val ** k)
    return val

FP = FP_coeffs(N + 2)
logFP = slog(FP, N + 2)

def extract_deg(pp, d):
    return {k: v for k, v in pp.items() if sum(k) == d}

# Xi_n = degree-(n+1) part of logFP[n]
Xi = [extract_deg(logFP[n], n + 1) for n in range(N + 2)]

# At u_3 = 0: (u_1, u_2) = (2, 3), so s = 5, p = 6? Wait but our sv, pv = (2, 3) as E_1, E_2.
# Hmm, we need to pick (u_1, u_2) such that E_1 = u_1 + u_2 = sv and E_2 = u_1 u_2 = pv.
# For sv=2, pv=3: roots of t^2 - 2t + 3 = 0 → t = 1 ± i sqrt(2). Complex — annoying but sympy OK.
# Alternative: pick (u_1, u_2) directly.

u1, u2 = 2, 3  # so E_1 = 5, E_2 = 6
# Recompute with sv=5, pv=6
sv_int, pv_int = 5, 6

# Actually let me redo with (sv, pv) = (5, 6) from the start.
# Simpler: use symbolic s, p and eval later. Restart with clean approach.

sv, pv = sp.Rational(5), sp.Rational(6)

Y_ser = [sp.S(0)] * (N + 2)
for n in range(1, N + 2):
    s = sp.S(0)
    if n == 1: s += 1
    s += sv * Y_ser[n-1]
    for k in range(n):
        s += pv * Y_ser[k] * Y_ser[n-1-k]
    Y_ser[n] = s
Y_expr = sum(Y_ser[n] * T**n for n in range(N + 2))

q_series = [sp.S(0)] * (N + 2)
q_series[0] = sp.S(1)
q_series[1] = -sv
# q^2 = 1 - 2sT + (s^2 - 4p) T^2, so q = sqrt(q^2). Use series.
q2 = [sp.S(0)] * (N + 2)
q2[0] = 1; q2[1] = -2*sv; q2[2] = sv**2 - 4*pv
for n in range(1, N + 2):
    s = q2[n] if n <= 2 else sp.S(0)
    for k in range(1, n):
        s -= q_series[k] * q_series[n-k]
    q_series[n] = s / 2
q_expr = sum(q_series[n] * T**n for n in range(N + 2))
q_expr = sp.expand(q_expr)

# Verify q = 1 - sT - 2pTY:
q_check = 1 - sv*T - 2*pv*T*Y_expr
# Check first several coefs
for n in range(min(6, N + 2)):
    a = sp.expand(q_expr).coeff(T, n)
    b = sp.expand(q_check).coeff(T, n)
    assert a == b, f"n={n}: {a} vs {b}"
print(f"q = 1 - sT - 2pTY matches for {sv_int, pv_int}")

# Now compute Route A numerically (as coefficient of u_3^2 in Xi, at u_1=2, u_2=3)
routeA_coefs = [None] * (N + 2)
for n in range(N + 2):
    total = Fr(0)
    for (i, j, k), c in Xi[n].items():
        if k == 2:
            total += c * (u1 ** i) * (u2 ** j)
    routeA_coefs[n] = total

# Compute R^{(-1)} directly (from Day 162 closed form)
R1R2 = 1 - T**2*(sv**2 - 4*pv)
inner = pv*Y_expr**2 * ((q_expr+1)**2 - sv*T) + sp.Rational(1,2)*(q_expr + R1R2)
Rmn = T * inner / q_expr**3
Rmn_ser = sp.series(Rmn, T, 0, N + 2).removeO()

# Compute L_0, L_{-1}, I_Delta
K0_expr = (pv*Y_expr*(2*q_expr+1) + sv*q_expr) / q_expr**2
K0_series = sp.series(K0_expr, T, 0, N + 2).removeO()
theta_K0 = sum(n * K0_series.coeff(T, n) * T**n for n in range(N + 2))
L0_num = 1 + 3*T*K0_series + T**2*K0_series**2 + T*theta_K0
L0_series = sp.series(L0_num / q_expr, T, 0, N + 2).removeO()

H_expr = pv * Y_expr / T
K_expr = -pv * Y_expr / q_expr**2

H_series = sp.series(H_expr, T, 0, N + 2).removeO()
Hp_series = sp.diff(H_series, T)
Hpp_series = sp.diff(Hp_series, T)
K_series = sp.series(K_expr, T, 0, N + 2).removeO()
Kp_series = sp.diff(K_series, T)

R3 = -T**2 * q_expr**2
R2 = q_expr**2 * (1 - sv*T)
coef_Hp = -11*T + 14*sv*T**2 + (12*pv - 3*sv**2) * T**3
coef_H = 1 + 12*sv*T + (5*pv - sv**2) * T**2
coef_H2 = 23*T**2 + sv*T**3
coef_K = -sv + (2*sv**2 + 10*pv) * T + (4*pv*sv - sv**3) * T**2

SOURCE = (
    R3 * Hpp_series
    + coef_Hp * Hp_series
    + coef_H * H_series
    + 3 * R3 * (H_series * Kp_series + K_series * Hp_series)
    + coef_H2 * H_series**2
    + 18 * T**3 * H_series * Hp_series
    + T**4 * H_series**3
    + 3 * R3 * H_series * K_series**2
    + R2 * Kp_series
    + coef_Hp * 2 * H_series * K_series
    + coef_K * K_series
    + R2 * K_series**2
)
SOURCE_ser = sp.series(SOURCE, T, 0, N + 2).removeO()
q3H = q_expr**3 * H_series
Lm1_series = sp.series(-SOURCE_ser / q3H, T, 0, N + 2).removeO()

Delta = sp.series(Lm1_series - L0_series, T, 0, N + 2).removeO()
I_Delta = sum(Delta.coeff(T, n) * T**(n+1) / (n+1) for n in range(N + 2))
I_Delta = sp.series(I_Delta, T, 0, N + 2).removeO()

print(f"\n=== Identity check at (u_1, u_2) = ({u1}, {u2}), so E_1={sv_int}, E_2={pv_int} ===")
print(f"Comparing R^(-1)[n] (Day 162 closed) vs Route_A[n] - I_Delta[n]")
all_ok = True
for n in range(1, N + 1):
    lhs = sp.Rational(Rmn_ser.coeff(T, n))
    rhs_A = sp.Rational(routeA_coefs[n].numerator, routeA_coefs[n].denominator) if routeA_coefs[n] != 0 else sp.S(0)
    rhs_I = sp.Rational(I_Delta.coeff(T, n))
    diff = lhs - (rhs_A - rhs_I)
    ok = (diff == 0)
    all_ok = all_ok and ok
    print(f"  n={n:2d}: R^(-1)={lhs}, RouteA-I_Delta={rhs_A - rhs_I}, diff={diff}  {'OK' if ok else 'FAIL'}")
print(f"\nOverall to N={N}: {'PASS' if all_ok else 'FAIL'}")
