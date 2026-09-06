"""Day 170 Step 1 — Verify the full identity to high N.

Identity to check:
  R^{(-1)}(T) = Route_A(T) - I_Delta(T)

where:
  R^{(-1)}      = Day 162 closed form = T[E_2 Y^2((q+1)^2 - E_1 T) + (q + R_1R_2)/2]/q^3
  Route_A       = (1/2) d^2_{u_3} Xi|_{u_3=0} = Day 167 chain-rule closed form (uses xi_0, xi_1, xi_2)
  I_Delta       = int_0^T (L_{-1}(T') - L_0(T')) dT'   with L_0 (Day 168) and L_{-1} (Day 169)

Everything computed as sympy series in T.
"""
import sympy as sp

N = 12  # series precision — we'll check up through this order

E1s, E2s, T = sp.symbols('E1 E2 T')

# ----- Basic series -----
Y_expr = sp.S(0)
for _ in range(N+2):
    Y_expr = sp.expand(T*(1 + E1s*Y_expr + E2s*Y_expr**2))
    Y_expr = sp.series(Y_expr, T, 0, N+3).removeO()

q_sp = sp.expand(1 - E1s*T - 2*T*E2s*Y_expr)
phi_sp = 1 + E1s*Y_expr + E2s*Y_expr**2
log_q = sp.series(sp.log(q_sp), T, 0, N+2).removeO()
log_YoverTq = sp.series(sp.log(phi_sp / q_sp), T, 0, N+2).removeO()

# ----- xi_0 = int_0^T E_2 Y/T dT'  (i.e. sum_n E_2 Y_n/n T^n) -----
xi_0_series = sum(E2s * Y_expr.coeff(T, n) * T**n / n for n in range(1, N+2))
d_E1_xi0 = sp.diff(xi_0_series, E1s)
d_E2_xi0 = sp.diff(xi_0_series, E2s)
d_E1E1_xi0 = sp.diff(xi_0_series, E1s, 2)
d_E1E2_xi0 = sp.diff(xi_0_series, E1s, E2s)
d_E2E2_xi0 = sp.diff(xi_0_series, E2s, 2)

# ----- xi_1|_{E_3=0} closed form: E_2 xi_1 = -log q - d_E1 xi_0 - E_1 d_E2 xi_0 -----
E2xi1 = sp.expand(-log_q - d_E1_xi0 - E1s * d_E2_xi0)
xi_1_closed = sp.series(sp.expand(E2xi1 / E2s), T, 0, N+2).removeO()
d_E1_xi1 = sp.diff(xi_1_closed, E1s)
d_E2_xi1 = sp.diff(xi_1_closed, E2s)

# ----- [E_3^1] log W: E_2 [E_3^1]logW = T(q+R_1R_2)/q^3 - d_E1 log(Y/(Tq)) - E_1 d_E2 log(Y/(Tq)) -----
R1R2 = 1 - T**2*(E1s**2 - 4*E2s)
q3_inv = sp.series(1/q_sp**3, T, 0, N+2).removeO()
Tqinv_series = sp.series(T*(q_sp + R1R2)*q3_inv, T, 0, N+2).removeO()
d_E1_logW = sp.diff(log_YoverTq, E1s)
d_E2_logW = sp.diff(log_YoverTq, E2s)
E2_E31_logW = sp.expand(Tqinv_series - d_E1_logW - E1s * d_E2_logW)
E31_logW_closed = sp.series(sp.expand(E2_E31_logW / E2s), T, 0, N+2).removeO()

# ----- xi_2|_{E_3=0} closed form: 2 E_2 xi_2 = [E_3^1]logW - 3 d_E1 xi_1 - 2 E_1 d_E2 xi_1 -----
E2xi2 = sp.expand((E31_logW_closed - 3*d_E1_xi1 - 2*E1s*d_E2_xi1) / 2)
xi_2_closed = sp.series(sp.expand(E2xi2 / E2s), T, 0, N+2).removeO()

# ----- Route A = (1/2) d^2_{u_3} Xi|_{u_3=0} -----
# = (1/2) [xi_0_{E1E1} + 2 E1 xi_0_{E1E2} + E1^2 xi_0_{E2E2}
#          + 2 E2 (xi_1_{E1} + E1 xi_1_{E2})
#          + 2 E2^2 xi_2 ]
routeA = sp.series(
    (d_E1E1_xi0 + 2*E1s*d_E1E2_xi0 + E1s**2*d_E2E2_xi0
     + 2*E2s*(d_E1_xi1 + E1s*d_E2_xi1)
     + 2*E2s**2*xi_2_closed) / 2,
    T, 0, N+2).removeO()

# ----- R^{(-1)} = T[E_2 Y^2((q+1)^2 - E_1 T) + (q + R_1R_2)/2]/q^3  (Day 162) -----
inner = E2s*Y_expr**2 * ((q_sp+1)**2 - E1s*T) + sp.Rational(1,2)*(q_sp + R1R2)
Rmn_series = sp.series(sp.expand(T * inner * q3_inv), T, 0, N+2).removeO()

# ----- L_0 = (1 + 3 T K + T^2 K^2 + T theta K)/q  where K = [pY(2q+1) + s q]/q^2 (Day 168) -----
K0_expr = (E2s*Y_expr*(2*q_sp+1) + E1s*q_sp) / q_sp**2  # Day 158
K0_series = sp.series(K0_expr, T, 0, N+2).removeO()
theta_K0 = sum(n * K0_series.coeff(T, n) * T**n for n in range(N+2))
L0_num = 1 + 3*T*K0_series + T**2 * K0_series**2 + T * theta_K0
L0_series = sp.series(L0_num / q_sp, T, 0, N+2).removeO()

# ----- L_{-1} needs care. From Day 169:
# H = pY/T (same as day158 H_0), K_{-1} = -pY/q^2.
# SOURCE = R3 H'' + [-11T+14sT^2+(12p-3s^2)T^3] H' + [1+12sT+(5p-s^2)T^2] H
#        + 3 R3 (H K' + K H') + [23T^2+sT^3] H^2 + 18 T^3 H H' + T^4 H^3 + 3 R3 H K^2
#        + R2 K' + [-11T+14sT^2+(12p-3s^2)T^3]*2*H*K + [-s+(2s^2+10p)T+(4ps-s^3)T^2] K
#        + R2 K^2
# where R2 = q^2(1-sT), R3 = -T^2 q^2.
# L_{-1} = -SOURCE / (q^3 H)
H_expr = E2s * Y_expr / T  # = pY/T
K_expr = -E2s * Y_expr / q_sp**2  # = -pY/q^2

# For derivatives, compute as series and differentiate directly
H_series = sp.series(H_expr, T, 0, N+2).removeO()
Hp_series = sp.diff(H_series, T)
Hpp_series = sp.diff(Hp_series, T)

K_series = sp.series(K_expr, T, 0, N+2).removeO()
Kp_series = sp.diff(K_series, T)

R3 = -T**2 * q_sp**2
R2 = q_sp**2 * (1 - E1s*T)

coef_Hp = -11*T + 14*E1s*T**2 + (12*E2s - 3*E1s**2) * T**3
coef_H  = 1 + 12*E1s*T + (5*E2s - E1s**2) * T**2
coef_H2 = 23*T**2 + E1s*T**3
coef_K  = -E1s + (2*E1s**2 + 10*E2s) * T + (4*E2s*E1s - E1s**3) * T**2

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
SOURCE_ser = sp.series(SOURCE, T, 0, N+2).removeO()
q3_H = q_sp**3 * H_series
Lm1_series = sp.series(-SOURCE_ser / q3_H, T, 0, N+2).removeO()

# ----- I_Delta = int_0^T (L_{-1} - L_0) dT' -----
Delta = Lm1_series - L0_series
Delta = sp.series(Delta, T, 0, N+2).removeO()
I_Delta = sp.S(0)
for n in range(N+2):
    c = Delta.coeff(T, n)
    if c:
        I_Delta += c * T**(n+1) / (n+1)
I_Delta = sp.series(I_Delta, T, 0, N+2).removeO()

# ----- Identity check: R^{(-1)} == Route_A - I_Delta ? -----
print(f"=== Identity check to N={N} ===")
print(f"R^(-1)  [Day 162 closed]  vs  Route_A [Day 167] - I_Delta [Day 169]")
print()
all_ok = True
for n in range(1, N+1):
    lhs = sp.expand(Rmn_series.coeff(T, n))
    rhs = sp.expand((routeA - I_Delta).coeff(T, n))
    diff = sp.expand(lhs - rhs)
    ok = (diff == 0)
    all_ok = all_ok and ok
    print(f"  n={n}: {'OK' if ok else 'FAIL — diff=' + str(diff)[:200]}")
print()
print(f"Overall to N={N}: {'PASS' if all_ok else 'FAIL'}")
