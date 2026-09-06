"""Day 170 Step 10 — Directly evaluate symbolic L_{-1} at series level.

If symbolic L_{-1} matches L_actual at series level, then the derivation is correct.
Then I know the bug is in the ring-reduction step.
"""
import sympy as sp
import sys
sys.path.insert(0, '/home/agent/projects/scratch/day152')

T, s, p, Y, q = sp.symbols('T s p Y q')

phi_Y = 1 + s*Y + p*Y**2
Yp = phi_Y / q
qp = -(s*(1-s*T) + 4*p*T) / q

def d_dT(expr):
    return sp.diff(expr, T) + sp.diff(expr, Y) * Yp + sp.diff(expr, q) * qp

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

Lm1_raw = -SOURCE / (q**3 * H)

# Evaluate at s=2, p=3 with Y and q as SERIES in T
sv, pv = sp.Rational(2), sp.Rational(3)
Nplus = 10
N = Nplus - 2

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

Y_ser_list = Y_series(Nplus)
Y_expr = sum(Y_ser_list[n]*T**n for n in range(Nplus+1))
q_expr = 1 - sv*T - 2*pv*T*Y_expr  # exact algebraic identity
q_expr = sp.expand(q_expr)

# Substitute Y, q into Lm1_raw as series
Lm1_sub = Lm1_raw.subs({s: sv, p: pv, Y: Y_expr, q: q_expr})
Lm1_ser = sp.series(Lm1_sub, T, 0, N+1).removeO()

L_actual = [
    sp.S(0), sp.S(0),
    -10*pv,
    -49*sv*pv,
    -145*sv**2*pv - 95*pv**2,
    -335*sv**3*pv - 658*sv*pv**2,
    -665*sv**4*pv - 2611*sv**2*pv**2 - 644*pv**3,
    -1190*sv**5*pv - 7784*sv**3*pv**2 - 5758*sv*pv**3,
    -1974*sv**6*pv - 19362*sv**4*pv**2 - 28638*sv**2*pv**3 - 3777*pv**4,
]

print("=== Symbolic L_{-1} vs L_actual ===")
for n in range(min(len(L_actual), N+1)):
    comp = sp.expand(Lm1_ser.coeff(T, n))
    act = sp.expand(L_actual[n])
    diff = sp.expand(comp - act)
    print(f"  n={n}: comp={comp}, act={act}, {'OK' if diff==0 else 'FAIL diff='+str(diff)}")
