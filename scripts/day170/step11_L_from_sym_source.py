"""Day 170 Step 11 — Compute L_{-1} using sym SOURCE directly, no substitution shortcuts.

Test: substitute Y_ser and q_ser INTO -SOURCE/(q^3 H) SEPARATELY (numerator and denominator),
then divide as series.
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

# Substitute Y, q as series
sv, pv = sp.Rational(2), sp.Rational(3)
Nplus = 12
N = 10

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

Y_list = Y_series(Nplus)
Y_expr = sum(Y_list[n]*T**n for n in range(Nplus+1))
q_expr = 1 - sv*T - 2*pv*T*Y_expr
q_expr = sp.expand(q_expr)

# Substitute Y_expr and q_expr into SOURCE
SOURCE_sub = SOURCE.subs({s: sv, p: pv, Y: Y_expr, q: q_expr})
# Now SOURCE_sub has 1/T, 1/q, 1/Y from H, Hp, Hpp, K, Kp terms
# But it should be a power series in T (positive powers only, starting at T^2).
SOURCE_ser = sp.series(SOURCE_sub, T, 0, N+1).removeO()

# Print SOURCE coefficients
print("=== SOURCE (from sym after substitution) ===")
for n in range(N+1):
    print(f"  [T^{n}] = {sp.expand(SOURCE_ser.coeff(T, n))}")

# Compute L = -SOURCE / (q^3 * H) = -SOURCE * T / (p * Y * q^3)
# At series level: divide by pY (starts at T^1 with p=3) and q^3 (starts at 1).

# Method 1: substitute directly and use sp.series on the full ratio
denom_expr = pv * Y_expr * q_expr**3
denom_expr = sp.expand(denom_expr)

L_sub = -SOURCE_sub * T / denom_expr
L_ser = sp.series(L_sub, T, 0, N+1).removeO()

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

print("\n=== L_{-1} from sym ===")
for n in range(min(len(L_actual), N+1)):
    comp = sp.expand(L_ser.coeff(T, n))
    act = sp.expand(L_actual[n])
    diff = sp.expand(comp - act)
    print(f"  n={n}: comp={comp}, act={act}, {'OK' if diff==0 else 'FAIL diff='+str(diff)}")
