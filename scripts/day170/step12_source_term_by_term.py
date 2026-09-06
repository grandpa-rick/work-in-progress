"""Day 170 Step 12 — Compare SOURCE term by term. Find which subterm is wrong."""
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

TERMS = {
    'R3_Hpp': R3 * Hpp,
    'coefHp_Hp': coef_Hp * Hp,
    'coefH_H': coef_H * H,
    '3R3_HKp+KHp': 3 * R3 * (H * Kp + K * Hp),
    'coefH2_H^2': coef_H2 * H**2,
    '18T3_HHp': 18 * T**3 * H * Hp,
    'T4_H^3': T**4 * H**3,
    '3R3_HK^2': 3 * R3 * H * K**2,
    'R2_Kp': R2 * Kp,
    '2coefHp_HK': coef_Hp * 2 * H * K,
    'coefK_K': coef_K * K,
    'R2_K^2': R2 * K**2,
}

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
q_expr = sp.expand(1 - sv*T - 2*pv*T*Y_expr)

def to_ser_at_n(expr, n):
    e = expr.subs({s: sv, p: pv, Y: Y_expr, q: q_expr})
    return sp.series(e, T, 0, n+1).removeO().coeff(T, n)

# Also compute the same terms via SERIES arithmetic (step 31c style)
def series_mul(f, g, N):
    return [sum(f[k]*g[n-k] for k in range(n+1)) for n in range(N+1)]
def series_inv(f, N):
    inv = [sp.S(0)]*(N+1); inv[0] = 1/f[0]
    for n in range(1, N+1):
        s0 = sp.S(0)
        for k in range(1, n+1):
            s0 -= f[k]*inv[n-k]
        inv[n] = s0/f[0]
    return inv
def series_deriv(f, N):
    return [(k+1)*f[k+1] if k+1 < len(f) else sp.S(0) for k in range(N+1)]

def q_series_(N):
    q2 = [sp.S(0)] * (N+1); q2[0] = 1; q2[1] = -2*sv; q2[2] = sv**2 - 4*pv
    qq = [sp.S(0)] * (N+1); qq[0] = 1
    for n in range(1, N+1):
        s0 = q2[n] if n <= 2 else sp.S(0)
        for k in range(1, n):
            s0 -= qq[k]*qq[n-k]
        qq[n] = s0/2
    return qq

qs = q_series_(Nplus)
Ys = Y_list
q2s = series_mul(qs, qs, Nplus)
q2_inv = series_inv(q2s, Nplus)
Hs = [pv*Ys[k+1] if k+1 <= Nplus else sp.S(0) for k in range(Nplus+1)]
Ks = series_mul([-pv*Ys[k] for k in range(Nplus+1)], q2_inv, Nplus)
Hps = series_deriv(Hs, Nplus)
Hpps = series_deriv(Hps, Nplus)
Kps = series_deriv(Ks, Nplus)

R3s_ser = [sp.S(0), sp.S(0), -1, 2*sv, 4*pv-sv**2] + [sp.S(0)]*(Nplus-3)
R2s_ser = [1, -3*sv, 3*sv**2-4*pv, 4*pv*sv-sv**3] + [sp.S(0)]*(Nplus-2)
coefHp_ser = [0, -11, 14*sv, 12*pv-3*sv**2] + [sp.S(0)]*(Nplus-2)
coefH_ser = [1, 12*sv, 5*pv-sv**2] + [sp.S(0)]*(Nplus-1)
coefH2_ser = [0, 0, 23, sv] + [sp.S(0)]*(Nplus-2)
coefK_ser = [-sv, 2*sv**2+10*pv, 4*pv*sv-sv**3] + [sp.S(0)]*(Nplus-1)

def sh(f, k):
    return [sp.S(0)]*k + list(f)

R3_Hpp_ser = series_mul(R3s_ser, Hpps, Nplus)
coefHp_Hp_ser = series_mul(coefHp_ser, Hps, Nplus)
coefH_H_ser = series_mul(coefH_ser, Hs, Nplus)
HKp_ser = series_mul(Hs, Kps, Nplus)
KHp_ser = series_mul(Ks, Hps, Nplus)
R3_HKp_KHp_ser = series_mul(R3s_ser, [a+b for a, b in zip(HKp_ser, KHp_ser)], Nplus)
H2 = series_mul(Hs, Hs, Nplus)
H3 = series_mul(H2, Hs, Nplus)
HK = series_mul(Hs, Ks, Nplus)
K2 = series_mul(Ks, Ks, Nplus)
HK2 = series_mul(Hs, K2, Nplus)
HHp = series_mul(Hs, Hps, Nplus)
coefH2_H2_ser = series_mul(coefH2_ser, H2, Nplus)
T3_HHp_ser = sh(HHp, 3)[:Nplus+1]
T4_H3_ser = sh(H3, 4)[:Nplus+1]
R3_HK2_ser = series_mul(R3s_ser, HK2, Nplus)
R2_Kp_ser = series_mul(R2s_ser, Kps, Nplus)
coefHp_HK_ser = series_mul(coefHp_ser, HK, Nplus)
coefK_K_ser = series_mul(coefK_ser, Ks, Nplus)
R2_K2_ser = series_mul(R2s_ser, K2, Nplus)

TERMS_ser = {
    'R3_Hpp': R3_Hpp_ser,
    'coefHp_Hp': coefHp_Hp_ser,
    'coefH_H': coefH_H_ser,
    '3R3_HKp+KHp': [3*x for x in R3_HKp_KHp_ser],
    'coefH2_H^2': coefH2_H2_ser,
    '18T3_HHp': [18*x for x in T3_HHp_ser],
    'T4_H^3': T4_H3_ser,
    '3R3_HK^2': [3*x for x in R3_HK2_ser],
    'R2_Kp': R2_Kp_ser,
    '2coefHp_HK': [2*x for x in coefHp_HK_ser],
    'coefK_K': coefK_K_ser,
    'R2_K^2': R2_K2_ser,
}

for name, expr in TERMS.items():
    ser = TERMS_ser[name]
    sym_val = sp.expand(to_ser_at_n(expr, 4))
    ser_val = sp.expand(ser[4])
    diff = sp.expand(sym_val - ser_val)
    print(f"{name}: sym[T^4]={sym_val}, ser[T^4]={ser_val}, diff={diff}")
