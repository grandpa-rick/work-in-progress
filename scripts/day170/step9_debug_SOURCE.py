"""Day 170 Step 9 — Debug SOURCE mismatch. Compare symbolic-SOURCE to series-SOURCE."""
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

# Numerical values
sv, pv = sp.Rational(2), sp.Rational(3)
N = 8
Nplus = N + 2

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

def q_series(N):
    q2 = [sp.S(0)] * (N+1)
    q2[0] = 1; q2[1] = -2*sv; q2[2] = sv**2 - 4*pv
    q = [sp.S(0)] * (N+1); q[0] = 1
    for n in range(1, N+1):
        s0 = q2[n] if n <= 2 else sp.S(0)
        for k in range(1, n):
            s0 -= q[k]*q[n-k]
        q[n] = s0/2
    return q

Y_ser = Y_series(Nplus)
q_ser = q_series(Nplus)
Y_expr = sum(Y_ser[n]*T**n for n in range(N+1))
q_expr = sum(q_ser[n]*T**n for n in range(N+1))

# Verify q_expr matches 1 - sT - 2pTY:
q_alt = 1 - sv*T - 2*pv*T*Y_expr
for k in range(N+1):
    a = sp.expand(q_expr).coeff(T, k)
    b = sp.expand(q_alt).coeff(T, k)
    assert a == b, f"k={k}: q_ser={a} vs q_alt={b}"

# Substitute Y, q into H, Hp, Hpp, K, Kp
def to_ser(expr):
    e = expr.subs({s: sv, p: pv, Y: Y_expr, q: q_expr})
    return sp.series(e, T, 0, N+1).removeO()

H_s = to_ser(H)
Hp_s = to_ser(Hp)
Hpp_s = to_ser(Hpp)
K_s = to_ser(K)
Kp_s = to_ser(Kp)

# Now do the same via series arithmetic (like step 31c)
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

q2s = series_mul(q_ser, q_ser, N)
q2_inv = series_inv(q2s, N)

# H = pY/T
H_ser2 = [pv * Y_ser[n+1] for n in range(N+1)]  # (pY/T)[n] = p*Y[n+1]
Hp_ser2 = series_deriv(H_ser2, N)
Hpp_ser2 = series_deriv(Hp_ser2, N)

# K = -pY/q^2. Y as series starts at T^1 with 1.
K_ser2 = series_mul([-pv * Y_ser[n] for n in range(N+1)], q2_inv, N)  # -pY * 1/q^2
Kp_ser2 = series_deriv(K_ser2, N)

# Compare
def cmp(name, e1, l2):
    print(f"\n{name}:")
    for n in range(min(N+1, 6)):
        a = sp.expand(e1.coeff(T, n))
        b = sp.expand(l2[n])
        diff = sp.expand(a - b)
        print(f"  [T^{n}] sym={a}, ser={b}, {'OK' if diff==0 else 'DIFF='+str(diff)}")

cmp('H', H_s, H_ser2)
cmp('Hp', Hp_s, Hp_ser2)
cmp('Hpp', Hpp_s, Hpp_ser2)
cmp('K', K_s, K_ser2)
cmp('Kp', Kp_s, Kp_ser2)

# Compare SOURCE full assembly
print("\n\n=== SOURCE comparison ===")
SOURCE_sym = to_ser(SOURCE)

# Now do SOURCE via series arithmetic:
def series_scal(f, c):
    return [c*x for x in f]
def series_add(*fs):
    N0 = max(len(f) for f in fs)
    return [sum((f[k] if k < len(f) else sp.S(0)) for f in fs) for k in range(N0)]
def series_shift(f, k):
    return [sp.S(0)]*k + list(f)
def cap(f, N):
    return f[:N+1] + [sp.S(0)]*(N+1-len(f))

R3s = cap([sp.S(0), sp.S(0), -1, 2*sv, 4*pv - sv**2] + [sp.S(0)]*(N-3), N)
R2s = cap([1, -3*sv, 3*sv**2 - 4*pv, 4*pv*sv - sv**3] + [sp.S(0)]*(N-2), N)

coef_Hp_pol = -11*T + 14*sv*T**2 + (-3*sv**2 + 12*pv) * T**3
coef_H_pol = 1 + 12*sv*T + (-sv**2 + 5*pv) * T**2
coef_H2_pol = 23*T**2 + sv*T**3
coef_K_pol = -sv + (2*sv**2 + 10*pv) * T + (-sv**3 + 4*pv*sv) * T**2

def poly_to_ser(pl, N):
    if pl == 0:
        return [sp.S(0)]*(N+1)
    coefs = sp.Poly(pl, T).all_coeffs()[::-1]
    return coefs + [sp.S(0)]*(N+1 - len(coefs))

Hp_ser2 = Hp_ser2[:N+1]
Hpp_ser2 = Hpp_ser2[:N+1]
Kp_ser2 = Kp_ser2[:N+1]
H_ser2 = H_ser2[:N+1]
K_ser2 = K_ser2[:N+1]

HK = series_mul(H_ser2, K_ser2, N)
HHp = series_mul(H_ser2, Hp_ser2, N)
KHp = series_mul(K_ser2, Hp_ser2, N)
HKp = series_mul(H_ser2, Kp_ser2, N)
K2 = series_mul(K_ser2, K_ser2, N)
H2 = series_mul(H_ser2, H_ser2, N)
H3 = series_mul(H2, H_ser2, N)
HK2 = series_mul(H_ser2, K2, N)

R3_Hpp = series_mul(R3s, Hpp_ser2, N)
coef_Hp_ser = poly_to_ser(coef_Hp_pol, N)
coefHp_Hp = series_mul(coef_Hp_ser, Hp_ser2, N)
coef_H_ser = poly_to_ser(coef_H_pol, N)
coefH_H = series_mul(coef_H_ser, H_ser2, N)
R3_HKp_KHp = series_mul(R3s, series_add(HKp, KHp), N)
coef_H2_ser = poly_to_ser(coef_H2_pol, N)
coefH2_H2 = series_mul(coef_H2_ser, H2, N)
T3_HHp = cap(series_shift(HHp, 3), N)
T4_H3 = cap(series_shift(H3, 4), N)
R3_HK2 = series_mul(R3s, HK2, N)
R2_Kp = series_mul(R2s, Kp_ser2, N)
R2_K2 = series_mul(R2s, K2, N)
coefHp_HK = series_mul(coef_Hp_ser, HK, N)
coef_K_ser = poly_to_ser(coef_K_pol, N)
coefK_K = series_mul(coef_K_ser, K_ser2, N)

SOURCE_ser = series_add(
    R3_Hpp,
    coefHp_Hp,
    coefH_H,
    series_scal(R3_HKp_KHp, 3),
    coefH2_H2,
    series_scal(T3_HHp, 18),
    T4_H3,
    series_scal(R3_HK2, 3),
    R2_Kp,
    series_scal(coefHp_HK, 2),
    coefK_K,
    R2_K2,
)

for n in range(N+1):
    sym_val = sp.expand(SOURCE_sym.coeff(T, n))
    ser_val = sp.expand(SOURCE_ser[n])
    diff = sp.expand(sym_val - ser_val)
    print(f"  n={n}: sym={sym_val}, ser={ser_val}, {'OK' if diff==0 else 'DIFF='+str(diff)}")
