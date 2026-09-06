"""Day 170 Step 13 — Correct SOURCE with the missing 18 T^3 H^2 K term.

Verify L_{-1} = -SOURCE/(q^3 H) matches L_actual now.
Then reduce to normal form in the ring.
"""
import sympy as sp
import time

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
    + 18 * T**3 * H**2 * K       # <-- MISSING TERM (P_3 G^3 e=1 layer 1)
)

Lm1_raw = -SOURCE / (q**3 * H)

# Series verification first
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

print("=== Series check with corrected SOURCE ===")
all_ok = True
for n in range(min(len(L_actual), N+1)):
    comp = sp.expand(Lm1_ser.coeff(T, n))
    act = sp.expand(L_actual[n])
    diff = sp.expand(comp - act)
    ok = (diff == 0)
    all_ok = all_ok and ok
    print(f"  n={n}: comp={comp}, act={act}, {'OK' if ok else 'FAIL diff='+str(diff)}")
print(f"\nOverall: {'PASS' if all_ok else 'FAIL'}")

if all_ok:
    print("\n=== Now reducing to normal form ===")
    print("Step 1: sp.cancel...")
    t0 = time.time()
    Lm1_rat = sp.cancel(Lm1_raw)
    print(f"  cancel: {time.time()-t0:.1f}s")
    num = sp.expand(sp.numer(Lm1_rat))
    den = sp.expand(sp.denom(Lm1_rat))
    print(f"  num deg Y = {sp.Poly(num, Y).degree() if Y in num.free_symbols else 0}, terms = {len(num.args) if num.is_Add else 1}")
    print(f"  den = {sp.factor(den)}")

    # Reduce using Groebner
    Y_rel = p*T*Y**2 - (1-s*T)*Y + T
    q_rel = q**2 - ((1-s*T)**2 - 4*p*T**2)
    print("\nStep 2: Groebner reduce mod (Y-rel, q^2-rel)...")
    t0 = time.time()
    gb = sp.groebner([Y_rel, q_rel], Y, q, order='lex', domain=sp.QQ.frac_field(T, s, p))
    print(f"  gb: {time.time()-t0:.1f}s")
    print(f"  gb basis:")
    for g in gb.polys:
        print(f"    {g.as_expr()}")

    num_red = gb.reduce(sp.Poly(num, Y, q, domain=sp.QQ.frac_field(T, s, p)))[1].as_expr()
    den_red = gb.reduce(sp.Poly(den, Y, q, domain=sp.QQ.frac_field(T, s, p)))[1].as_expr()
    num_red = sp.expand(num_red)
    den_red = sp.expand(den_red)
    print(f"\n  Reduced num: {num_red}")
    print(f"  Reduced den: {sp.factor(den_red)}")

    Lm1_reduced = num_red / den_red
    Lm1_final = sp.cancel(Lm1_reduced)
    num_pol = sp.Poly(sp.numer(Lm1_final), Y, q)
    print(f"\n=== L_{{-1}} normal form (poly in Y, q) ===")
    for (dY, dq), c in sorted(num_pol.terms()):
        cf = sp.factor(c)
        print(f"  Y^{dY} q^{dq}: {cf}")
    print(f"\nDenom: {sp.factor(sp.denom(Lm1_final))}")
