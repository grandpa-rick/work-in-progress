"""Day 170 Step 8 — Treat Y and q as INDEPENDENT symbols with their own derivative rules.

Y = T*phi(Y), phi = 1 + sY + pY^2  → Y' = phi/q
q^2 = (1-sT)^2 - 4pT^2  → q' = -[s(1-sT) + 4pT]/q

Later reduce using q = 1 - sT - 2pTY (algebraic identity from linking Y and q).

Compute L_{-1} symbolically, then reduce.
"""
import sympy as sp
import time

T, s, p, Y, q = sp.symbols('T s p Y q')

phi_Y = 1 + s*Y + p*Y**2
Yp = phi_Y / q                                # Y'
qp = -(s*(1-s*T) + 4*p*T) / q                 # q'

def d_dT(expr):
    """Total T-derivative treating Y, q as independent functions of T."""
    return sp.diff(expr, T) + sp.diff(expr, Y) * Yp + sp.diff(expr, q) * qp

# Verify: d/dT of q^2 should equal 2*q*qp = -2[s(1-sT)+4pT] = derivative of (1-sT)^2 - 4pT^2
check = sp.expand(d_dT(q**2))
expected = sp.diff((1-s*T)**2 - 4*p*T**2, T)
diff = sp.simplify(check - expected)
print(f"d/dT(q^2) - expected = {sp.simplify(diff * q)/q}")  # divide out q

# H = pY/T
H = p*Y/T
Hp = d_dT(H)
Hpp = d_dT(Hp)

# K_{-1} = -pY/q^2
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

# Simplify
print("\nStep 1: sp.together / cancel SOURCE fraction...")
t0 = time.time()
Lm1 = sp.together(Lm1_raw)
print(f"  together: {time.time()-t0:.1f}s")

t0 = time.time()
num = sp.expand(sp.numer(Lm1))
den = sp.expand(sp.denom(Lm1))
print(f"  expand: {time.time()-t0:.1f}s")
print(f"  num degree Y = {sp.Poly(num, Y).degree() if Y in num.free_symbols else 0}, q = {sp.Poly(num, q).degree() if q in num.free_symbols else 0}, terms = {len(num.args) if num.is_Add else 1}")
print(f"  den = {sp.factor(den)}")

# The denominator should factor to something like T^k * q^m * Y^0 or so.
# Now reduce mod (Y-relation, q^2-relation).
# Y-rel: pTY^2 - (1-sT)Y + T = 0
# q^2-rel: q^2 - ((1-sT)^2 - 4pT^2) = 0
# We keep them as an ideal.

Y_relation = p*T*Y**2 - (1-s*T)*Y + T
q_relation = q**2 - ((1-s*T)**2 - 4*p*T**2)

print("\nStep 2: reduce num, den using Groebner basis of (Y-rel, q^2-rel) over Q(T,s,p)[Y,q]...")

# Use lex order with Y > q, so Y^2 gets reduced first, then q^2.
t0 = time.time()
# Compute Groebner basis
gb = sp.groebner([Y_relation, q_relation], Y, q, order='lex', domain=sp.QQ.frac_field(T, s, p))
print(f"  gb: {time.time()-t0:.1f}s")
print(f"  basis: {[str(g)[:80] for g in gb.polys]}")

t0 = time.time()
num_r = gb.reduce(sp.Poly(num, Y, q, domain=sp.QQ.frac_field(T, s, p)))[1]
den_r = gb.reduce(sp.Poly(den, Y, q, domain=sp.QQ.frac_field(T, s, p)))[1]
print(f"  reduce: {time.time()-t0:.1f}s")

print(f"\n  Reduced numerator:")
print(f"  {sp.expand(num_r.as_expr())}")
print(f"\n  Reduced denominator:")
print(f"  {sp.expand(den_r.as_expr())}")

# Denominator should now be poly in T, s, p (no Y, q).
Lm1_reduced = sp.expand(num_r.as_expr()) / sp.expand(den_r.as_expr())
Lm1_final = sp.cancel(Lm1_reduced)
print(f"\nL_{{-1}} = {Lm1_final}")

# Try to simplify: get as A + BY + Cq + DYq form
Lm1_num = sp.numer(Lm1_final)
Lm1_den = sp.denom(Lm1_final)
Lm1_num_pol = sp.Poly(sp.expand(Lm1_num), Y, q)
print(f"\nAs poly in Y, q:")
for (dY, dq), c in sorted(Lm1_num_pol.terms()):
    print(f"  Y^{dY} q^{dq}: {sp.factor(c)}")
print(f"\nDenom: {sp.factor(Lm1_den)}")

# Sanity check
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

N = 8
sv, pv = sp.Rational(2), sp.Rational(3)
Ys = Y_series(N, sv, pv)
Y_ser = sum(Ys[n] * T**n for n in range(N+1))
q_ser_val = 1 - sv*T - 2*pv*T*Y_ser  # = q as series
q_ser_val = sp.series(q_ser_val, T, 0, N).removeO()

Lm1_val = Lm1_final.subs({s: sv, p: pv, Y: Y_ser, q: q_ser_val})
Lm1_ser = sp.series(Lm1_val, T, 0, N).removeO()

L_actual = [
    sp.S(0), sp.S(0),
    -10*pv,
    -49*sv*pv,
    -145*sv**2*pv - 95*pv**2,
    -335*sv**3*pv - 658*sv*pv**2,
    -665*sv**4*pv - 2611*sv**2*pv**2 - 644*pv**3,
    -1190*sv**5*pv - 7784*sv**3*pv**2 - 5758*sv*pv**3,
]

print(f"\n=== Series verification ({N} terms) ===")
for n in range(min(N, len(L_actual))):
    computed = sp.expand(Lm1_ser.coeff(T, n))
    actual = sp.expand(L_actual[n])
    diff = sp.expand(computed - actual)
    print(f"  n={n}: comp={computed}, act={actual}, {'OK' if diff==0 else 'FAIL diff='+str(diff)}")
