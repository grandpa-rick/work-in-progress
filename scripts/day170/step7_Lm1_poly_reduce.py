"""Day 170 Step 7 — L_{-1} normal form using sp.Poly + polynomial division.

Better approach: extract numerator and denominator as sp.Poly in Y (over Q(T,s,p)),
then reduce mod I = (pTY^2 - (1-sT)Y + T) using polynomial rem.
Then rationalize denominator with conjugate.
"""
import sympy as sp
import time

T, s, p, Y = sp.symbols('T s p Y')

q = 1 - s*T - 2*p*T*Y

phi_Y = 1 + s*Y + p*Y**2
Yp = phi_Y / q

def d_dT(expr):
    return sp.diff(expr, T) + sp.diff(expr, Y) * Yp

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

Lm1 = -SOURCE / (q**3 * H)

# Use sp.cancel to get single rational form (may take a while)
print("Step 1: sp.cancel...")
t0 = time.time()
Lm1_rat = sp.cancel(Lm1)
print(f"  done in {time.time()-t0:.1f}s")

num = sp.numer(Lm1_rat)
den = sp.denom(Lm1_rat)

# Convert to sp.Poly in Y, treating T, s, p as symbols
# But we need FIELD Q(T, s, p), so we use extension or fractions.
# Let me use sp.Poly with modulus over Q(T,s,p).

# Actually simpler: reduce Y-powers using rem.
Y_relation = p*T*Y**2 - (1 - s*T)*Y + T  # = 0

def reduce_mod_Y_rel(expr):
    """Reduce mod pTY^2 - (1-sT)Y + T = 0. Result is deg <= 1 in Y."""
    pol = sp.Poly(expr, Y)
    rel = sp.Poly(Y_relation, Y)
    _, r = sp.div(pol, rel, Y)
    return sp.expand(r.as_expr())

print("\nStep 2: reduce num, den mod Y-rel...")
t0 = time.time()
num_red = reduce_mod_Y_rel(num)
den_red = reduce_mod_Y_rel(den)
print(f"  done in {time.time()-t0:.1f}s")
print(f"  num after reduction: {sp.Poly(num_red, Y).degree()}, terms={len(sp.expand(num_red).args) if sp.expand(num_red).is_Add else 1}")
print(f"  den after reduction: {sp.Poly(den_red, Y).degree()}, terms={len(sp.expand(den_red).args) if sp.expand(den_red).is_Add else 1}")

# Now den has degree <= 1 in Y. Rationalize.
print("\nStep 3: rationalize denominator...")
den_pol = sp.Poly(den_red, Y)
a = den_pol.nth(0)
b = den_pol.nth(1) if den_pol.degree() >= 1 else sp.S(0)
# conjugate: replace Y with Y_conj = (1-sT)/(pT) - Y (sum minus Y)
# (a + bY)(a + b·Y_conj): use Y*Y_conj = 1/p, Y + Y_conj = (1-sT)/(pT).
# But rather compute directly:
conj = a + b * ((1-s*T)/(p*T) - Y)
num_r2 = num_red * conj
den_r2 = den_red * conj

t0 = time.time()
num_r2_red = reduce_mod_Y_rel(num_r2)
den_r2_red = reduce_mod_Y_rel(den_r2)
print(f"  reduce again: {time.time()-t0:.1f}s")

den_pol2 = sp.Poly(den_r2_red, Y)
print(f"  den degree in Y: {den_pol2.degree()}")
assert den_pol2.degree() <= 0, f"still has Y: {den_r2_red}"

D_final = sp.expand(den_r2_red)
num_pol = sp.Poly(num_r2_red, Y)
A_num = sp.expand(num_pol.nth(0))
B_num = sp.expand(num_pol.nth(1) if num_pol.degree() >= 1 else 0)

print("\nStep 4: cancel common factors in A/D and B/D...")
t0 = time.time()
A = sp.cancel(A_num / D_final)
B = sp.cancel(B_num / D_final)
print(f"  done in {time.time()-t0:.1f}s")

print("\n=== L_{-1} = A + B*Y ===")
print(f"A = {A}")
print(f"\nB = {B}")

# Print factored forms
print("\n=== Factored ===")
print(f"A: num = {sp.factor(sp.numer(A))}")
print(f"A: den = {sp.factor(sp.denom(A))}")
print(f"B: num = {sp.factor(sp.numer(B))}")
print(f"B: den = {sp.factor(sp.denom(B))}")

# Save
with open('/home/agent/projects/scratch/day170/Lm1_normal_form.py', 'w') as f:
    f.write("# L_{-1} normal form in Q(T, s, p)[Y]/(pTY^2 - (1-sT)Y + T)\n")
    f.write(f"A = {A}\n")
    f.write(f"B = {B}\n")
    f.write("# L_{-1} = A + B*Y\n")

# Sanity check: at series level, does A + B*Y match L_{-1} closed form?
print("\n=== Sanity check at series ===")
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

A_val = A.subs({s: sv, p: pv})
B_val = B.subs({s: sv, p: pv})
A_ser = sp.series(A_val, T, 0, N).removeO()
B_ser = sp.series(B_val, T, 0, N).removeO()
Lm1_expr = A_ser + B_ser * Y_ser
Lm1_expr = sp.series(Lm1_expr, T, 0, N).removeO()

L_actual = [
    sp.S(0), sp.S(0),
    -10*pv,
    -49*sv*pv,
    -145*sv**2*pv - 95*pv**2,
    -335*sv**3*pv - 658*sv*pv**2,
    -665*sv**4*pv - 2611*sv**2*pv**2 - 644*pv**3,
    -1190*sv**5*pv - 7784*sv**3*pv**2 - 5758*sv*pv**3,
]

for n in range(min(N, len(L_actual))):
    computed = sp.expand(Lm1_expr.coeff(T, n))
    actual = sp.expand(L_actual[n])
    diff = sp.expand(computed - actual)
    print(f"  n={n}: computed={computed}, actual={actual}, {'OK' if diff==0 else 'FAIL diff='+str(diff)}")
