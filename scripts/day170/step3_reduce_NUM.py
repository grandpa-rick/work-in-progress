"""Day 170 Step 3 — Reduce L_{-1} to a compact closed form.

Two key insights:
1. q = 1 - sT - 2pTY (algebraic identity, not just series — verified below).
2. Ring reduces to Q(s,p,T)[Y]/(pTY^2 + (sT-1)Y + T). Every element has unique normal form A(T,s,p) + B(T,s,p) Y.

Strategy for L_{-1}:
 - L_{-1} = -SOURCE / (q^3 H) where H = pY/T, K = K_{-1} = -pY/q^2.
 - Compute SOURCE fully symbolically using H' = ..., H'' = ..., K' = ... rules.
 - Divide by q^3 H, substitute q = 1 - sT - 2pTY.
 - Reduce Y^2, get normal form.

Verify against L_{-1} series from FP.
"""
import sympy as sp
import sys
sys.path.insert(0, '/home/agent/projects/scratch/day152')

T, s, p, Y, q = sp.symbols('T s p Y q')

# Y satisfies Y = T*phi(Y) = T + sTY + pTY^2, i.e. pTY^2 = (1-sT)Y - T.
# q^2 = (1-sT)^2 - 4pT^2, i.e. q = 1 - sT - 2pTY (positive branch).

# Derivative rules (chain rule with Y = Y(T)):
# Y' = phi(Y)/(1 - T*phi'(Y)) = phi(Y)/q  (since 1 - T*phi'(Y) = 1 - T(s + 2pY) = q)
# q' = -[s(1-sT) + 4pT]/q  (from 2q q' = derivative of (1-sT)^2 - 4pT^2)

phi_Y = 1 + s*Y + p*Y**2       # phi(Y) = 1 + sY + pY^2 = Y/T (from Y = T*phi(Y))
Yp = phi_Y / q                   # Y' (as function of T, Y, q, s, p)
qp = -(s*(1-s*T) + 4*p*T) / q    # q'

# H = pY/T
H = p*Y/T
# H' = d/dT (pY/T) = pY'/T - pY/T^2 = p*[Y'/T - Y/T^2]
Hp = p*Yp/T - p*Y/T**2
# H'' = derivative once more
def d_dT(expr):
    """Total derivative w.r.t. T, treating Y and q as functions of T with derivatives Yp, qp."""
    e = sp.diff(expr, T)                                        # explicit T-partial
    e += sp.diff(expr, Y) * Yp                                  # + Y' * d/dY
    e += sp.diff(expr, q) * qp                                  # + q' * d/dq
    return sp.together(e)

Hpp = d_dT(Hp)

# K_{-1} = -pY/q^2
K = -p*Y/q**2
Kp = d_dT(K)

# Auxiliary polynomials (Day 169):
# R_3 = -T^2 q^2  (verified: R_3 = -T^2 + 2sT^3 + (4p - s^2)T^4)
# R_2 = q^2 (1 - sT)
# R_1 = -p + 2ps T + (4p^2 - ps^2) T^2  = q^2 * (-p) + ... actually let's just define via given.
R3 = -T**2 * q**2
R2 = q**2 * (1 - s*T)
# R_1 not directly used in SOURCE below because Day 169 formula has R_1 only in top-layer check.

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

# L_{-1} = -SOURCE / (q^3 * H)
Lm1 = -SOURCE / (q**3 * H)
Lm1 = sp.together(Lm1)

print("Raw L_{-1} (before ring reduction):")
print(f"  num size (# args if Add) ≈ {len(sp.expand(sp.numer(Lm1)).args) if sp.expand(sp.numer(Lm1)).is_Add else 1}")
print(f"  denom = {sp.factor(sp.denom(Lm1))}")

# Now reduce in the ring Q(s,p,T)[Y,q]/(q^2 - ((1-sT)^2 - 4pT^2), Y - T - sTY - pTY^2).
# Actually the two-relation approach: reduce q^2 -> (1-sT)^2 - 4pT^2 (poly in T, s, p; no Y),
# then reduce Y^2 -> ((1-sT)Y - T)/(pT).

Lm1_expanded = sp.expand(Lm1)
Lm1_num = sp.expand(sp.numer(Lm1_expanded))
Lm1_den = sp.expand(sp.denom(Lm1_expanded))

print(f"\nNumerator num-terms (Add args): {len(Lm1_num.args) if Lm1_num.is_Add else 1}")
print(f"Denominator: {sp.factor(Lm1_den)}")

# Step A: reduce q^2 in Lm1_num using q^2 = (1-sT)^2 - 4pT^2
q2_val = (1-s*T)**2 - 4*p*T**2

def reduce_q2(expr):
    """Repeatedly substitute q^2 -> q2_val until no q^2 or higher remains."""
    e = sp.expand(expr)
    for _ in range(20):
        e_new = sp.expand(e.subs(q**2, q2_val))
        if e_new == e:
            break
        e = e_new
    return e

Lm1_num_r = reduce_q2(Lm1_num)
Lm1_den_r = reduce_q2(Lm1_den)

# Step B: reduce Y^2 using pTY^2 = (1-sT)Y - T, i.e. Y^2 = ((1-sT)Y - T)/(pT)
Y2_val = ((1-s*T)*Y - T)/(p*T)

def reduce_Y2(expr):
    e = sp.expand(expr)
    for _ in range(20):
        e_new = sp.expand(e.subs(Y**2, Y2_val))
        if e_new == e:
            break
        e = e_new
    return e

Lm1_num_rr = reduce_Y2(Lm1_num_r)
Lm1_den_rr = reduce_Y2(Lm1_den_r)

# Now numer and denom should have Y at most linear, q at most linear.
Lm1_num_pol = sp.Poly(Lm1_num_rr, Y, q)
Lm1_den_pol = sp.Poly(Lm1_den_rr, Y, q)

print("\n=== Reduced numerator (poly in Y, q) ===")
for (dY, dq), c in sorted(Lm1_num_pol.terms()):
    print(f"  Y^{dY} q^{dq}: {sp.factor(sp.expand(c))}")

print("\n=== Reduced denominator (poly in Y, q) ===")
for (dY, dq), c in sorted(Lm1_den_pol.terms()):
    print(f"  Y^{dY} q^{dq}: {sp.factor(sp.expand(c))}")

# Try to simplify: rationalize denominator if it has q terms.
# If denom = A(T,s,p) + B(T,s,p) q, multiply num & denom by (A - Bq), get A^2 - B^2 q^2 pure poly.
# Similarly for Y.
