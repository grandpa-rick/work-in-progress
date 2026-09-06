"""Day 170 Step 4 — Cleaner reduction of L_{-1}.

Approach: Represent everything as A + B*Y + C*q + D*Y*q with A, B, C, D rational in T, s, p.
Then use rewriting rules Y^2 = ((1-sT)Y - T)/(pT), q^2 = (1-sT)^2 - 4pT^2, and
Y*q reducible via q = 1-sT-2pTY.

But simpler: work over Q(T,s,p)[Y] since q = 1-sT-2pTY is an algebraic identity.
Then every rational in Y (with Y^2 reducible) has normal form A + B*Y.
"""
import sympy as sp
import sys
sys.path.insert(0, '/home/agent/projects/scratch/day152')

T, s, p, Y = sp.symbols('T s p Y')

# q as a function of Y (algebraic identity)
q_of_Y = 1 - s*T - 2*p*T*Y

# Y satisfies pTY^2 = (1-sT)Y - T, so Y^2 = ((1-sT)Y - T)/(pT)
Y2_val = ((1-s*T)*Y - T) / (p*T)

def reduce_Y_poly(expr):
    """Reduce Y^2 -> ((1-sT)Y - T)/(pT). Result is degree <= 1 in Y."""
    e = sp.expand(expr)
    for _ in range(60):
        e_new = sp.expand(e.subs(Y**2, Y2_val))
        if e_new == e:
            break
        e = e_new
    return e

def to_AB_form(expr):
    """Reduce expression to A + B*Y form (over Q(T,s,p)). Handles rationals in Y via cross-multiply."""
    # Cross-multiply if needed
    expr_together = sp.together(sp.simplify(expr))
    num = sp.expand(sp.numer(expr_together))
    den = sp.expand(sp.denom(expr_together))
    # Reduce Y^2 in num and den
    num_r = reduce_Y_poly(num)
    den_r = reduce_Y_poly(den)
    # If denominator has Y term, rationalize
    den_poly = sp.Poly(den_r, Y)
    if den_poly.degree() >= 1:
        # den_r = a + b Y (a, b in Q(T, s, p))
        a_coef = den_poly.nth(0)
        b_coef = den_poly.nth(1)
        # Multiply top and bottom by conjugate: use Y satisfies X^2 = ((1-sT)Y - T)/(pT).
        # Actually, we want to eliminate Y from denom. Since Y^2 - ((1-sT)/(pT)) Y + 1/p = 0.
        # So (a + bY)(a + b·Y_conj) = a^2 + a·b·(Y+Y_conj) + b^2·Y·Y_conj
        # If Y, Y_conj are roots: sum = (1-sT)/(pT), product = 1/p.
        # So (a + bY)(a + bY_conj) = a^2 + ab(1-sT)/(pT) + b^2/p.
        # Which is polynomial in T, s, p (no Y).
        conj_factor = a_coef + b_coef * ((1-s*T)/(p*T) - Y)  # a + b*(sum - Y) = a + b*Y_conj
        num_new = sp.expand(num_r * conj_factor)
        den_new = sp.expand(den_r * conj_factor)
        num_new = reduce_Y_poly(num_new)
        den_new = reduce_Y_poly(den_new)
        # den_new should now be Y-free
        den_poly2 = sp.Poly(den_new, Y)
        if den_poly2.degree() >= 1:
            # Iterate
            num_new = reduce_Y_poly(num_new)
            den_new = reduce_Y_poly(den_new)
            den_poly2 = sp.Poly(den_new, Y)
            assert den_poly2.degree() == 0, f"still Y in denom: {den_new}"
        # Now num_new / den_new with den_new poly in T, s, p only.
        num_pol = sp.Poly(num_new, Y)
        A = sp.together(num_pol.nth(0) / den_new)
        B = sp.together(num_pol.nth(1) / den_new)
        return A, B
    else:
        # den is poly in T, s, p only
        num_pol = sp.Poly(num_r, Y)
        A = sp.together(num_pol.nth(0) / den_r)
        B = sp.together(num_pol.nth(1) / den_r)
        return A, B

# --- Test the reduction on Y itself ---
A, B = to_AB_form(Y)
print(f"Y = {A} + ({B})*Y")

A, B = to_AB_form(q_of_Y)  # q = 1 - sT - 2pTY
print(f"q = {A} + ({B})*Y")

A, B = to_AB_form(1/q_of_Y)
print(f"1/q = {A} + ({B})*Y")
print(f"   check by multiplying: {sp.expand((A + B*Y) * q_of_Y - 1)}")
print(f"   check reduced: {reduce_Y_poly(sp.expand((A + B*Y) * q_of_Y - 1))}")
