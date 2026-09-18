"""
Day 201: Take Rick's Lemma 1 closed form for p_2(Y).e_r(X), substitute
q = t^alpha, expand in eps = t - 1 to O(eps^2), and see what operator
appears at each order.

Rick's Lemma 1:
    p_2(Y) . e_r(X) = A0(q,t) * e_(r,1,1)  +  A1(q,t) * e_(r,2)
                     + A2(q,t) * e_(r+1,1)  +  tau_r(q,t) * e_(r+2)

    A0 = 1/q^3
    A1 = -(q*t - q + t + 1)/q^3
    A2 = (q^2 - 1)/q^3
    tau_r * q^3 = A + B t^r + C t^(2r)
      A = -(q^2-1)(q-t-1)/(t^2-1)
      B = t(q^2-1)(q-t)/(t-1)
      C = -q t^3 (q^2-1)/(t^2-1)
"""
import sympy as sp

q, t, alpha, eps = sp.symbols('q t alpha epsilon', positive=True)

def A_coef(): return -(q**2 - 1)*(q - t - 1)/(t**2 - 1)
def B_coef(): return t*(q**2 - 1)*(q - t)/(t - 1)
def C_coef(): return -q * t**3 * (q**2 - 1)/(t**2 - 1)

def tau(r):
    return (A_coef() + B_coef() * t**r + C_coef() * t**(2*r))/q**3

def A0(): return 1/q**3
def A1(): return -(q*t - q + t + 1)/q**3
def A2(): return (q**2 - 1)/q**3

# Substitute q = t^alpha, expand around t = 1 + eps.
def jack_limit_coef(coef, order):
    """Take a rational function of q,t, substitute q = (1+eps)^alpha and
    t = 1 + eps, and return the coefficient of eps^order."""
    expr = coef.subs(q, (1 + eps)**alpha).subs(t, 1 + eps)
    # Series expansion around eps = 0
    ser = sp.series(expr, eps, 0, order + 1).removeO()
    poly = sp.Poly(sp.expand(ser), eps)
    return sp.simplify(poly.coeff_monomial(eps**order))

def report_coef(name, coef, r_val=None):
    if r_val is not None:
        coef = coef.subs('r', r_val) if hasattr(coef, 'subs') else coef
    print(f"  {name}:")
    for order in range(3):
        c = jack_limit_coef(coef, order)
        print(f"    eps^{order}: {sp.simplify(c)}")

print("="*72)
print("Jack limit expansion: q = t^alpha, t = 1 + eps")
print("="*72)

print("\n--- Coefficient of e_(r,1,1):  A0 = 1/q^3 ---")
for order in range(3):
    print(f"  eps^{order}: {sp.simplify(jack_limit_coef(A0(), order))}")

print("\n--- Coefficient of e_(r,2):  A1 = -(qt-q+t+1)/q^3 ---")
for order in range(3):
    print(f"  eps^{order}: {sp.simplify(jack_limit_coef(A1(), order))}")

print("\n--- Coefficient of e_(r+1,1):  A2 = (q^2-1)/q^3 ---")
for order in range(3):
    print(f"  eps^{order}: {sp.simplify(jack_limit_coef(A2(), order))}")

# For tau_r we do it separately by r_val since t^r must be expanded too.
for r_val in [2, 3, 4]:
    print(f"\n--- Coefficient of e_(r+2)  (r={r_val}):  tau_r ---")
    tau_r = tau(r_val)
    for order in range(3):
        # substitute q, t, then extract
        expr = tau_r.subs(q, (1 + eps)**alpha).subs(t, 1 + eps)
        ser = sp.series(expr, eps, 0, order + 1).removeO()
        poly = sp.Poly(sp.expand(ser), eps)
        c = poly.coeff_monomial(eps**order)
        print(f"    eps^{order}: {sp.simplify(c)}")

