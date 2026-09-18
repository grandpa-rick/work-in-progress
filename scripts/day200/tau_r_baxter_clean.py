"""
Day 200 refinement: try to expose the Baxter atom more cleanly.

Two closed forms so far:

(1) THREE-MONOMIAL form (crispest, `computed`):
    tau_r(q,t) * q^3 = A + B t^r + C t^(2r)
      A = -(q^2-1)(q-t-1)/(t^2-1)
      B =  t(q^2-1)(q-t)/(t-1)
      C = -q t^3 (q^2-1)/(t^2-1)

(2) THREE-TERM BAXTER form (uses Clio's atom B_r):
    tau_r(q,t) * q^3 = A0 + M B_r + N t^(2r)
      B_r = (q-t) + t^(r-1)(1-qt)
    verifies at r=5,6,7 fitting on r=2,3,4.

I'll rewrite (1) directly in terms of B_r to expose the Baxter structure.

We have
  A + B t^r + C t^(2r) = ?
Using   B_r = (q-t) + t^(r-1)(1-qt) = (q-t) + t^(r-1) - q t^r

Multiply out:  q * B_r = q(q-t) + q t^(r-1) - q^2 t^r
Try to rewrite B t^r as combination:
  B t^r =  t(q^2-1)(q-t)/(t-1) * t^r
        =  (q^2-1)(q-t)/(t-1) * t^(r+1)

Hmm.  Let's try:  what is (something) * B_r * (t^r-piece)?

Actually the cleanest observation:
  tau_r * q^3 * (t-1) = -(q^2-1) * [ (q-t-1)/(t+1) - t(q-t)(t+1)*t^r/(t+1) + qt^3*t^(2r)/(t+1)]

Let me just factor out -(q^2-1)/((t-1)(t+1)) globally.
"""
import sympy as sp

q, t = sp.symbols('q t')


def A_coef():
    return -(q**2 - 1) * (q - t - 1) / (t**2 - 1)

def B_coef():
    return t * (q**2 - 1) * (q - t) / (t - 1)

def C_coef():
    return -q * t**3 * (q**2 - 1) / (t**2 - 1)


def tau_closed(rval):
    return (A_coef() + B_coef() * t**rval + C_coef() * t**(2*rval)) / q**3


def Baxter_r(rval):
    return (q - t) + t**(rval - 1) * (1 - q*t)


def try_baxter_rewrite():
    """
    Compute (tau_r * q^3) - K1(q,t) * B_r(q,t) - K2(q,t) * B_r(q,t) * t^r
    for various trial K1, K2; look for a clean residual.
    """
    print("=" * 72)
    print("Try: tau_r*q^3 =?= K1(q,t) * B_r + K2(q,t) * B_r * t^r + Residual(t^(2r))?")
    print("=" * 72)
    # Fit K1, K2 as symbols, on r = 2, 3
    K1, K2, R2 = sp.symbols('K1 K2 R2')  # R2 = t^(2r) coefficient
    eqns = []
    for rval in [2, 3, 4]:
        tau_val = sp.expand(tau_closed(rval) * q**3)
        B_val = sp.expand(Baxter_r(rval))
        eqns.append(K1 * B_val + K2 * B_val * t**rval + R2 * t**(2*rval) - tau_val)
    sol = sp.solve(eqns, [K1, K2, R2], dict=True)
    if not sol:
        print("  no solution")
        return
    sol = sol[0]
    K1v = sp.factor(sp.simplify(sol[K1]))
    K2v = sp.factor(sp.simplify(sol[K2]))
    R2v = sp.factor(sp.simplify(sol[R2]))
    print(f"  K1 = {K1v}")
    print(f"  K2 = {K2v}")
    print(f"  R2 = {R2v}")

    for rval in [5, 6, 7]:
        tau_val = sp.expand(tau_closed(rval) * q**3)
        B_val = sp.expand(Baxter_r(rval))
        pred = sp.simplify(K1v * B_val + K2v * B_val * t**rval + R2v * t**(2*rval))
        diff = sp.simplify(sp.expand(pred - tau_val))
        print(f"  r = {rval}: diff = {diff}  {'OK' if diff == 0 else 'FAIL'}")


def try_pure_baxter():
    """
    Ansatz: tau_r * q^3 = K1 * B_r + K2 * B_r * t^r (no t^(2r) residual).
    """
    print()
    print("=" * 72)
    print("Try: tau_r*q^3 =?= K1(q,t) * B_r + K2(q,t) * B_r * t^r  (pure double-Baxter)")
    print("=" * 72)
    K1, K2 = sp.symbols('K1 K2')
    eqns = []
    for rval in [2, 3]:
        tau_val = sp.expand(tau_closed(rval) * q**3)
        B_val = sp.expand(Baxter_r(rval))
        eqns.append(K1 * B_val + K2 * B_val * t**rval - tau_val)
    sol = sp.solve(eqns, [K1, K2], dict=True)
    if not sol:
        print("  no solution")
        return
    sol = sol[0]
    K1v = sp.factor(sp.simplify(sol[K1]))
    K2v = sp.factor(sp.simplify(sol[K2]))
    print(f"  K1 = {K1v}")
    print(f"  K2 = {K2v}")
    for rval in [4, 5, 6]:
        tau_val = sp.expand(tau_closed(rval) * q**3)
        B_val = sp.expand(Baxter_r(rval))
        pred = sp.simplify(K1v * B_val + K2v * B_val * t**rval)
        diff = sp.simplify(sp.expand(pred - tau_val))
        print(f"  r = {rval}: diff = {diff}  {'OK' if diff == 0 else 'FAIL'}")


def try_prod_ansatz():
    """
    Ansatz: tau_r * q^3 = -(q^2-1) * X_r(q,t) * B_r(q,t) for some rational X_r.
    Compute X_r := tau_r * q^3 / [-(q^2-1) * B_r] for each r, and inspect
    for a simple form.
    """
    print()
    print("=" * 72)
    print("Compute X_r := tau_r*q^3 / [-(q^2-1) * B_r]:")
    print("=" * 72)
    for rval in [2, 3, 4, 5, 6]:
        tau_val = sp.expand(tau_closed(rval) * q**3)
        B_val = sp.expand(Baxter_r(rval))
        X = sp.cancel(sp.together(tau_val / (-(q**2 - 1) * B_val)))
        print(f"  r = {rval}: X_r = {sp.factor(X)}")


def try_bar_form():
    """
    Note that the closed form can be rewritten using B_r as follows.
    Compute:  qtau_r := tau_r * q^3 * (t-1) / [-(q^2-1)]

    Then explicit substitution gives:
      qtau_r = (q-t-1) - t(q-t)(t+1)*t^r/(t+1) ... hmm.

    Let me just print the (t-1) * tau_r * q^3 / [-(q^2-1)] expanded and factor.
    """
    print()
    print("=" * 72)
    print("Structure: (t^2-1) * tau_r * q^3 / [-(q^2-1)] = polynomial in q, t:")
    print("=" * 72)
    for rval in [2, 3, 4, 5, 6]:
        expr = sp.expand(tau_closed(rval) * q**3 * (t**2 - 1) / (-(q**2 - 1)))
        print(f"  r = {rval}: {expr}")
        # Try to express as sum of Baxter atoms times t^k:
        # Formal Baxter atom (q-t) + t^(r-1)(1-qt) = q-t + t^(r-1) - qt^r
        Br = sp.expand(Baxter_r(rval))
        # try dividing (t^2-1)*tau_r*q^3/(-(q^2-1)) by (t+1) * B_r
        quot = sp.cancel(sp.together(expr / ((t + 1) * Br)))
        print(f"    ratio / ((t+1) * B_r) = {sp.factor(quot)}")


if __name__ == "__main__":
    try_baxter_rewrite()
    try_pure_baxter()
    try_prod_ansatz()
    try_bar_form()
