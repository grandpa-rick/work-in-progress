"""
Day 200: Analyze the Baxterised structure of the Day 200 closed form.

We have:
  tau_r(q,t) * q^3 = A + B t^r + C t^(2r)
    A = -(q^2-1)(q-t-1)/(t^2-1)
    B =  t(q^2-1)(q-t)/(t-1)
    C = -q t^3 (q^2-1)/(t^2-1)

Clio's Baxter atom:  B_r := (q-t) + t^(r-1)(1-qt).

Question:  does tau_r * q^3 factor through B_r somehow?
"""
import sympy as sp

q, t, r = sp.symbols('q t r')


def tau_r_qcube(rval):
    return (-(q**2 - 1) * (q - t - 1) / (t**2 - 1)
            + t * (q**2 - 1) * (q - t) / (t - 1) * t**rval
            + (-q * t**3 * (q**2 - 1) / (t**2 - 1)) * t**(2*rval))


def Baxter_r(rval):
    return (q - t) + t**(rval - 1) * (1 - q*t)


# ------------------------------------------------------------------
# Test 1: Does B_r divide tau_r * q^3 as a polynomial in q?
# ------------------------------------------------------------------
def test_divides_by_baxter():
    print("=" * 72)
    print("Does B_r = (q-t) + t^(r-1)(1-qt) divide tau_r*q^3 in Q(t)[q]?")
    print("=" * 72)
    for rval in [2, 3, 4, 5, 6, 7]:
        tau_val = sp.expand(tau_r_qcube(rval))
        B_val = sp.expand(Baxter_r(rval))
        # Both are polynomials in q, coefficients in Q(t)
        tau_poly = sp.Poly(tau_val, q)
        B_poly = sp.Poly(B_val, q)
        quo, rem = sp.div(tau_poly, B_poly)
        rem_expr = sp.expand(rem.as_expr())
        rem_simplified = sp.simplify(rem_expr)
        divides = (rem_simplified == 0)
        print(f"  r = {rval}: divides? {divides}")
        if divides:
            print(f"    quotient (in q) = {sp.factor(quo.as_expr())}")


# ------------------------------------------------------------------
# Test 2: 2-term Baxterised ansatz
#   tau_r * q^3 = M(q,t) * B_r(q,t) + N(q,t) * t^(2r) ??
# ------------------------------------------------------------------
def test_baxter_plus_correction():
    print()
    print("=" * 72)
    print("Ansatz: tau_r*q^3 = M(q,t) * B_r + N(q,t) t^(2r)  (2-term).")
    print("Fit on r = 2, 3 and predict r = 4, 5.")
    print("=" * 72)
    M, N = sp.symbols('M N')
    eqns = []
    for rval in [2, 3]:
        tau_val = sp.expand(tau_r_qcube(rval))
        B_val = sp.expand(Baxter_r(rval))
        eqns.append(M * B_val + N * t**(2*rval) - tau_val)
    sol = sp.solve(eqns, [M, N], dict=True)
    if not sol:
        print("  No solution.")
        return
    sol = sol[0]
    Mv = sp.factor(sp.simplify(sol[M]))
    Nv = sp.factor(sp.simplify(sol[N]))
    print(f"  M(q,t) = {Mv}")
    print(f"  N(q,t) = {Nv}")
    for rval in [4, 5, 6]:
        tau_val = sp.expand(tau_r_qcube(rval))
        B_val = sp.expand(Baxter_r(rval))
        pred = sp.simplify(Mv * B_val + Nv * t**(2*rval))
        diff = sp.simplify(sp.expand(pred - tau_val))
        print(f"  r = {rval}: diff = {diff}  MATCH? {diff == 0}")


# ------------------------------------------------------------------
# Test 3: 3-term Baxterised ansatz
#   tau_r * q^3 = A' + M(q,t) * B_r + N(q,t) t^(2r) ??
# ------------------------------------------------------------------
def test_three_term_baxter():
    print()
    print("=" * 72)
    print("Ansatz: tau_r*q^3 = A0(q,t) + M(q,t) * B_r + N(q,t) t^(2r).")
    print("Fit on r = 2, 3, 4 and predict r = 5.")
    print("=" * 72)
    A0, M, N = sp.symbols('A0 M N')
    eqns = []
    for rval in [2, 3, 4]:
        tau_val = sp.expand(tau_r_qcube(rval))
        B_val = sp.expand(Baxter_r(rval))
        eqns.append(A0 + M * B_val + N * t**(2*rval) - tau_val)
    sol = sp.solve(eqns, [A0, M, N], dict=True)
    if not sol:
        print("  No solution (system inconsistent).")
        return
    sol = sol[0]
    A0v = sp.factor(sp.simplify(sol[A0]))
    Mv = sp.factor(sp.simplify(sol[M]))
    Nv = sp.factor(sp.simplify(sol[N]))
    print(f"  A0(q,t) = {A0v}")
    print(f"  M(q,t)  = {Mv}")
    print(f"  N(q,t)  = {Nv}")

    for rval in [5, 6, 7]:
        tau_val = sp.expand(tau_r_qcube(rval))
        B_val = sp.expand(Baxter_r(rval))
        pred = sp.simplify(A0v + Mv * B_val + Nv * t**(2*rval))
        diff = sp.simplify(sp.expand(pred - tau_val))
        print(f"  r = {rval}: diff = {diff}  MATCH? {diff == 0}")


# ------------------------------------------------------------------
# Test 4: Try tau_r * q^3 = M(q,t) * B_r + N(q,t) * B_r^2 (double-Baxter!)
# ------------------------------------------------------------------
def test_double_baxter():
    print()
    print("=" * 72)
    print("Ansatz: tau_r*q^3 = M(q,t) * B_r + N(q,t) * B_r^2.")
    print("=" * 72)
    M, N = sp.symbols('M N')
    eqns = []
    for rval in [2, 3]:
        tau_val = sp.expand(tau_r_qcube(rval))
        B_val = sp.expand(Baxter_r(rval))
        eqns.append(M * B_val + N * B_val**2 - tau_val)
    sol = sp.solve(eqns, [M, N], dict=True)
    if not sol:
        print("  no solution")
        return
    sol = sol[0]
    Mv = sp.factor(sp.simplify(sol[M]))
    Nv = sp.factor(sp.simplify(sol[N]))
    print(f"  M(q,t) = {Mv}")
    print(f"  N(q,t) = {Nv}")
    for rval in [4, 5]:
        tau_val = sp.expand(tau_r_qcube(rval))
        B_val = sp.expand(Baxter_r(rval))
        pred = sp.simplify(Mv * B_val + Nv * B_val**2)
        diff = sp.simplify(sp.expand(pred - tau_val))
        print(f"  r = {rval}: diff = {diff}  MATCH? {diff == 0}")


# ------------------------------------------------------------------
# Test 5: (1-t)*tau_r*q^3 vs B'_r := q [r]_t - t [r-2]_t  (Clio's LHS)
# ------------------------------------------------------------------
def test_bprime_atom():
    print()
    print("=" * 72)
    print("Alternative Baxter atom (Clio's LHS):  B'_r := q [r]_t - t [r-2]_t.")
    print("Test: tau_r*q^3 = M(q,t) * B'_r + N(q,t) * (B'_r)^2 ?")
    print("=" * 72)
    def Bprime(rval):
        return q * sum(t**i for i in range(rval)) - t * sum(t**i for i in range(rval - 2))
    M, N = sp.symbols('M N')
    eqns = []
    for rval in [2, 3]:
        tau_val = sp.expand(tau_r_qcube(rval))
        Bp = sp.expand(Bprime(rval))
        eqns.append(M * Bp + N * Bp**2 - tau_val)
    sol = sp.solve(eqns, [M, N], dict=True)
    if not sol:
        print("  no solution")
        return
    sol = sol[0]
    Mv = sp.factor(sp.simplify(sol[M]))
    Nv = sp.factor(sp.simplify(sol[N]))
    print(f"  M(q,t) = {Mv}")
    print(f"  N(q,t) = {Nv}")
    for rval in [4, 5]:
        tau_val = sp.expand(tau_r_qcube(rval))
        Bp = sp.expand(Bprime(rval))
        pred = sp.simplify(Mv * Bp + Nv * Bp**2)
        diff = sp.simplify(sp.expand(pred - tau_val))
        print(f"  r = {rval}: diff = {diff}  MATCH? {diff == 0}")


if __name__ == "__main__":
    test_divides_by_baxter()
    test_baxter_plus_correction()
    test_three_term_baxter()
    test_double_baxter()
    test_bprime_atom()
