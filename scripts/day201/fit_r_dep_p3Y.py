"""
Day 201: Fit r-dependent closed forms for p_3(Y) . e_r.

Empirical data at r = 2, 3, 4 shows 5 r-INDEPENDENT coefficients and
2 r-DEPENDENT coefficients:
  r-indep: e_(r, 1, 1, 1), e_(r, 2, 1), e_(r+1, 1, 1), e_(r, 3), e_(r+1, 2).
  r-dep:   e_(r+2, 1), e_(r+3).

r-indep closed forms (q^{-6} normalization):
  c_(r, 1^3)     = 1
  c_(r, 2, 1)    = -(q^2 t - q^2 + q t - q + t + 2)
  c_(r+1, 1, 1)  = q^3 - 1
  c_(r, 3)       = q^3 t^3 - q^3 t^2 - q^3 t + q^3 + q^2 t^3 - q^2 + q t^3 - q + t^2 + t + 1
  c_(r+1, 2)     = -(q^3 - 1)(q t^2 + t - q + 1)

r-dep closed forms (conjectural):
  c_(r+2, 1) · q^6 = -(q^3 - 1)(q^2 t^{r+1} - q^2 - q t^{r+2} + q t + 1)   [BAXTER-2 shape]
  c_(r+3) · q^6   = A + B t^r + C t^{2r} + D t^{3r}?  [BAXTER-4 shape, unknown A,B,C,D]

Fit Baxter-4 using r = 2, 3, 4 data (3 eqns, 4 unknowns => underdetermined).
Attempt Baxter-3 (3 unknowns): A + B t^r + C t^{2r}?
"""
import sympy as sp

q, t = sp.symbols('q t')

# c_(r+3) data at r = 2, 3, 4 (times q^6):
data = {
    2: (q**3 - 1) * (t**4 + t**3 + t**2 + t + 1) *
       (q**3 * t**5 - q**3 * t**4 - q**3 * t + q**3
        + q**2 * t**3 - q**2 + q*t - q + 1),
    3: (q**3 - 1) * (t + 1) * (t**2 - t + 1) *
       (q**3 * t**9 - q**3 * t**5 - q**3 * t**4 + q**3
        + q**2 * t**6 + q**2 * t**5 + q**2 * t**4
        - q**2 * t**2 - q**2 * t - q**2
        + q*t**3 - q + t**2 + t + 1),
    4: (q**3 - 1) * (t**6 + t**5 + t**4 + t**3 + t**2 + t + 1) *
       (q**3 * t**9 - q**3 * t**8 + q**3 * t**6 - q**3 * t**5
        - q**3 * t**4 + q**3 * t**3 - q**3 * t + q**3
        + q**2 * t**5 - q**2 + q*t - q + 1),
}

# Also verify c_(r+2, 1) Baxter-2 conjecture:
def c_rp2_1_conj(r_val):
    return -(q**3 - 1) * (q**2 * t**(r_val+1) - q**2 - q * t**(r_val+2) + q*t + 1)

def c_rp2_1_data(r_val):
    if r_val == 2:
        return -(q**3 - 1) * (q**2 * t**3 - q**2 - q * t**4 + q*t + 1)
    if r_val == 3:
        return -(q**3 - 1) * (q**2 * t**4 - q**2 - q * t**5 + q*t + 1)
    if r_val == 4:
        return -(q**3 - 1) * (q**2 * t**5 - q**2 - q * t**6 + q*t + 1)

def verify_baxter2_c_rp2_1():
    print("Verifying c_(r+2, 1) · q^6 = -(q^3-1)(q^2 t^{r+1} - q^2 - q t^{r+2} + qt + 1)")
    for r_val in [2, 3, 4]:
        diff = sp.simplify(c_rp2_1_conj(r_val) - c_rp2_1_data(r_val))
        print(f"  r = {r_val}: diff = {diff}")

def try_baxter_fit(order):
    """Fit c_(r+3) · q^6 = sum_{j=0}^{order-1} A_j(q,t) t^{jr} using data
    at r = 2, 3, ..., 2 + order - 1.  Return the fitted A_j's."""
    if order > 3:
        print(f"Baxter-{order} fit needs at least {order} data points; we have 3 (r=2,3,4).")
        return None
    r_vals = [2, 3, 4][:order]
    A = sp.symbols(f'A0:{order}')
    eqs = []
    for r_val in r_vals:
        expr = sum(A[j] * t**(j * r_val) for j in range(order))
        eqs.append(expr - sp.expand(data[r_val]))
    # Solve
    try:
        sol = sp.solve(eqs, A, dict=True)
    except Exception as e:
        print(f"solve failed: {e}")
        return None
    if not sol:
        print(f"Baxter-{order} fit: no solution")
        return None
    sol = sol[0]
    print(f"Baxter-{order} fit:")
    for j in range(order):
        val = sp.simplify(sp.factor(sol[A[j]]))
        print(f"  A_{j} = {val}")
    return sol

def verify_fit(sol, order):
    """Given fit at r = 2..order+1, verify at higher r (if available)."""
    if sol is None:
        return
    for r_val in [2, 3, 4]:
        A = sp.symbols(f'A0:{order}')
        pred = sum(sol[A[j]] * t**(j * r_val) for j in range(order))
        diff = sp.simplify(sp.expand(pred - data[r_val]))
        print(f"  verify r = {r_val}: diff = {diff}")

if __name__ == "__main__":
    print("=" * 72)
    print("Day 201: r-dependent closed-form fits for p_3(Y) . e_r")
    print("=" * 72)
    print()
    verify_baxter2_c_rp2_1()
    print()
    print("Trying Baxter-3 fit for c_(r+3):")
    sol3 = try_baxter_fit(3)
    if sol3:
        verify_fit(sol3, 3)
    print()
    print("Baxter-3 fit will consume all r = 2, 3, 4 data. For Baxter-4 we need r = 5.")
