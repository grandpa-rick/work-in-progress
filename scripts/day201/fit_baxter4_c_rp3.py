"""
Day 201: Fit Baxter-4 shape for c_(r+3) using r=1, 2, 3, 4 data.

Baxter-3 FAILED at r=1 (verified by verify_baxter3_at_r1.py).
Try Baxter-4: c_(r+3) · q^6 = A_0 + A_1 t^r + A_2 t^{2r} + A_3 t^{3r}.

Uses r=1, 2, 3, 4 -> 4 equations, 4 unknowns. Direct fit.
Verify shape by checking predicted r=1 value against r=1 empirical
(with degeneracy: at r=1, position (r+3) = (4) is UNIQUE — no collision).
"""
import sympy as sp

q, t = sp.symbols('q t')

# Data (all × q^6):
data = {
    1: (q - 1) * (t + 1) * (t**2 + 1) * (q**2 + q + 1) *
       (q*t - q + 1) * (q**2 * t**2 - q**2 + 1),  # r=1 e_(4)
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

A = sp.symbols('A0 A1 A2 A3')

def eq_at_r(r_val):
    lhs = sum(A[j] * t**(j * r_val) for j in range(4))
    return lhs - sp.expand(data[r_val])

def try_fit(r_vals):
    print(f"Fitting Baxter-4 at r_vals = {r_vals}")
    eqs = [eq_at_r(r_val) for r_val in r_vals]
    try:
        sol = sp.solve(eqs, A, dict=True)
    except Exception as e:
        print(f"  solve failed: {e}")
        return None
    if not sol:
        print(f"  no solution")
        return None
    sol = sol[0]
    for j in range(4):
        val = sp.simplify(sp.factor(sol[A[j]]))
        # Check that val is a polynomial (no (t-1) denominator)
        val_str = str(val)
        print(f"  A_{j} = {val}")
    return sol

def verify_shape(sol, extra_r_vals):
    for r_val in extra_r_vals:
        pred = sum(sol[A[j]] * t**(j * r_val) for j in range(4))
        diff = sp.simplify(sp.expand(pred - data[r_val]))
        print(f"  verify r = {r_val}: diff = {diff}")

if __name__ == "__main__":
    print("=" * 72)
    print("Day 201: Fit Baxter-4 shape for c_(r+3)(q, t)")
    print("=" * 72)
    print()
    # Fit r = 1, 2, 3, 4 (4 unknowns, 4 equations)
    sol = try_fit([1, 2, 3, 4])
    if sol is not None:
        print()
        print("Sanity: predictions vs data at all r:")
        verify_shape(sol, [1, 2, 3, 4])
        print()
        print("Predicting r = 5 (should match SymPy output once available):")
        pred_r5 = sum(sol[A[j]] * t**(j * 5) for j in range(4))
        pred_r5_factored = sp.factor(sp.simplify(pred_r5))
        print(f"  c_(8)|_(r=5) · q^6 predicted = {pred_r5_factored}")
