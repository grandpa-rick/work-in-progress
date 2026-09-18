"""
Day 201: Verify Baxter-3 fit for c_(r+3) against r=1 data (independent check).

Baxter-3 fitted from r=2, 3, 4 gives:
  A_0 = (q-1)(q^2+q+1)(q^3 t^{13} - q^3 t^{12} + q^3 t^{10} - q^3 t^9 + q^3 t^7
       - q^3 t^6 + q^3 t^4 - q^3 t^3 + q^3 t - q^3 + q^2 - qt + q - 1)/(t-1)
  A_1 = -t(q-1)(q^2+q+1)(q^3 t^{10} - q^3 + q^2 t^2 + q^2 - qt^3 + qt^2 - t^2)/(t-1)
  A_2 = q^2 t^3 (q-1)(q^2+q+1)(qt^5 - q + t)/(t-1)

Prediction at r=1: A_0 + A_1 t + A_2 t^2.
Empirical r=1 (from log): (q-1)(t+1)(t^2+1)(q^2+q+1)(qt-q+1)(q^2 t^2 - q^2 + 1).
"""
import sympy as sp

q, t = sp.symbols('q t')

A_0 = (q-1)*(q**2+q+1) * (
    q**3 * t**13 - q**3 * t**12 + q**3 * t**10 - q**3 * t**9
    + q**3 * t**7 - q**3 * t**6 + q**3 * t**4 - q**3 * t**3
    + q**3 * t - q**3 + q**2 - q*t + q - 1
) / (t - 1)

A_1 = -t * (q-1) * (q**2+q+1) * (
    q**3 * t**10 - q**3 + q**2 * t**2 + q**2 - q*t**3 + q*t**2 - t**2
) / (t - 1)

A_2 = q**2 * t**3 * (q-1) * (q**2+q+1) * (q * t**5 - q + t) / (t - 1)

def baxter3_predict(r_val):
    return A_0 + A_1 * t**r_val + A_2 * t**(2*r_val)

# r=1 empirical:
r1_empirical = (q-1) * (t+1) * (t**2+1) * (q**2+q+1) * (q*t - q + 1) * (q**2 * t**2 - q**2 + 1)

pred_r1 = baxter3_predict(1)
diff_r1 = sp.simplify(sp.expand(pred_r1 - r1_empirical))
print(f"Baxter-3 prediction at r=1 vs empirical: diff = {diff_r1}")

if diff_r1 == 0:
    print("MATCH at r=1: Baxter-3 fit for c_(r+3) is VERIFIED at r=1.")
else:
    print("MISMATCH at r=1: Baxter-3 fit needs to be corrected.")

# Sanity: verify at r=2, 3, 4 (should trivially match since we fit these):
r_data = {
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
print()
print("Sanity checks (should be zero):")
for r_val in [2, 3, 4]:
    diff = sp.simplify(sp.expand(baxter3_predict(r_val) - r_data[r_val]))
    print(f"  r = {r_val}: diff = {diff}")
