"""
Day 201: Verify all conjectured closed forms against r=5 SymPy data.
"""
import sympy as sp
q, t = sp.symbols('q t')

def q_int(n):
    return sum(q**i for i in range(n))

# k-uniform r-indep formulas at k=3:
def c_r_111(): return 1 / q**6
def c_r_21(): return -(q**2*t - q**2 + q*t - q + t + 2) / q**6
def c_rp1_11(): return (q**3 - 1) / q**6
def c_r_3(): return (q**3*t**3 - q**3*t**2 - q**3*t + q**3 + q**2*t**3 - q**2 + q*t**3 - q + t**2 + t + 1) / q**6
def c_rp1_2(): return -(q**3 - 1)*(q*t**2 + t - q + 1) / q**6

# r-dep formulas:
def c_rp2_1(r_val): return -(q**3 - 1)*(q**2 * t**(r_val+1) - q**2 - q * t**(r_val+2) + q*t + 1) / q**6

# Baxter-4 formula for c_(r+3):
def baxter4_A(j, r_val):
    """Baxter-4 fit from r=1..4."""
    # Copied from fit_baxter4 output
    if j == 0:
        return -(q - 1)*(q**2 + q + 1)*(q**3 - q**2*t**2 - q**2*t - q**2 + q*t**3 - q + t**2 + t + 1) / ((t - 1)*(t**2 + t + 1))
    if j == 1:
        return -t*(q - 1)*(q - t)*(q**2 + q + 1)*(-q**2 + q*t**2 - q*t + q + t) / (t - 1)
    if j == 2:
        return q**2*t**3*(-q**4 + q**3*t + q - t) / (t - 1)
    if j == 3:
        return q**3*t**6*(q**3 - 1) / (t**3 - 1)

def c_rp3(r_val):
    return sum(baxter4_A(j, r_val) * t**(j * r_val) for j in range(4)) / q**6

# r=5 empirical data:
r5_data = {
    (5, 1, 1, 1): sp.Integer(1) / q**6,
    (5, 2, 1): -(q**2*t - q**2 + q*t - q + t + 2) / q**6,
    (6, 1, 1): (q - 1)*(q**2 + q + 1) / q**6,
    (5, 3): (q**3*t**3 - q**3*t**2 - q**3*t + q**3 + q**2*t**3 - q**2 + q*t**3 - q + t**2 + t + 1) / q**6,
    (6, 2): -(q - 1)*(t + 1)*(q**2 + q + 1)*(q*t - q + 1) / q**6,
    (7, 1): -(q - 1)*(q**2 + q + 1)*(q**2*t**6 - q**2 - q*t**7 + q*t + 1) / q**6,
    (8,): (q - 1)*(t + 1)*(t**2 + 1)*(t**4 + 1)*(q**2 + q + 1)*(q**3*t**11 - q**3*t**10 + q**3*t**8 - q**3*t**7 - q**3*t**4 + q**3*t**3 - q**3*t + q**3 + q**2*t**6 - q**2 + q*t - q + 1) / q**6,
}

r_val = 5
predictions = {
    (r_val, 1, 1, 1): c_r_111(),
    (r_val, 2, 1): c_r_21(),
    (r_val + 1, 1, 1): c_rp1_11(),
    (r_val, 3): c_r_3(),
    (r_val + 1, 2): c_rp1_2(),
    (r_val + 2, 1): c_rp2_1(r_val),
    (r_val + 3,): c_rp3(r_val),
}

print("=" * 72)
print("Day 201: Verify all p_3(Y).e_5 closed forms against SymPy at r=5, m=8")
print("=" * 72)
print()
ok = 0
fail = 0
for mu, pred in predictions.items():
    emp = r5_data.get(mu)
    if emp is None:
        print(f"  mu = {mu}: NO EMPIRICAL DATA")
        continue
    diff = sp.simplify(sp.expand(pred - emp))
    status = "OK" if diff == 0 else f"FAIL (diff = {diff})"
    print(f"  mu = {mu}: {status}")
    if diff == 0:
        ok += 1
    else:
        fail += 1

print()
print(f"Total: {ok} PASS, {fail} FAIL out of {ok+fail}")
