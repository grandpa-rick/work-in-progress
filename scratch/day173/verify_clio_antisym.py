"""Day 173: Independent verification of Clio's Xi-free antisymmetric claim.

Claim:
    R^{(-1)}_n  =  (1/(2c)) * [deg_{u1,u2} = n-1] ( [T^n] log( F_c / F_{-c} ) )
for any c != 0, where F_c := F_P|_{u_3 = c}.

Rick's ground-truth R^{(-1)}_n (Day 167, def line 36):
    R^{(-1)}_n := ∂_{u_3} X^{(0)}_n |_{u_3 = 0}
    where X^{(0)}_n = degree-n homogeneous layer of [T^n] log F_P (u-degree = n).
    Equivalently: R^{(-1)}_n = coeff of u_3^1 in X^{(0)}_n (a poly in u_1, u_2 of
    total degree n-1). This matches step6_extended_verify.py's
    R_m1_n = coeff_u3(X0_layer[n], 1) formulation.

We verify Clio's identity for c in {1, 2, -1, 3} at n = 2..N. Exact Fractions.

Note: (1/(2c)) log(F_c/F_{-c}) is manifestly ODD in c since log(F_{-c}/F_c) =
-log(F_c/F_{-c}). So c and -c should give the same RHS; we check both signs
for redundancy (c=1 vs c=-1 must agree; c=2 must also equal the c=1 value).
"""
import sys, time
sys.path.insert(0, '/home/agent/projects/scratch/day152')
from lib import FP_coeffs, pmul, padd, psub, pscal, pconst, sdiv
from fractions import Fraction as Fr

# ------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------
def restrict_u3(p, val):
    """Set u_3 := val in the poly-dict p."""
    val = Fr(val)
    out = {}
    for (i, j, k), c in p.items():
        contribution = c * (val ** k)
        key = (i, j, 0)
        w = out.get(key, Fr(0)) + contribution
        if w:
            out[key] = w
        elif key in out:
            del out[key]
    return {k: v for k, v in out.items() if v}

def slog_series(f, N):
    """log of series f with f[0]=1, coefficients are u-poly dicts."""
    g = [{} for _ in range(N+1)]
    for n in range(1, N+1):
        s = pscal(f[n], n)
        for k in range(1, n):
            s = psub(s, pscal(pmul(g[k], f[n-k]), k))
        g[n] = pscal(s, Fr(1, n))
    return g

def extract_deg(p, d):
    """Full 3-var poly-dict: keep monomials with sum(exp) == d."""
    return {k: v for k, v in p.items() if sum(k) == d}

def extract_2var_deg(p, d):
    """2-var poly-dict (u_3 exponent = 0): keep i+j == d."""
    return {(i, j, 0): v for (i, j, k), v in p.items() if k == 0 and (i+j) == d}

def coeff_u3(p, k):
    """Coefficient of u_3^k in p, returned as 2-var dict (u_3 exp cleared)."""
    return {(i, j, 0): v for (i, j, kk), v in p.items() if kk == k}

def dict_scal(p, c):
    c = Fr(c)
    if c == 0:
        return {}
    return {k: v * c for k, v in p.items()}

def dict_diff(a, b):
    keys = set(a) | set(b)
    d = {k: a.get(k, Fr(0)) - b.get(k, Fr(0)) for k in keys}
    return {k: v for k, v in d.items() if v}

def fmt_poly(p, maxterms=4):
    if not p:
        return "0"
    items = sorted(p.items())
    if len(items) <= maxterms:
        return " + ".join(f"{v}*u1^{k[0]}*u2^{k[1]}" for k, v in items)
    head = items[:maxterms]
    return " + ".join(f"{v}*u1^{k[0]}*u2^{k[1]}" for k, v in head) + f" + ... ({len(items)} terms)"

# ------------------------------------------------------------------
# Compute
# ------------------------------------------------------------------
N = 10           # verify n = 2..N
CS = [1, 2, -1, 3, Fr(1, 2)]

t0 = time.time()
print(f"Computing FP_coeffs to N={N}...")
FP = FP_coeffs(N)
print(f"  done in {time.time()-t0:.2f}s")

# Ground-truth R^{(-1)}_n via log F_P (3-var), layer decomposition.
t0 = time.time()
print(f"Computing 3-var log F_P to N={N}...")
logFP = slog_series(FP, N)
print(f"  done in {time.time()-t0:.2f}s")

X0_layer = [extract_deg(logFP[n], n) for n in range(N+1)]
Rm1_true = [coeff_u3(X0_layer[n], 1) for n in range(N+1)]  # R^{(-1)}_n

# ------------------------------------------------------------------
# Verify Clio's Xi-free identity for each c
# ------------------------------------------------------------------
print()
print("Ground-truth R^{(-1)}_n (Day 167 def) computed.")
print("Now testing Clio's claim for c in", CS)
print()

all_ok = True
fail_records = []

for c in CS:
    print(f"--- c = {c} ---")
    t0 = time.time()
    Fc  = [restrict_u3(FP[n], c)  for n in range(N+1)]
    Fmc = [restrict_u3(FP[n], -c) for n in range(N+1)]
    # log(F_c / F_{-c}) as series
    Ratio = sdiv(Fc, Fmc, N)
    log_ratio = slog_series(Ratio, N)
    print(f"  built log(F_c/F_{{-c}}) in {time.time()-t0:.2f}s")

    for n in range(2, N+1):
        deg_piece = extract_2var_deg(log_ratio[n], n-1)
        RHS = dict_scal(deg_piece, Fr(1, 2*c))
        diff = dict_diff(Rm1_true[n], RHS)
        ok = not diff
        tag = "PASS" if ok else "FAIL"
        print(f"    n={n}: {tag}")
        if not ok:
            all_ok = False
            fail_records.append((c, n, Rm1_true[n], RHS, diff))

print()
print("=" * 60)
if all_ok:
    print(f"ALL {(N-1)*len(CS)} instances PASS: c in {CS}, n = 2..{N}.")
    print("Clio's Xi-free antisymmetric identity is confirmed on the")
    print("tested range against Rick's ground-truth R^{(-1)}_n.")
else:
    print(f"FAILURES: {len(fail_records)}")
    for c, n, lhs, rhs, diff in fail_records:
        print(f"\n  c={c}, n={n}:")
        print(f"    LHS (R^{{-1}}_{n}):   {fmt_poly(lhs)}")
        print(f"    RHS (Clio / 2c):    {fmt_poly(rhs)}")
        print(f"    LHS - RHS:          {fmt_poly(diff)}")
print("=" * 60)
