#!/usr/bin/env python3
"""
Day 183 — Compute Witt-formula primitives p_k for the free graded Lie algebra
underlying U(L(a)) with dim U(L(a))_k = b_k.

  b_k = 3, 27, 417, 7851, 164124, ...
  a_k = INVERTi(b_k) = 3, 18, 282, 5268, ...

By Andrews-Gagnon-Gélinas-Schlums-Zabrocki (2505.06941), b_k = dim U(L(a))_k
where L(a) is the free graded Lie algebra on a_n generators in degree n.

Witt/PLog formula for graded Lie dim from generator counts:
    sum_n p_n t^n = sum_{d>=1} (mu(d)/d) log B(t^d),
where B(t) = 1 + sum b_k t^k = 1/(1 - A(t)).
"""

from fractions import Fraction
from sympy import symbols, Poly, log, series, Rational, sympify, prod, factorial
from sympy import Symbol, Integer

# ------------------------------------------------------------------
# Data
# ------------------------------------------------------------------

a = [0, 3, 18, 282, 5268, 109647, 2438928, 56758176, 1364824620,
     33643660620, 845633502606, 21590775239850, 558411335278644]

b = [1, 3, 27, 417, 7851, 164124, 3661389, 85384566, 2056373739,
     50751637140, 1276862920140, 32626363346505, 844375375808301]

N = 12  # want p_1 .. p_12

# ------------------------------------------------------------------
# Power series arithmetic over Q, truncated mod t^{N+1}
# ------------------------------------------------------------------

def ps_add(P, Q, n):
    return [P[i] + Q[i] for i in range(n+1)]

def ps_scalar(P, c, n):
    return [c * P[i] for i in range(n+1)]

def ps_mul(P, Q, n):
    R = [Fraction(0)] * (n+1)
    for i in range(n+1):
        if P[i] == 0:
            continue
        for j in range(n+1-i):
            R[i+j] += P[i]*Q[j]
    return R

def ps_log(P, n):
    """log(P) where P[0] == 1. Uses log(P)' = P'/P, so if L = log P,
       L[k] for k>=1 satisfies P[0]*k*L[k] = k*P[k] - sum_{j=1}^{k-1} j*L[j]*P[k-j].
       Since P[0]=1: k*L[k] = k*P[k] - sum_{j=1}^{k-1} j*L[j]*P[k-j].
    """
    assert P[0] == 1
    L = [Fraction(0)] * (n+1)
    for k in range(1, n+1):
        s = Fraction(k) * P[k]
        for j in range(1, k):
            s -= Fraction(j) * L[j] * P[k-j]
        L[k] = s / k
    return L

def ps_substitute_td(P, d, n):
    """Given P(t), return P(t^d) truncated to degree n."""
    Q = [Fraction(0)] * (n+1)
    for i in range(n+1):
        if i*d > n:
            break
        Q[i*d] = P[i]
    return Q

# Moebius function
def mobius(k):
    if k == 1:
        return 1
    # factor
    n = k
    result = 1
    p = 2
    while p*p <= n:
        if n % p == 0:
            n //= p
            if n % p == 0:
                return 0
            result = -result
        p += 1
    if n > 1:
        result = -result
    return result

# ------------------------------------------------------------------
# Build B(t) = 1 + b_1 t + ... + b_N t^N
# ------------------------------------------------------------------

B = [Fraction(b[i]) for i in range(N+1)]

# Verify B(t) = 1/(1-A(t)) via a_k = INVERTi(b_k)
A = [Fraction(0)] * (N+1)
for k in range(1, N+1):
    A[k] = Fraction(a[k])

# check B * (1 - A) == 1
one_minus_A = [Fraction(1)] + [-A[k] for k in range(1, N+1)]
prod_check = ps_mul(B, one_minus_A, N)
assert prod_check[0] == 1
for k in range(1, N+1):
    assert prod_check[k] == 0, f"INVERTi consistency fails at k={k}: got {prod_check[k]}"
print("INVERTi consistency B*(1-A) = 1: PASS")

# ------------------------------------------------------------------
# Compute sum_n p_n t^n = sum_{d>=1} (mu(d)/d) log B(t^d)
# ------------------------------------------------------------------

P_series = [Fraction(0)] * (N+1)
for d in range(1, N+1):
    mu = mobius(d)
    if mu == 0:
        continue
    Btd = ps_substitute_td(B, d, N)
    logBtd = ps_log(Btd, N)
    P_series = ps_add(P_series, ps_scalar(logBtd, Fraction(mu, d), N), N)

# Check all p_n are integers
p = [None] * (N+1)
for k in range(1, N+1):
    val = P_series[k]
    assert val.denominator == 1, f"p_{k} not integer: {val}"
    p[k] = int(val)

print()
print("Witt-formula primitives p_k = dim L(a)_k:")
for k in range(1, N+1):
    print(f"  p_{k:2d} = {p[k]}")

print()
print("Mod 3 pattern:")
print("  " + ", ".join(str(p[k] % 3) for k in range(1, N+1)))

print()
print("Mod 9 pattern:")
print("  " + ", ".join(str(p[k] % 9) for k in range(1, N+1)))

# ------------------------------------------------------------------
# SANITY CHECK: PBW ==> B(t) = prod_{n>=1} (1 - t^n)^{-p_n}
# ------------------------------------------------------------------

# Build prod_{n=1}^{N} (1 - t^n)^{-p_n} mod t^{N+1}
# For each n, (1 - t^n)^{-p_n} = sum_{j>=0} C(p_n + j - 1, j) t^{nj}.

def binom(a_val, k):
    # a_val is a plain int (possibly large), k >= 0
    if k < 0:
        return 0
    r = 1
    for i in range(k):
        r = r * (a_val - i) // (i + 1)
    return r

pbw = [Fraction(0)] * (N+1)
pbw[0] = Fraction(1)

running = [Fraction(1)] + [Fraction(0)] * N
for n in range(1, N+1):
    pn = p[n]
    if pn == 0:
        continue
    # factor for exponent pn:
    factor_poly = [Fraction(0)] * (N+1)
    jmax = N // n
    for j in range(0, jmax+1):
        factor_poly[n*j] = Fraction(binom(pn + j - 1, j))
    running = ps_mul(running, factor_poly, N)

print()
print("PBW sanity check: does prod (1-t^n)^(-p_n) equal B(t)?")
ok = True
for k in range(N+1):
    lhs = running[k]
    rhs = B[k]
    if lhs != rhs:
        ok = False
        print(f"  MISMATCH at t^{k}: prod = {lhs}, B = {rhs}")
print("  Sanity: " + ("PASS" if ok else "FAIL"))

# ------------------------------------------------------------------
# OEIS lookup
# ------------------------------------------------------------------

print()
print("OEIS lookup:")
p_vals = [p[k] for k in range(1, N+1)]
query = ",".join(str(x) for x in p_vals)
print(f"  Query: {query}")

try:
    import urllib.request
    url = f"https://oeis.org/search?q={query}&fmt=text"
    req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
    with urllib.request.urlopen(req, timeout=30) as resp:
        text = resp.read().decode('utf-8', errors='replace')
    # Look for A-numbers in results
    import re
    matches = re.findall(r'^%I (A\d+)', text, re.MULTILINE)
    if not matches:
        # try alternative pattern
        matches = re.findall(r'A\d{6}', text)
    if matches:
        seen = set()
        uniq = [m for m in matches if not (m in seen or seen.add(m))]
        print(f"  OEIS hits: {uniq[:10]}")
        # print first few lines of response
        lines = text.split('\n')
        for line in lines[:60]:
            if line.startswith('%') or 'No results' in line or 'showing' in line.lower():
                print(f"    {line}")
    else:
        # inspect for "no results"
        if 'no results' in text.lower() or 'not in' in text.lower() or 'showing 0' in text.lower():
            print("  OEIS reports NO MATCH — sequence appears NEW.")
        else:
            print("  No A-numbers extracted; raw excerpt:")
            for line in text.split('\n')[:20]:
                print(f"    {line}")
except Exception as e:
    print(f"  OEIS lookup FAILED: {e}")
    print(f"  Manual check: paste {query} into https://oeis.org")
