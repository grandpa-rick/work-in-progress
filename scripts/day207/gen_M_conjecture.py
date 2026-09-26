"""Test conjectural closed form for M_{jb}^{(k)} against engine output, k=2,3,4 (and 5 if pickle exists)."""
import sympy as sp, pickle, os
s, t = sp.symbols('s t')
def tp(n): return sp.Mul(*[1 - t**i for i in range(1, n+1)])
def sp_(n): return sp.Mul(*[1 - s*t**i for i in range(n)])   # (s;t)_n
def Mconj(k, j, b):
    R = s*(t**(k-b) - t**(k-b-j) + 1) - t**j
    return sp.Mul(*[s - t**i for i in range(j)]) * s**b * R * sp_(k-j-b) / (1 - s) \
        * tp(k) / (tp(k-j-b) * tp(j)) / t**(k*j)
for k in (2, 3, 4, 5):
    f = f'/home/agent/projects/scripts/day207/k{k}_gf_coeffs.pkl'
    if not os.path.exists(f): continue
    C, M, coeffs = pickle.load(open(f, 'rb'))
    ratios = {key: sp.factor(sp.cancel(M[key] / Mconj(k, *key))) for key in M}
    miss = [(j, b) for j in range(k+1) for b in range(k+1-j) if (j, b) not in M]
    print(k, 'support = {j+b<=k}:', not miss and all(j+b <= k for j, b in M), ' ratios:', sorted(set(ratios.values()), key=str), {kk: v for kk, v in ratios.items()})

# Implied Pieri conjecture: coef of e_{r+n} e_b in t^{-C(k,2)} e_k(Y).e_r (k = b+n) equals s^b F_n(t^{r-b}),
# F_n(w) = 1/(s-1) sum_{j=0}^n (s;t)_{n-j} prod_{i<j}(s-t^i) (s(t^n - t^{n-j} + 1) - t^j) w^j / ((t;t)_{n-j} (t;t)_j)
u = sp.Symbol('u')
def F(n, w):
    return sum(sp_(n-j) * sp.Mul(*[s - t**i for i in range(j)]) * (s*(t**n - t**(n-j) + 1) - t**j) * w**j / (tp(n-j) * tp(j)) for j in range(n+1)) / (s - 1)
for k in (2, 3, 4, 5):
    f = f'/home/agent/projects/scripts/day207/k{k}_gf_coeffs.pkl'
    if not os.path.exists(f): continue
    _, _, coeffs = pickle.load(open(f, 'rb'))
    good = all(sp.cancel(coeffs[(n, k-n)] - s**(k-n) * F(n, u * t**(-(k-n)))) == 0 for n in range(k+1)) and len(coeffs) == k+1
    print(k, 'coef(e_{r+n}e_{k-n}) == s^{k-n} F_n(t^{r-k+n}) for all n:', good)
