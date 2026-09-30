# engine convention check: E_k e_r (subset formula) vs 207b closed form sum_b s^b F_{k-b}(t^{r-b}) e_b e_{r+k-b}, symbolic s,t
import sympy as sp
from ek_subset_engine import *
s, t = sp.symbols('s t')
def tt(n): return sp.prod([1-t**i for i in range(1, n+1)])
def sp_(a, n): return sp.prod([1-a*t**i for i in range(n)])
def alpha(j): return 0 if j < 0 else sp.prod([s-t**i for i in range(1, j+1)])/tt(j)
def c(n, j): return 0 if (j < 0 or j > n) else sp_(s, n-j)/tt(n-j)*(alpha(j)-s*t**(n-j)*alpha(j-1))
def Fn(n, w): return sum(c(n, j)*w**j for j in range(n+1))
ok = True
for k in range(1, 4):
    for r in range(0, 4):
        m = k + r; xs = setup(m)
        lhs = Ek(e(r, xs), k, xs, s, t)
        rhs = sum(s**b*Fn(k-b, t**(r-b))*e(b, xs)*e(r+k-b, xs) for b in range(0, k+1))
        d = sp.simplify(sp.together(sp.expand(lhs) - sp.expand(rhs)))
        print(k, r, d == 0, flush=True); ok &= d == 0
print('ALL OK' if ok else 'FAIL')
