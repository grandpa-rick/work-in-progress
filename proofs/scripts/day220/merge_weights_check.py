"""Check computed W_k(J) against closed form (-1)^p [n]_t prod [k]_{t^j}/[k]_t for n<=MAXN."""
import sys
from sympy import symbols, cancel
from merge_weights import W, parts
t = symbols('t')
def q(m, x): return sum(x**i for i in range(m))
MAXN = int(sys.argv[1]); tot = ok = 0
for n in range(2, MAXN+1):
    for k in range(1, n):
        for J in parts(n-k):
            cf = (-1)**len(J)*q(n, t)/q(k, t)
            for j in J: cf *= q(k, t**j)
            good = cancel(W(k, J) - cf) == 0; tot += 1; ok += good
            print(n, k, J, 'OK' if good else 'MISMATCH', flush=True)
    print(f'through n={n}: {ok}/{tot}', flush=True)
