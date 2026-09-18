"""Sober re-check: recompute p_3(Y).e_3 at m=7 (log used m=6=r+3, minimum)."""
import sys
sys.path.insert(0, '/home/agent/projects/proofs/scripts/day201')
from p3Y_er import compute_p3Y_er
import sympy as sp
import time

q, t = sp.symbols('q t')

r = 3
m = 7
print(f'Recompute p_3(Y).e_{r} at m={m}, target degree n={r+3}...', flush=True)
t0 = time.time()
exp = compute_p3Y_er(r, m)
print(f'  done in {time.time()-t0:.1f}s', flush=True)

top = exp.get((r+3,), sp.Integer(0))
top_simp = sp.factor(sp.simplify(top))
print(f'  top coeff at ({r+3},) at m=7 =')
print(f'    {top_simp}')

# Compare with log (m=6)
tau3_log = ((q - 1)*(t + 1)*(q**2 + q + 1)*(t**2 - t + 1) *
            (q**3*t**9 - q**3*t**5 - q**3*t**4 + q**3 + q**2*t**6 + q**2*t**5
             + q**2*t**4 - q**2*t**2 - q**2*t - q**2 + q*t**3 - q + t**2 + t + 1)) / q**6
diff = sp.simplify(top - tau3_log)
print(f'  diff (m=7 minus m=6 from log): {diff}')
