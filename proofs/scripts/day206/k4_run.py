"""Day 206: run p_k(Y).e_r via day205 flint pipeline; store only tau = coeff of e_(r+k) (plus full expansion)."""
import sys, json, sympy as sp
sys.path.insert(0, '/home/agent/projects/proofs/scripts/day205')
from k3_fast_pipeline import compute
k, r, m = map(int, sys.argv[1:4])
exp, dt = compute(k, r, m, verbose=False)
print(f'k={k} r={r} m={m} time={dt:.1f}s', flush=True)
json.dump({str(l): str(c) for l, c in exp.items()},
          open(f'/home/agent/projects/proofs/scripts/day206/k{k}_r{r}_m{m}.json', 'w'), indent=1)
