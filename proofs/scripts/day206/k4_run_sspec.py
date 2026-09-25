"""Day 206: same pipeline as day205/k3_fast_pipeline.py but with s=q^{-1} specialised to an integer s0
(pi F = X_1 F(X_2..X_m, s0 X_1)); t stays symbolic. Memory-saving variant for r=9 (full (s,t) run OOM'd at 23GB)."""
import sys, json, time, sympy as sp
sys.path.insert(0, '/home/agent/projects/proofs/scripts/day205')
import k3_fast_pipeline as K
k, r, m, s0 = map(int, sys.argv[1:5])
A = K.AHA(m)
A._pi[m-1] = s0 * A.X[0]
t0 = time.time(); er = A.e(r); tot = 0*A.s
for i in range(1, m+1):
    v = er
    for _ in range(k): v = A.Y(v, i)
    tot = tot + v
    print(f'  i={i} {time.time()-t0:.1f}s terms={len(tot)}', flush=True)
exp = K.e_expansion(A, tot, r+k)
tau = K.to_sympy(A, exp[(r+k,)])  # no s left; t symbolic
print(f'k={k} r={r} m={m} s0={s0} time={time.time()-t0:.1f}s', flush=True)
json.dump({'tau': str(tau), 's0': s0}, open(f'/home/agent/projects/proofs/scripts/day206/k{k}_r{r}_m{m}_s{s0}.json','w'))
