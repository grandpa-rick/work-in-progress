"""Day 208 (Sep 29): factored e-basis tables for t^{-C(k,2)} e_k(Y).(e_a e_b), with m-stability check (m=n, n+1)."""
from ek_two_col import compute, sp
import sys
cases = {2: [(1,1),(2,1),(2,2),(3,1),(3,2),(3,3)], 3: [(1,1),(2,1),(2,2)]}
for k, L in cases.items():
    for (a, b) in L:
        n = a + b + k
        r1, dt1 = compute(k, a, b, n); r2, dt2 = compute(k, a, b, n + 1)
        stable = all(sp.expand(r1.get(l, 0) - r2.get(l, 0)) == 0 for l in set(r1) | set(r2))
        print(f'k={k} (a,b)=({a},{b}) m={n},{n+1} stable={stable} maxlen={max(len(l) for l in r1)} [{dt1+dt2:.1f}s]')
        for l in sorted(r1, reverse=True):
            print('    e', l, ':', sp.factor(r1[l]))
        sys.stdout.flush()
