"""Day 220: check Lead_{lam mu}(t) = sum over tight merge histories of prod W_k(J), W from closed form
W_k(J) = (-1)^p [n]_t prod_i [k]_{t^{j_i}} / [k]_t  (computed n<=6 in merge_weights_n6.log), against logged leads (mu a coarsening of lam)."""
import sys, re, functools
from fractions import Fraction as Fr
from itertools import combinations
from kappa import kappa
def qi(m, x): return sum(x**i for i in range(m))
def W(k, J, t):
    n = k+sum(J); r = Fr((-1)**len(J))*qi(n, t)/qi(k, t)
    for j in J: r *= qi(k, t**j)
    return r
def lead(lam, mu, t):
    # process lam[l-1], ..., lam[0]; state = sorted tuple of block sums (multiset); weight sums
    states = {(): Fr(1)}
    for k in reversed(lam):
        new = {}
        for st, w in states.items():
            idx = range(len(st))
            for p in range(0, len(st)+1):
                for S in combinations(idx, p):
                    J = [st[i] for i in S]; rest = [st[i] for i in idx if i not in S]
                    ww = w * (W(k, J, t) if p else 1)
                    key = tuple(sorted(rest+[k+sum(J)], reverse=True)); new[key] = new.get(key, 0) + ww
        states = new
    return states.get(tuple(mu), 0)
tot = ok = 0
for f, t in [('val_n6_t3_5.log', Fr(3, 5)), ('val_n6_t7_3.log', Fr(7, 3)), ('val_n6_t-2.log', Fr(-2))]:
    for line in open(f):
        m = re.match(r'\s*(\([\d, ]*\)) -> (\([\d, ]*\)): v=(\d+) .* lead=(\S+)', line)
        if not m: continue
        lam, mu, v, ld = eval(m.group(1)), eval(m.group(2)), int(m.group(3)), Fr(m.group(4))
        if kappa(lam, mu) != len(mu) or len(mu) == len(lam): continue  # coarsenings only
        L = lead(lam, mu, t); tot += 1
        good = (v == len(lam)-len(mu) and L == ld) or (L == 0 and v > len(lam)-len(mu))
        ok += good
        print(f, lam, mu, 'v', v, 'logged lead', ld, 'history lead', L, 'OK' if good else 'MISMATCH')
print('history-lead check', ok, '/', tot)
print('lead at t=0 (= (-1)^m #histories): (1^6)->(6):', lead((1,)*6, (6,), Fr(0)), '; (1,1,1)->(3):', lead((1,1,1),(3,),Fr(0)))
