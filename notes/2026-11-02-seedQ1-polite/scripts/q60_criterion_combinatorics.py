#!/usr/bin/env python3
"""Audit 2026-11-02 (Q60). Pure root-system combinatorics, no crystal computation.
For A_n (n=2,3): enumerate Weyl chambers of theta in Q_R^* via S = {alpha>0 : theta(alpha)>0}
(biconvex sets), and for each j classify (theta, j):
  LOC  := no positive root a with theta(a)<=0<theta(a+alpha_j)            (Q60 criterion)
  ALT  := alpha_j notin S  or  T u {alpha_j} biconvex   (T = complement of S)
  CASE := 'neg'  (alpha_j in T: covered by TW Prop 3.13 / Kamnitzer Thm 3.5 deduction)
          'adj'  (alpha_j in S, T u {alpha_j} biconvex: LOC true, NOT covered by literature)
          'far'  (alpha_j in S, not adj: LOC false; non-locality claim, NOT covered)
Checks LOC == ALT on every chamber. Grade: computed."""
import itertools, random
def roots(n):
    return [tuple(1 if i <= k <= j else 0 for k in range(n)) for i in range(n) for j in range(i, n)]
def add(a, b): return tuple(x + y for x, y in zip(a, b))
def closed(X, R):
    return all(add(a, b) not in R or add(a, b) in X for a in X for b in X)
for n in (2, 3):
    R = roots(n); simple = [tuple(int(k == j) for k in range(n)) for j in range(n)]
    chambers = {}
    random.seed(1)
    for _ in range(20000):
        th = [random.gauss(0, 1) for _ in range(n)]
        S = frozenset(a for a in R if sum(t * x for t, x in zip(th, a)) > 0)
        chambers.setdefault(S, th)
    cnt = {'neg': 0, 'adj': 0, 'far': 0}; agree = 0; tot = 0
    for S, th in chambers.items():
        T = frozenset(R) - S
        assert closed(S, R) and closed(T, R)
        for j, aj in enumerate(simple):
            loc = not any(a not in S and add(a, aj) in S for a in R if add(a, aj) in R)
            Tj = T | {aj}
            alt = (aj not in S) or (closed(Tj, R) and closed(frozenset(R) - Tj, R))
            case = 'neg' if aj not in S else ('adj' if alt else 'far')
            cnt[case] += 1; tot += 1; agree += (loc == alt)
            assert (case == 'far') == (not loc)
    print(f"A{n}: chambers found = {len(chambers)}; (theta,j) pairs = {tot}; LOC==ALT on {agree}/{tot}; "
          f"cases: neg={cnt['neg']} adj={cnt['adj']} far={cnt['far']}")
