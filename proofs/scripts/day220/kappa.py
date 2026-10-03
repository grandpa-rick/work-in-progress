"""kappa(lam,mu) = max #blocks r of a set partition of the parts of lam into blocks C, with mu = union of nu^C, nu^C a partition of |lam_C|, nu^C dominates lam_C.
Predicted (Day 220): v_{(s-1)} c_{lam mu} = l(lam) - kappa(lam,mu) for generic t (lower bound PROVED)."""
import sys, re, functools
from itertools import combinations
def dom(a, b):
    sa = sb = 0
    for i in range(max(len(a), len(b))):
        sa += a[i] if i < len(a) else 0; sb += b[i] if i < len(b) else 0
        if sa < sb: return False
    return True
def subms(m):  # sub-multisets of a sorted tuple (as sorted tuples), nonempty
    out = set()
    for r in range(1, len(m)+1):
        for c in combinations(range(len(m)), r): out.add(tuple(m[i] for i in c))
    return out
def minus(m, s):
    m = list(m)
    for x in s: m.remove(x)
    return tuple(m)
@functools.lru_cache(None)
def kappa(lam, mu):  # lam, mu sorted desc tuples, same size; -inf if impossible
    if not lam and not mu: return 0
    if not lam or not mu: return -10**9
    best = -10**9; a = lam[0]  # block containing the first part of lam
    for L in subms(lam[1:]) | {()}:
        C = (a,)+L; Cs = tuple(sorted(C, reverse=True)); n = sum(C)
        for M in subms(mu):
            if sum(M) == n and dom(tuple(sorted(M, reverse=True)), Cs):
                best = max(best, 1 + kappa(minus(lam, C), minus(mu, M)))
    return best
if __name__ == '__main__':
    tot = ok = 0
    for f in sys.argv[1:]:
        for line in open(f):
            m = re.match(r'\s*(\([\d, ]*\)) (?:->|) ?(\([\d, ]*\)): v=(\d+)', line) or re.match(r'(\([\d, ]*\)) (\([\d, ]*\)) v= (\d+)', line)
            if not m: continue
            lam, mu, v = eval(m.group(1)), eval(m.group(2)), int(m.group(3))
            k = kappa(tuple(sorted(lam, reverse=True)), tuple(sorted(mu, reverse=True)))
            tot += 1; ok += (v == len(lam)-k)
            if v != len(lam)-k: print(f, lam, mu, 'v=', v, 'l-kappa=', len(lam)-k)
        print(f, 'running total', ok, '/', tot)
