# Day 235: printed Thm 3.4 says: at t=0 the lead [(s-1)^{l-kappa}] c_{lam,mu} is (-1)^{l-kappa} N(lam,mu), N>=1
# "a count of minimal raising-operator configurations". Proposed one-phrase definition (from the printed proof idea,
# Macdonald III (2.15') raising operators R_ij, i<j):
#   N(lam,mu) = #{ (k_ij)_{1<=i<j<=l} in Z_{>=0} : exactly l-kappa entries nonzero,
#                  lam + sum k_ij (eps_i - eps_j) is a rearrangement of mu (padded with zeros) }.
# Ground truth c_{lam,mu}(s,0): scripts/day234 subset-formula table (n<=7), used as a library only.
# NEGATIVE CONTROL (breaks the definition, not a relabelling): count arrays with exactly l-kappa+1 nonzero entries
# (the next order); the claimed identity must then FAIL.
import sys, os, itertools
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'day234'))
from common import *   # load_table, kappa, dom, coeff_sm1, partitions, R, s, t

from functools import lru_cache
def count(lam, mu, edges):
    l = len(lam)
    if len(mu) > l: return 0
    pairs = tuple((i, j) for i in range(l) for j in range(i + 1, l))
    target = tuple(sorted(list(mu) + [0] * (l - len(mu))))
    n = sum(lam)
    @lru_cache(None)
    def rec(idx, vec, rem):
        if rem < 0 or len(pairs) - idx < rem: return 0
        if idx == len(pairs):
            return 1 if tuple(sorted(vec)) == target else 0
        i, j = pairs[idx]
        tot = rec(idx + 1, vec, rem)          # k_ij = 0
        if rem > 0:
            for k in range(1, n + 1):
                v = list(vec); v[i] += k; v[j] -= k
                if v[i] > n or v[j] < -n: break
                tot += rec(idx + 1, tuple(v), rem - 1)
        return tot
    return rec(0, tuple(lam), edges)

T = load_table()
ok = bad = nc_fire = nc_tot = 0
for n in range(2, 8):
    for lam in partitions(n):
        for mu in partitions(n):
            if mu == lam or not dom(mu, lam): continue
            l = len(lam); kap = kappa(lam, mu); m = l - kap
            c0 = c(T, lam, mu).compose(t, R(0))
            lead = coeff_sm1(c0, m)                     # element of R with no s,t left
            L = int(lead.LC) if lead != 0 else 0
            N = count(lam, mu, m)
            good = (L == (-1) ** m * N) and N >= 1
            ok += good; bad += not good
            if not good: print('FAIL', lam, mu, 'lead', L, 'N', N)
            Nc = count(lam, mu, m + 1); nc_tot += 1
            nc_fire += (L != (-1) ** m * Nc)
print(f"Thm 3.4 t=0 lead = (-1)^(l-kappa) N, N as defined, n<=7: {ok} pass, {bad} fail")
print(f"negative control (count l-kappa+1 edges instead): identity fails on {nc_fire}/{nc_tot}")
