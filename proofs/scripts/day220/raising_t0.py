"""Day 220: at t=0, e*_lam = omega Htilde_lam(x;s) (Day 217e Thm B), so c_{lam mu}(s,0) = [h_mu] s^{n(lam)} prod_{i<j}(1-R_ij)/(1-R_ij/s) h_lam.
Expand (1-R)/(1-qR) = 1 - (1-q) sum_{m>=1} q^{m-1} R^m.  Compute c(s,0) exactly; compare v and lead with val_n*_t0.log;
also check v = E_min = l - kappa and lead = (-1)^{E_min} #(minimal configs)."""
import re, itertools
from sympy import symbols, expand, Poly, cancel, Rational
from kappa import kappa
s = symbols('s'); q = 1/s
def coeffs(lam):
    l = len(lam); pairs = [(i, j) for i in range(l) for j in range(i+1, l)]; n = sum(lam)
    res = {}; cnt = {}
    def rec(idx, vec, w, E):
        if idx == len(pairs):
            if min(vec) < 0: return
            mu = tuple(sorted([x for x in vec if x > 0], reverse=True))
            res[mu] = res.get(mu, 0) + w
            cnt.setdefault(mu, {}); cnt[mu][E] = cnt[mu].get(E, 0) + 1
            return
        i, j = pairs[idx]
        rec(idx+1, vec, w, E)
        for m in range(1, n+1):
            v2 = list(vec); v2[i] += m; v2[j] -= m
            if v2[j] < -n: break
            rec(idx+1, v2, w*(-(1-q))*q**(m-1), E+1)
    rec(0, list(lam), s**sum(i*x for i, x in enumerate(lam)), 0)
    return {mu: cancel(c) for mu, c in res.items()}, cnt
tot = ok = 0
for f in ['val_n4_t0.log', 'val_n5_t0.log', 'val_n6_t0.log']:
    cache = {}
    for line in open(f):
        m = re.match(r'\s*(\([\d, ]*\)) -> (\([\d, ]*\)): v=(\d+) .* lead=(\S+)', line)
        if not m: continue
        lam, mu, v, ld = eval(m.group(1)), eval(m.group(2)), int(m.group(3)), Rational(m.group(4))
        if lam not in cache: cache[lam] = coeffs(lam)
        C, cnt = cache[lam]; c = C.get(mu, 0)
        pc = Poly(expand(c*s**60), s); vv = 0
        while pc.eval(1) == 0: pc = Poly(cancel(pc.as_expr()/(s-1)), s); vv += 1
        lead = pc.eval(1)  # s^60 factor is 1 at s=1
        Emin = min(e for e, k in cnt[mu].items()); Nmin = cnt[mu][Emin]
        good = (vv == v and lead == ld and Emin == len(lam)-kappa(lam, mu) and lead == (-1)**Emin*Nmin)
        tot += 1; ok += good
        if not good: print('BAD', lam, mu, v, vv, ld, lead, Emin, Nmin)
print('raising-operator t=0 check:', ok, '/', tot)
