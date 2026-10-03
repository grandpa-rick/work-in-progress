"""Day 220: at t=0, e*_lam = omega Htilde_lam(x;s) (Day 217e Thm B), so c_{lam mu}(s,0) = [h_mu] s^{n(lam)} prod_{i<j}(1-R_ij)/(1-R_ij/s) h_lam,
(1-R)/(1-qR) = 1 - (1-q) sum_{m>=1} q^{m-1} R^m, q=1/s. eps = s-1, truncated series; E = number of edges.
Pairs processed column by column (j descending): entry j is raised only by pairs (j,k), k>j (earlier columns), and lowered
only in column j, so entries never need to go negative (prune). Compare v and lead with val_n*_t0.log (subset-formula engine);
check v = E_min = l - kappa and lead = (-1)^E_min * #(min configs)."""
import re
from fractions import Fraction as Fr
from kappa import kappa
def run(lam):
    l = len(lam); n = sum(lam); P = l + 2
    def mul(a, b):
        r = [Fr(0)]*P
        for i, x in enumerate(a):
            if x:
                for j, y in enumerate(b[:P-i]): r[i+j] += x*y
        return r
    qser = [Fr((-1)**k) for k in range(P)]
    qp = [[Fr(1)]+[Fr(0)]*(P-1)]
    for _ in range(n): qp.append(mul(qp[-1], qser))
    omq = [Fr(0)] + [-x for x in qser[1:]]
    spow = [Fr(1)]+[Fr(0)]*(P-1)
    for _ in range(sum(i*x for i, x in enumerate(lam))): spow = mul(spow, [Fr(1), Fr(1)]+[Fr(0)]*(P-2))
    states = {(tuple(lam), 0): spow}
    order = [(i, j) for j in range(l-1, 0, -1) for i in range(j)]
    for (i, j) in order:
        new = {}
        def add(key, w):
            new[key] = [a+b for a, b in zip(new.get(key, [Fr(0)]*P), w)]
        for (vec, E), w in states.items():
            add((vec, E), w)
            for m in range(1, vec[j]+1):
                v2 = list(vec); v2[i] += m; v2[j] -= m
                ww = [-x for x in mul(mul(w, omq), qp[m-1])]
                add((tuple(v2), E+1), ww)
        states = new
    res = {}
    for (vec, E), w in states.items():
        mu = tuple(sorted([x for x in vec if x > 0], reverse=True))
        d = res.setdefault(mu, {}); d[E] = [a+b for a, b in zip(d.get(E, [Fr(0)]*P), w)]
    return res
if __name__ == '__main__':
    tot = ok = 0
    for f in ['val_n4_t0.log', 'val_n5_t0.log', 'val_n6_t0.log']:
        cache = {}
        for line in open(f):
            m = re.match(r'\s*(\([\d, ]*\)) -> (\([\d, ]*\)): v=(\d+) .* lead=(\S+)', line)
            if not m: continue
            lam, mu, v, ld = eval(m.group(1)), eval(m.group(2)), int(m.group(3)), Fr(m.group(4))
            if lam not in cache: cache[lam] = run(lam)
            byE = cache[lam].get(mu, {})
            ser = [sum(byE[E][k] for E in byE) for k in range(len(lam)+2)]
            vv = next(k for k, x in enumerate(ser) if x != 0); lead = ser[vv]
            Emin = min(E for E in byE if any(byE[E])); Nmin = byE[Emin][Emin]*(-1)**Emin
            good = (vv == v and lead == ld and Emin == len(lam)-kappa(lam, mu) == vv and lead == (-1)**Emin*Nmin)
            tot += 1; ok += good
            if not good: print('BAD', lam, mu, 'log v', v, 'raising v', vv, 'log lead', ld, 'raising lead', lead, 'Emin', Emin, 'Nmin', Nmin)
    print('raising-operator t=0 check:', ok, '/', tot)
