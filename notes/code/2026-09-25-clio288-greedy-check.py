"""Independent check of Clio note 1 (cylindric M-convexity), Section 5.
Bead model: x=(x_1..x_m), x_1<...<x_m<x_1+n, x_{i+m}=x_i+n.
Horizontal strip nu<<rho: nu_i <= rho_i < nu_{i+1} (i=1..m, nu_{m+1}=nu_1+n).
Chains enumerated from the RAW definition (1), not via Lemma 3.1."""
import itertools, sys
from collections import Counter

def ext(x, n, i):  # x_i for any integer i (1-indexed)
    m = len(x); q, r = divmod(i-1, m)
    return x[r] + q*n

def strips_from(x, lam, n):
    m = len(x)
    ranges = [range(x[i], min(ext(x, n, i+2)-1, lam[i])+1) for i in range(m)]
    for y in itertools.product(*ranges):
        if all(y[i] < y[i+1] for i in range(m-1)) and y[-1] < y[0]+n:
            yield y

def is_strip(nu, rho, n):
    m = len(nu)
    return all(nu[i] <= rho[i] < ext(nu, n, i+2) for i in range(m))

def weights(mu, lam, n, ell):
    # dict: weight vector -> multiplicity (number of chains)
    cur = {mu: Counter({(): 1})}
    for t in range(ell):
        nxt = {}
        for x, ws in cur.items():
            for y in strips_from(x, lam, n):
                assert is_strip(x, y, n)
                a = sum(y) - sum(x)
                c = nxt.setdefault(y, Counter())
                for w, k in ws.items():
                    c[w + (a,)] += k
        cur = nxt
    return cur.get(lam, Counter())

def greedy(mu, lam, n, T):
    g = [mu]
    for t in range(T):
        p = g[-1]
        g.append(tuple(min(ext(p, n, i+2)-1, lam[i]) for i in range(len(mu))))
    return g

def dom(a, b):  # a <= b in dominance (same size), padded
    L = max(len(a), len(b)); a = list(a)+[0]*(L-len(a)); b = list(b)+[0]*(L-len(b))
    sa = sb = 0
    for x, y in zip(a, b):
        sa += x; sb += y
        if sa > sb: return False
    return True

def srt(a): return tuple(sorted((v for v in a if v > 0), reverse=True))

def comps(d, ell):
    if ell == 1: yield (d,); return
    for a in range(d+1):
        for r in comps(d-a, ell-1): yield (a,)+r

def shapes(n, m, lo, hi):
    for x in itertools.combinations(range(lo, hi), m):
        if x[-1] < x[0]+n: yield x

stats = Counter(); fails = []
NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 6
DMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 7
LMAX = int(sys.argv[3]) if len(sys.argv) > 3 else 5
for n in range(2, NMAX+1):
    for m in range(1, n):
        for mu in shapes(n, m, 0, n):
            if mu[0] != 0: continue
            for lam in shapes(n, m, 0, n+DMAX+1):
                if not all(l >= u for l, u in zip(lam, mu)): continue
                d = sum(lam)-sum(mu)
                if d > DMAX: continue
                g = greedy(mu, lam, n, 200)
                l0 = next((t for t, x in enumerate(g) if x == lam), None)
                if l0 is None: fails.append(('l0 infinite', n, mu, lam)); continue
                stats['shapes'] += 1
                for ell in range(1, LMAX+1):
                    W = weights(mu, lam, n, ell)
                    if (len(W) > 0) != (ell >= l0):
                        fails.append(('5.4 nonempty iff ell>=l0', n, mu, lam, ell)); continue
                    if not W: continue
                    stats['instances'] += 1
                    gam = tuple(sum(g[t])-sum(g[t-1]) for t in range(1, ell+1))
                    lh = srt(gam)
                    # 5.2(1): steps are strips, g^t subset lam
                    if not all(is_strip(g[t-1], g[t], n) and all(a <= b for a, b in zip(g[t], lam)) for t in range(1, ell+1)):
                        fails.append(('5.2(1)', n, mu, lam, ell))
                    if gam not in W: fails.append(('5.2(2) gamma not in W', n, mu, lam, ell))
                    if not all(dom(srt(w), lh) for w in W): fails.append(('5.5 max', n, mu, lam, ell))
                    # prefix-sum domination in position order (the actual content of 5.5)
                    for w in W:
                        if any(sum(w[:r]) > sum(gam[:r]) for r in range(1, ell+1)):
                            fails.append(('prefix', n, mu, lam, ell, w)); break
                    # free corollary: gamma already weakly decreasing
                    if any(gam[t] < gam[t+1] for t in range(ell-1)):
                        fails.append(('gamma not sorted', n, mu, lam, ell, gam))
                    else: stats['gamma_sorted'] += 1
                    # Thm 7.1: W = {alpha : sort alpha <= lh}
                    target = {a for a in comps(d, ell) if dom(srt(a), lh)}
                    if set(W) != target: fails.append(('7.1 support', n, mu, lam, ell))
                    # symmetry of coefficients (Cor 3.4 strengthened)
                    if any(W[w] != W[tuple(sorted(w, reverse=True))] for w in W):
                        fails.append(('coef symmetry', n, mu, lam, ell))
                    # 5.4: lambda-hat independent of ell
                    if ell >= l0:
                        stats['ell_indep_checked'] += 1
                        base = srt(tuple(sum(g[t])-sum(g[t-1]) for t in range(1, l0+1)))
                        if base != lh: fails.append(('5.4 ell-indep', n, mu, lam, ell))
print(dict(stats)); print('failures:', len(fails)); print(fails[:10])
