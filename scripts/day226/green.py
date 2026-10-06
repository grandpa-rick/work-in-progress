# Green polynomials X^lam_rho(t) via Kostka-Foulkes (charge) and MN characters. grade: computed
# Convention: Macdonald III (7.1): p_rho = sum_lam X^lam_rho(t) P_lam(x;t);  X^lam_rho = sum_mu chi^mu_rho K_{mu lam}(t).
import sympy as sp, functools, itertools
t = sp.symbols('t')
def parts(n, m=None):
    if m is None: m = n
    if n == 0: yield (); return
    for k in range(min(n, m), 0, -1):
        for r in parts(n-k, k): yield (k,)+r
def nfun(l): return sum(i*x for i, x in enumerate(l))
def mult(l):
    d = {}
    for x in l: d[x] = d.get(x, 0)+1
    return d
def phi(m): return sp.prod([1-t**i for i in range(1, m+1)])
def b(l): return sp.prod([phi(m) for m in mult(l).values()])
def zee(l):
    from math import factorial
    r = 1
    for k, m in mult(l).items(): r *= k**m*factorial(m)
    return r
# --- MN characters via beta numbers
@functools.lru_cache(None)
def chi(lam, rho):
    if not rho: return 1 if not lam else 0
    r, rest = rho[0], rho[1:]
    L = len(lam); beta = [lam[i]+L-1-i for i in range(L)]
    tot = 0; S = set(beta)
    for bb in beta:
        if bb-r >= 0 and bb-r not in S:
            sign = (-1)**sum(1 for c in beta if bb-r < c < bb)
            nb = sorted([c for c in beta if c != bb]+[bb-r], reverse=True)
            L2 = len(nb); mu = tuple(x for x in (nb[i]-(L2-1-i) for i in range(L2)) if x > 0)
            tot += sign*chi(mu, rest)
    return tot
# --- SSYT of shape lam, content mu
def ssyt(lam, mu):
    rows = [[] for _ in lam]
    res = []
    def place(letter, remaining, shape_fill):
        # fill `remaining` copies of `letter` as a horizontal strip
        if letter > len(mu): res.append([r[:] for r in rows]); return
        if remaining == 0:
            nxt = letter+1
            place(nxt, mu[nxt-1] if nxt <= len(mu) else 0, None); return
        # choose counts per row: horizontal strip
        cur = [len(r) for r in rows]
        def rec(i, left, added):
            if i == len(lam):
                if left == 0:
                    for j, c in enumerate(added): rows[j].extend([letter]*c)
                    nxt = letter+1
                    place(nxt, mu[nxt-1] if nxt <= len(mu) else 0, None)
                    for j, c in enumerate(added):
                        for _ in range(c): rows[j].pop()
                return
            cap = lam[i]-cur[i]
            if i > 0: cap = min(cap, cur[i-1]-cur[i])  # strip: new cells in row i must lie under old cells of row i-1
            for c in range(min(cap, left), -1, -1): rec(i+1, left-c, added+[c])
        rec(0, remaining, [])
    place(1, mu[0], None)
    return res
def charge_word(w):
    w = list(w); n = max(w) if w else 0; tot = 0
    used = [False]*len(w)
    while not all(used):
        # extract standard subword: start with rightmost unused 1, then search leftwards cyclically for 2, ...
        pos = None; idx = 0; letters = []
        k = 1
        # find rightmost unused 1
        cand = [i for i in range(len(w)) if not used[i] and w[i] == 1]
        if not cand: break
        pos = cand[-1]; used[pos] = True; letters.append(pos)
        k = 2
        while True:
            c = [i for i in range(len(w)) if not used[i] and w[i] == k]
            if not c: break
            left = [i for i in c if i < pos]
            if left: pos = left[-1]
            else: pos = c[-1]; idx += 1
            used[pos] = True; tot += idx; k += 1
    return tot
@functools.lru_cache(None)
def KF(lam, mu):
    tot = 0
    for T in ssyt(lam, mu):
        word = [x for r in reversed(T) for x in r]
        tot += t**charge_word(word)
    return sp.expand(tot)
@functools.lru_cache(None)
def X(lam, rho):
    n = sum(lam)
    return sp.expand(sum(chi(mu, rho)*KF(mu, lam) for mu in parts(n)))
if __name__ == '__main__':
    bad = 0
    for n in range(1, 8):
        for lam in parts(n):
            l = len(lam)
            ref = sp.expand(t**nfun(lam)*sp.prod([1-t**(-i) for i in range(1, l)]))
            if sp.expand(X(lam, (n,))-ref) != 0: bad += 1; print('one-part BAD', lam)
            for rho in parts(n):
                if X(lam, rho).subs(t, 0) != chi(lam, rho): bad += 1; print('t=0 BAD', lam, rho)
    print('Green engine sanity n<=7: bad =', bad)
