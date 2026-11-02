#!/usr/bin/env python3
r"""Seed Q1 test case (2026-11-02): Sym coproduct vs. type-A tableau crystal.

Part A  gl_n (n=2,3), |lambda|<=4:
  A1  Levi split x_1..x_k | x_{k+1}..x_n : expand s_lambda(x,y) = sum d_{mu,nu} s_mu(x) s_nu(y)
      (pure Sym computation, sympy) and compare with the crystal B(lambda) restricted to the
      colours I\{k}: highest-weight elements for those colours, weights split as (mu,nu).
      Also count f_k-edges between Levi components and report the multiplicities d (the
      number of Levi-structure-preserving relabellings is prod d!).
  A2  Tensor product B(mu) (x) B(nu) with the signature rule: count highest-weight
      elements of weight lambda and compare with c^lambda_{mu nu} from s_mu s_nu (sympy).
Part B  sl_2: how many "local" normal crystal structures exist on the SET B(m) x B(n)
      (each f-edge changes exactly one factor by one step) with the right character?
      Versus the 2 signature conventions.
Part C  sl_2: uniform rules  f(b(x)b') = f b (x) b'  iff  C(phi(b), eps(b')) , else b (x) f b',
      C an arbitrary boolean function on {0..R}^2.  Which C give normal crystals on all
      B(m)(x)B(n), m,n<=R, AND make the identity map (B1(x)B2)(x)B3 -> B1(x)(B2(x)B3) a crystal
      isomorphism (strict associativity) for all m1,m2,m3 <= R3?
All outputs are 'computed' grade.
"""
import itertools, math
from collections import defaultdict, Counter
import sympy as sp

# ---------------------------------------------------------------- tableaux & crystal
def partitions(N, maxlen):
    def rec(n, mx, ln):
        if n == 0:
            yield ()
            return
        if ln == 0:
            return
        for p in range(min(n, mx), 0, -1):
            for rest in rec(n - p, p, ln - 1):
                yield (p,) + rest
    return list(rec(N, N, maxlen))

def ssyt(shape, n):
    """all SSYT of shape with entries 1..n, as tuple of row tuples"""
    cells = [(r, c) for r, L in enumerate(shape) for c in range(L)]
    res = []
    T = {}
    def rec(idx):
        if idx == len(cells):
            res.append(tuple(tuple(T[(r, c)] for c in range(shape[r])) for r in range(len(shape))))
            return
        r, c = cells[idx]
        lo = 1
        if c > 0: lo = max(lo, T[(r, c - 1)])
        if r > 0: lo = max(lo, T[(r - 1, c)] + 1)
        for v in range(lo, n + 1):
            T[(r, c)] = v
            rec(idx + 1)
        T.pop((r, c), None)
    rec(0)
    return res

def reading_word(T):
    w = []
    for row in reversed(T):
        w.extend(row)
    return w

def word_f(w, i):
    """LS bracket rule: cancel (i+1 ... i) pairs; f_i changes rightmost unmatched i."""
    stack = []          # positions of unmatched i+1
    unmatched_i = []
    for p, a in enumerate(w):
        if a == i + 1:
            stack.append(p)
        elif a == i:
            if stack:
                stack.pop()
            else:
                unmatched_i.append(p)
    if not unmatched_i:
        return None
    p = unmatched_i[-1]
    w2 = list(w); w2[p] = i + 1
    return w2

def word_e(w, i):
    unmatched_ip1 = []
    pend_i = 0
    # scan right to left: i's are '(' , i+1 are ')'
    for p in range(len(w) - 1, -1, -1):
        a = w[p]
        if a == i:
            pend_i += 1
        elif a == i + 1:
            if pend_i:
                pend_i -= 1
            else:
                unmatched_ip1.append(p)
    if not unmatched_ip1:
        return None
    p = unmatched_ip1[-1]   # leftmost unmatched i+1 (scanned last)
    w2 = list(w); w2[p] = i
    return w2

def tab_from_word(w, shape):
    rows = []
    pos = 0
    for L in reversed(shape):
        rows.append(tuple(w[pos:pos + L])); pos += L
    return tuple(reversed(rows))

def is_ssyt(T):
    for r, row in enumerate(T):
        for c in range(len(row)):
            if c and row[c] < row[c - 1]: return False
            if r and row[c] <= T[r - 1][c]: return False
    return True

def tab_f(T, i, shape):
    w = word_f(reading_word(T), i)
    if w is None: return None
    T2 = tab_from_word(w, shape); assert is_ssyt(T2); return T2

def tab_e(T, i, shape):
    w = word_e(reading_word(T), i)
    if w is None: return None
    T2 = tab_from_word(w, shape); assert is_ssyt(T2); return T2

def weight(entries, n):
    c = Counter(entries)
    return tuple(c[a] for a in range(1, n + 1))

# ---------------------------------------------------------------- Sym side (sympy)
def schur_poly(lam, xs):
    n = len(xs)
    lam = list(lam) + [0] * (n - len(lam))
    if len([p for p in lam if p > 0]) > n:
        return sp.Integer(0)
    M = sp.Matrix(n, n, lambda i, j: xs[i] ** (lam[j] + n - 1 - j))
    V = sp.Matrix(n, n, lambda i, j: xs[i] ** (n - 1 - j))
    return sp.cancel(M.det() / V.det())

def expand_in_schur(poly, xs, maxdeg):
    """peel leading (lex) monomials: returns {lam: coeff}"""
    n = len(xs)
    P = sp.Poly(sp.expand(poly), *xs)
    out = {}
    while not P.is_zero:
        mon, coeff = P.terms()[0]           # lex-leading term = dominant weight
        lam = tuple(m for m in mon if m > 0)
        assert list(mon) == sorted(mon, reverse=True), mon
        out[lam] = coeff
        P = P - sp.Poly(coeff * schur_poly(lam, xs), *xs)
    return out

def levi_coproduct(lam, n, k):
    """s_lam(x_1..x_n) = sum d s_mu(x_1..x_k) s_nu(x_{k+1}..x_n); peel by bi-leading terms."""
    xs = sp.symbols(f'x1:{n+1}')
    X, Y = xs[:k], xs[k:]
    P = sp.Poly(sp.expand(schur_poly(lam, xs)), *xs)
    out = {}
    while not P.is_zero:
        mon, coeff = P.terms()[0]
        mu = tuple(m for m in mon[:k] if m > 0); nu = tuple(m for m in mon[k:] if m > 0)
        assert list(mon[:k]) == sorted(mon[:k], reverse=True)
        assert list(mon[k:]) == sorted(mon[k:], reverse=True)
        out[(mu, nu)] = coeff
        P = P - sp.Poly(coeff * sp.expand(schur_poly(mu, X) * schur_poly(nu, Y)), *xs)
    return out

# ---------------------------------------------------------------- Part A
def partA():
    print("=== Part A1: Levi coproduct of Sym vs Levi restriction of B(lambda) ===")
    allok = True
    for n in (2, 3):
        for N in range(1, 5):
            for lam in partitions(N, n):
                B = ssyt(lam, n)
                for k in range(1, n):
                    sym = levi_coproduct(lam, n, k)
                    colours = [i for i in range(1, n) if i != k]
                    hw = [T for T in B if all(tab_e(T, i, lam) is None for i in colours)]
                    cr = Counter()
                    for T in hw:
                        wt = weight([a for r in T for a in r], n)
                        mu = tuple(a for a in wt[:k] if a > 0); nu = tuple(a for a in wt[k:] if a > 0)
                        cr[(mu, nu)] += 1
                    ok = dict(cr) == {key: int(v) for key, v in sym.items()}
                    allok &= ok
                    # components and f_k edges between components
                    comp = {}
                    for T in B:
                        U = T
                        while True:
                            moved = False
                            for i in colours:
                                V = tab_e(U, i, lam)
                                if V is not None: U = V; moved = True; break
                            if not moved: break
                        comp[T] = U
                    edges = sum(1 for T in B if tab_f(T, k, lam) is not None)
                    mult = [int(v) for v in sym.values()]
                    relabel = math.prod(math.factorial(m) for m in mult)
                    wtspaces = Counter(weight([a for r in T for a in r], n) for T in B)
                    amb = sum(1 for v in wtspaces.values() if v >= 2)
                    print(f"n={n} lam={lam} k={k}: |B|={len(B)} Sym-Delta={dict((key,int(v)) for key,v in sym.items())} "
                          f"crystal-match={ok}  f_k-edges={edges}  prod(d!)={relabel}  weight-spaces-dim>=2: {amb}")
    print("A1 all match:", allok)

    print("\n=== Part A2: highest weights of B(mu)(x)B(nu) (signature rule) vs c^lam_{mu nu} ===")
    allok = True
    for n in (2, 3):
        xs = sp.symbols(f'x1:{n+1}')
        for N1 in range(1, 3):
            for N2 in range(1, 3):
                for mu in partitions(N1, n):
                    for nu in partitions(N2, n):
                        sym = expand_in_schur(schur_poly(mu, xs) * schur_poly(nu, xs), xs, N1 + N2)
                        cr = Counter()
                        for T in ssyt(mu, n):
                            for U in ssyt(nu, n):
                                w = reading_word(T) + reading_word(U)
                                if all(word_e(w, i) is None for i in range(1, n)):
                                    cr[tuple(a for a in weight(w, n) if a > 0)] += 1
                        ok = dict(cr) == {key: int(v) for key, v in sym.items()}
                        allok &= ok
                        print(f"n={n} {mu}x{nu}: Sym={dict((k_,int(v)) for k_,v in sym.items())} crystal-hw={dict(cr)} {ok}")
    print("A2 all match:", allok)

# ---------------------------------------------------------------- Part B (sl_2 local structures)
def normal_sl2(V, f):
    """V: dict vertex->weight; f: dict vertex->vertex or None. normal iff f injective partial map,
    acyclic, and every string has weights k, k-2, ..., -k."""
    tgt = [f[v] for v in V if f[v] is not None]
    if len(tgt) != len(set(tgt)): return False
    has_pre = set(tgt)
    seen = set()
    for v in V:
        if v in has_pre: continue
        s = [v]
        while f[s[-1]] is not None:
            s.append(f[s[-1]])
            if len(s) > len(V): return False
        top = V[s[0]]
        if top < 0 or len(s) != top + 1: return False
        if any(V[x] != top - 2 * j for j, x in enumerate(s)): return False
        seen.update(s)
    return len(seen) == len(V)

def count_local_structures(m, n):
    """vertices (a,b), 0<=a<=m, 0<=b<=n (a = # f applied in B(m)); weight m+n-2a-2b.
    local: f(a,b) in {(a+1,b),(a,b+1),None}. Count normal structures (character auto-correct)."""
    V = {(a, b): m + n - 2 * a - 2 * b for a in range(m + 1) for b in range(n + 1)}
    verts = list(V)
    count = 0
    sols = []
    def opts(v):
        a, b = v
        o = [None]
        if a < m: o.append((a + 1, b))
        if b < n: o.append((a, b + 1))
        return o
    # brute force with pruning: vertices of nonpositive weight... just brute force small grids
    for choice in itertools.product(*[opts(v) for v in verts]):
        f = dict(zip(verts, choice))
        if normal_sl2(V, f):
            count += 1; sols.append(f)
    return count, sols

def signature_f(m, n, kashiwara=True):
    """sl_2 B(m)(x)B(n), element (a,b): phi(b1)=m-a, eps(b1)=a, phi(b2)=n-b, eps(b2)=b.
    Kashiwara convention: f acts on first iff phi(b1) > eps(b2). anti: on second iff phi(b2)>eps(b1)."""
    f = {}
    for a in range(m + 1):
        for b in range(n + 1):
            if kashiwara:
                if (m - a) > b: f[(a, b)] = (a + 1, b)
                else: f[(a, b)] = (a, b + 1) if b < n else None
            else:
                if (n - b) > a: f[(a, b)] = (a, b + 1)
                else: f[(a, b)] = (a + 1, b) if a < m else None
    return f

def partB():
    print("\n=== Part B: sl_2 local normal crystal structures on the set B(m) x B(n) ===")
    for m in range(1, 4):
        for n in range(1, 4):
            if (m + 1) * (n + 1) > 12: continue
            c, sols = count_local_structures(m, n)
            K = signature_f(m, n, True); A = signature_f(m, n, False)
            inK = any(s == K for s in sols); inA = any(s == A for s in sols)
            print(f"B({m})xB({n}): local normal structures = {c}; Kashiwara-rule among them: {inK}; anti-Kashiwara: {inA}; K==anti: {K==A}")

# ---------------------------------------------------------------- Part C (uniform rules)
class Cr:
    """finite sl_2 crystal: elements list, wt, eps, phi, f (dict or None)"""
    def __init__(self, elts, wt, f):
        self.elts = elts; self.wt = wt; self.f = f
        self.e = {v: None for v in elts}
        for v in elts:
            if f[v] is not None: self.e[f[v]] = v
        self.eps = {}; self.phi = {}
        for v in elts:
            k = 0; u = v
            while self.e[u] is not None: u = self.e[u]; k += 1
            self.eps[v] = k
            k = 0; u = v
            while self.f[u] is not None: u = self.f[u]; k += 1
            self.phi[v] = k

def Bm(m):
    elts = list(range(m + 1))
    return Cr(elts, {j: m - 2 * j for j in elts}, {j: (j + 1 if j < m else None) for j in elts})

def tensor(B1, B2, C):
    elts = [(x, y) for x in B1.elts for y in B2.elts]
    f = {}
    for (x, y) in elts:
        if C(B1.phi[x], B2.eps[y]):
            fx = B1.f[x]; f[(x, y)] = (fx, y) if fx is not None else None
        else:
            fy = B2.f[y]; f[(x, y)] = (x, fy) if fy is not None else None
    tgt = [v for v in f.values() if v is not None]
    if len(tgt) != len(set(tgt)): return None
    wt = {(x, y): B1.wt[x] + B2.wt[y] for (x, y) in elts}
    if not normal_sl2(wt, f): return None
    return Cr(elts, wt, f)

def partC(R=3, R3=2):
    print(f"\n=== Part C: uniform rules C(phi(b),eps(b')) on {{0..{R}}}^2 ===")
    cells = [(p, e) for p in range(R + 1) for e in range(R + 1)]
    normal_rules = []; assoc_rules = []
    Bs = {m: Bm(m) for m in range(R + 1)}
    for bits in range(1 << len(cells)):
        tab = {cells[j]: bool(bits >> j & 1) for j in range(len(cells))}
        C = lambda p, e: tab[(p, e)]
        ok = True
        for m in range(R + 1):
            for n in range(R + 1):
                if m + n == 0: continue
                if tensor(Bs[m], Bs[n], C) is None: ok = False; break
            if not ok: break
        if not ok: continue
        normal_rules.append(tab)
        # associativity on triples with m_i <= R3 (phi/eps of pairs stay <= 2*R3 <= R needed)
        good = True
        for m1, m2, m3 in itertools.product(range(R3 + 1), repeat=3):
            if m1 + m2 > R or m2 + m3 > R: continue
            L = tensor(tensor(Bs[m1], Bs[m2], C), Bs[m3], C)
            Rr = tensor(Bs[m1], tensor(Bs[m2], Bs[m3], C), C)
            if L is None or Rr is None: good = False; break
            for ((x, y), z) in L.elts:
                a = L.f[((x, y), z)]
                b = Rr.f[(x, (y, z))]
                a2 = None if a is None else (a[0][0], (a[0][1], a[1]))
                if a2 != b: good = False; break
            if not good: break
        if good: assoc_rules.append(tab)
    def show(tab):
        return "\n".join("   " + " ".join('F' if tab[(p, e)] else 's' for e in range(R + 1)) + f"   phi={p}" for p in range(R + 1))
    print(f"rules giving normal B(m)(x)B(n) for all m,n<= {R}: {len(normal_rules)}")
    print(f"... of which strictly associative on triples (m_i<= {R3}, pair sums <= {R}): {len(assoc_rules)}")
    sig = {(p, e): p > e for (p, e) in cells}
    for tab in assoc_rules:
        print(" rule (rows phi(b)=0..R, cols eps(b')=0..R; F = act on first):  signature? ", tab == sig)
        print(show(tab))
    # how many normal rules agree with signature on the 'forced' region?
    diff = Counter()
    for tab in normal_rules:
        for c in cells:
            if tab[c] != sig[c]: diff[c] += 1
    print(" cells where some normal rule differs from signature:", sorted(diff.items()))

if __name__ == "__main__":
    partA(); partB(); partC()
