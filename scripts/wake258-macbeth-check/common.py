"""Shared finite-structure helpers for the uid456 check pass (Rick, Wake 258).
Everything is re-implemented from the printed definitions in MacBeth's PDFs (content 6df0043)."""
from itertools import product

# ---------- monoids with normalizer: definitions as printed (mirror-LA Def 1.1, Prop 1.2, Def 1.3, Sec 5) ----------
def is_assoc(t, n):
    return all(t[t[a][b]][c] == t[a][t[b][c]] for a in range(n) for b in range(n) for c in range(n))

def is_normalizer(t, n, e, K, require_LA=True):
    if any(K[K[a]] != K[a] for a in range(n)): return False          # (K1)
    if K[e] != e: return False                                         # (U)
    if require_LA and any(K[t[K[x]][c]] != K[t[x][c]] for x in range(n) for c in range(n)): return False  # (LA)
    return True

def image(K, n): return [a for a in range(n) if K[a] == a]

def norm_assoc(t, n, K):
    I = image(K, n); s = lambda f, g: K[t[f][g]]
    return all(s(s(f, g), h) == s(f, s(g, h)) for f in I for g in I for h in I)

def RAsharp(t, n, K):
    I = image(K, n)
    return all(K[t[f][K[t[g][h]]]] == K[t[t[f][g]][h]] for f in I for g in I for h in I)

def mLA_full(t, n, K):
    return all(K[t[c][K[x]]] == K[t[c][x]] for c in range(n) for x in range(n))

def mLA_restricted(t, n, K):          # Prop 5.1 (2): all c in A, x in ImK.ImK
    I = image(K, n); X = {t[g][h] for g in I for h in I}
    return all(K[t[c][K[x]]] == K[t[c][x]] for c in range(n) for x in X)

def theta_explicit(t, n, K):          # Sec 5 printed formula: a th b iff K(x a y) = K(x b y) for all x,y
    return [[all(K[t[t[x][a]][y]] == K[t[t[x][b]][y]] for x in range(n) for y in range(n)) for b in range(n)] for a in range(n)]

def largest_congruence_in(t, n, K):   # independent algorithm: refine ~K until stable under left/right mult
    cls = [K[a] for a in range(n)]
    while True:
        sig = [(cls[a],) + tuple(cls[t[x][a]] for x in range(n)) + tuple(cls[t[a][y]] for y in range(n)) for a in range(n)]
        lab = {}; new = [lab.setdefault(s, len(lab)) for s in sig]
        if len(set(new)) == len(set(cls)): return [[new[a] == new[b] for b in range(n)] for a in range(n)]
        cls = new

def p51_three(t, n, K, th):           # (3): K(gh) theta_K gh for g,h in Im K
    I = image(K, n)
    return all(th[K[t[g][h]]][t[g][h]] for g in I for h in I)

def closed(t, n, K):
    I = set(image(K, n)); return all(t[f][g] in I for f in I for g in I)

def monoids_unit0(n):
    """All labelled monoid tables on {0..n-1} with unit 0 (backtracking with associativity pruning)."""
    t = [[None] * n for _ in range(n)]
    for a in range(n): t[0][a] = a; t[a][0] = a
    cells = [(a, b) for a in range(1, n) for b in range(1, n)]
    out = []
    def ok():
        for a in range(n):
            for b in range(n):
                ab = t[a][b]
                if ab is None: continue
                for c in range(n):
                    bc = t[b][c]
                    if bc is None: continue
                    l = t[ab][c]; r = t[a][bc]
                    if l is not None and r is not None and l != r: return False
        return True
    def rec(i):
        if i == len(cells): out.append([row[:] for row in t]); return
        a, b = cells[i]
        for v in range(n):
            t[a][b] = v
            if ok(): rec(i + 1)
        t[a][b] = None
    rec(0); return out

# ---------- finite groups as multiplication tables ----------
def table_from(elems, mul):
    idx = {x: i for i, x in enumerate(elems)}
    return [[idx[mul(a, b)] for b in elems] for a in elems]

def cyclic(m): return table_from(list(range(m)), lambda a, b: (a + b) % m)
def dihedral(m):   # (k,e): r^k s^e, s r = r^-1 s ; order 2m
    el = [(k, e) for e in (0, 1) for k in range(m)]
    return table_from(el, lambda x, y: ((x[0] + (y[0] if x[1] == 0 else -y[0])) % m, (x[1] + y[1]) % 2))
def dicyclic(m):   # order 4m: a^{2m}=1, x^2=a^m, x a x^-1 = a^-1
    el = [(k, e) for e in (0, 1) for k in range(2 * m)]
    def mul(x, y):
        k, e = x; l, f = y
        if e == 0: return ((k + l) % (2 * m), f)
        if f == 0: return ((k - l) % (2 * m), 1)
        return ((k - l + m) % (2 * m), 0)
    return table_from(el, mul)
def direct(A, B):
    na, nb = len(A), len(B); el = [(i, j) for i in range(na) for j in range(nb)]
    return table_from(el, lambda x, y: (A[x[0]][y[0]], B[x[1]][y[1]]))
def perm_group(gens, deg):
    e = tuple(range(deg)); els = {e}; frontier = [e]
    comp = lambda p, q: tuple(q[p[i]] for i in range(deg))   # p then q
    while frontier:
        nf = []
        for p in frontier:
            for g in gens:
                r = comp(p, g)
                if r not in els: els.add(r); nf.append(r)
        frontier = nf
    els = sorted(els); return table_from(els, comp)
def alt4(): return perm_group([(1, 2, 0, 3), (1, 0, 3, 2)], 4)
def sym(n):
    return perm_group([tuple([1, 0] + list(range(2, n))), tuple(list(range(1, n)) + [0])], n)

def ident(G): return next(e for e in range(len(G)) if all(G[e][x] == x for x in range(len(G))))
def inverses(G, e): return [next(y for y in range(len(G)) if G[x][y] == e) for x in range(len(G))]

def subgroups(G):
    n = len(G); e = ident(G)
    def gen(S):
        H = {e} | set(S); fr = list(H)
        while fr:
            nf = []
            for a in fr:
                for b in list(H):
                    for c in (G[a][b], G[b][a]):
                        if c not in H: H.add(c); nf.append(c)
            fr = nf
        return frozenset(H)
    subs = {frozenset([e])}; fr = [frozenset([e])]
    while fr:
        nf = []
        for H in fr:
            for g in range(n):
                if g not in H:
                    K = gen(set(H) | {g})
                    if K not in subs: subs.add(K); nf.append(K)
        fr = nf
    return sorted(subs, key=lambda s: (len(s), sorted(s)))

def is_normal(G, H, inv):
    return all(G[G[g][h]][inv[g]] in H for g in range(len(G)) for h in H)
def core(G, H, inv):
    return frozenset(h for h in H if all(G[G[g][h]][inv[g]] in H for g in range(len(G))))
def right_cosets(G, H, e):
    seen = set(); cos = []
    for g in [e] + [x for x in range(len(G)) if x != e]:
        if g in seen: continue
        c = frozenset(G[h][g] for h in H); seen |= c; cos.append(c)
    return cos   # cos[0] = H
def transversals(G, H, e):
    cos = right_cosets(G, H, e)
    for choice in product(*[sorted(c) for c in cos[1:]]):
        yield (e,) + choice
def rep_map(G, H, T):
    K = [None] * len(G)
    for t in T:
        for h in H: K[G[h][t]] = t
    return K
