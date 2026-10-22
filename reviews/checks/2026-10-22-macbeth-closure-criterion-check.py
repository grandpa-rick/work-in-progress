"""Rick's check of MacBeth closure-criterion note (4ab502e), Thm 2. Exact finite enumeration."""
import itertools
from collections import Counter

# ---------- Part A: monoids of order <=3, normalizers (K1),(U),(LA) ----------
def monoids(n):
    els = range(n)  # 0 = unit
    others = [(i, j) for i in range(1, n) for j in range(1, n)]
    for vals in itertools.product(els, repeat=len(others)):
        T = [[0]*n for _ in els]
        for i in els: T[0][i] = i; T[i][0] = i
        for (i, j), v in zip(others, vals): T[i][j] = v
        if all(T[T[a][b]][c] == T[a][T[b][c]] for a in els for b in els for c in els):
            yield T

def normalizers(T):
    n = len(T)
    for K in itertools.product(range(n), repeat=n):
        if K[0] != 0: continue
        if any(K[K[a]] != K[a] for a in range(n)): continue
        if any(K[T[K[x]][c]] != K[T[x][c]] for x in range(n) for c in range(n)): continue
        yield K

def analyse(T, K):
    n = len(T); I = [a for a in range(n) if K[a] == a]
    closed = all(K[T[f][g]] == T[f][g] for f in I for g in I)
    sh = lambda f, g: K[T[f][g]]
    assoc = all(sh(sh(f, g), h) == sh(f, sh(g, h)) for f in I for g in I for h in I)
    return closed, assoc

print("== Part A: order-3 witness ==")
# e=0,a=1,b=2 ; a*a=b, a*b=a, b*a=a, b*b=b
M = [[0,1,2],[1,2,1],[2,1,2]]
assert all(M[M[x][y]][z] == M[x][M[y][z]] for x in range(3) for y in range(3) for z in range(3))
K = (0,1,1)
assert K in list(normalizers(M))
print("witness (closed, assoc) =", analyse(M, K), " expect (False, True)")
for n in (1, 2, 3):
    tabs = list(monoids(n)); wit = 0; tabsw = 0; handed = 0
    for T in tabs:
        w = [K for K in normalizers(T) if analyse(T, K) == (False, True)]
        handed += sum(1 for K in normalizers(T) if not analyse(T, K)[1])
        wit += len(w); tabsw += bool(w)
    print(f"order {n}: {len(tabs)} labelled monoid tables, {wit} non-closed-but-assoc normalizers on {tabsw} tables, handed={handed}")

# ---------- Part B: groups as permutation groups ----------
def mul(p, q):  # diagrammatic: p then q
    return tuple(q[p[i]] for i in range(len(p)))
def inv(p):
    r = [0]*len(p)
    for i, x in enumerate(p): r[x] = i
    return tuple(r)
def gen(gens, n):
    e = tuple(range(n)); G = {e}; fr = [e]
    while fr:
        nf = []
        for g in fr:
            for s in gens:
                h = mul(g, s)
                if h not in G: G.add(h); nf.append(h)
        fr = nf
    return frozenset(G)
def subgroups(G, n):
    Gl = sorted(G); S = set()
    for a in Gl:
        for b in Gl:
            S.add(gen([a, b], n))
    # add 3-generated for safety (small groups)
    S2 = set(S)
    for H in S:
        for c in Gl: S2.add(gen(list(H) + [c], n))
    return S2
def is_normal(H, G):
    return all(mul(mul(inv(g), h), g) in H for g in G for h in H)

def perm_from_cycles(n, cycles):
    p = list(range(n))
    for c in cycles:
        for i in range(len(c)): p[c[i]] = c[(i+1) % len(c)]
    return tuple(p)

def direct(g1, n1, g2, n2):
    out = []
    for p in g1: out.append(tuple(list(p) + [n1 + i for i in range(n2)]))
    for q in g2: out.append(tuple(list(range(n1)) + [n1 + x for x in q]))
    return out, n1 + n2

groups = {}
c = lambda n, *cy: perm_from_cycles(n, cy)
groups["C4"] = ([c(4, (0,1,2,3))], 4)
groups["C2xC2"] = ([c(4, (0,1)), c(4, (2,3))], 4)
groups["S3"] = ([c(3, (0,1,2)), c(3, (0,1))], 3)
groups["D4"] = ([c(4, (0,1,2,3)), c(4, (0,2))], 4)
# Q8 regular rep: elements 1,i,j,k,-1,-i,-j,-k as 0..7
def qmul_table():
    import numpy as np
    basis = {'1': np.eye(2, dtype=complex), 'i': np.array([[1j,0],[0,-1j]]), 'j': np.array([[0,1],[-1,0]], dtype=complex)}
    basis['k'] = basis['i'] @ basis['j']
    els = [basis['1'], basis['i'], basis['j'], basis['k']]; els = els + [-x for x in els]
    def idx(m):
        for t, x in enumerate(els):
            if abs(m - x).max() < 1e-9: return t
    return els, idx
els, idx = qmul_table()
qi = tuple(idx(x @ els[1]) for x in els); qj = tuple(idx(x @ els[2]) for x in els)  # right mult = "x then i"
groups["Q8"] = ([qi, qj], 8)
groups["C8"] = ([c(8, (0,1,2,3,4,5,6,7))], 8)
groups["C2xC4"] = ([c(6, (0,1)), c(6, (2,3,4,5))], 6)
groups["A4"] = ([c(4, (0,1,2)), c(4, (0,1), (2,3))], 4)
s3 = [c(3, (0,1,2)), c(3, (0,1))]; c2 = [c(2, (0,1))]
groups["S3xC2"] = direct(s3, 3, c2, 2)
groups["S4"] = ([c(4, (0,1,2,3)), c(4, (0,1))], 4)

def right_cosets(H, G):
    seen = set(); cos = []
    for g in sorted(G):
        if g in seen: continue
        Hg = frozenset(mul(h, g) for h in H)  # h then g : coset H g
        seen |= Hg; cos.append(Hg)
    return cos

print("\n== Part B: right transversals T (1 in T) of all subgroups H ==")
print("columns: group | #(H,T) | (LA) fails | assoc&T subgroup | assoc&not subgroup [normal H / NON-normal H] | non-assoc")
nonnormal_examples = []
crit_mismatch = 0
for name, (gens, n) in groups.items():
    G = gen(gens, n); e = tuple(range(n))
    tot = la_fail = sub = cx_n = cx_nn = nonassoc = 0
    per_H = Counter()
    for H in subgroups(G, n):
        normal = is_normal(H, G)
        cos = right_cosets(H, G)
        choices = [[e] if e in C else sorted(C) for C in cos]
        for Tt in itertools.product(*choices):
            T = set(Tt); rep = {}
            for C, t in zip(cos, Tt):
                for g in C: rep[g] = t
            r = lambda g: rep[g]
            tot += 1
            if any(r(mul(r(x), y)) != r(mul(x, y)) for x in G for y in G): la_fail += 1
            o = lambda a, b: r(mul(a, b))
            assoc = all(o(o(a, b), d) == o(a, o(b, d)) for a in T for b in T for d in T)
            issub = all(mul(a, b) in T for a in T for b in T)
            # exact criterion derived from (LA)+(RA#): f phi(g,h) f^-1 in H, phi = gh * r(gh)^-1
            phis = {mul(mul(a, b), inv(r(mul(a, b)))) for a in T for b in T}
            crit = all(mul(mul(f, p), inv(f)) in H for f in T for p in phis)
            if crit != assoc: crit_mismatch += 1
            if not assoc: nonassoc += 1
            elif issub: sub += 1
            elif normal:
                cx_n += 1; per_H[len(H)] += 1
            else:
                cx_nn += 1
                if len(nonnormal_examples) < 3 or (name == "S3xC2" and len(nonnormal_examples) < 6):
                    core = frozenset.intersection(*[frozenset(mul(mul(inv(g), h), g) for h in H) for g in G])
                    nonnormal_examples.append((name, len(H), len(core), all(p in core for p in phis)))
    print(f"{name:6s} | {tot:5d} | {la_fail} | {sub:4d} | {cx_n:4d} / {cx_nn:4d} | {nonassoc}   (normal-H counterex. by |H|: {dict(per_H)})")
print("criterion [assoc <=> f*phi*f^-1 in H for all f in T, phi in defects] mismatches:", crit_mismatch)
print("non-normal-H associative non-subgroup examples (group,|H|,|core H|, defects in core?):", nonnormal_examples)

# ---------- Part C: the paper's Q8 instance and the C4 instance; non-split check ----------
print("\n== Part C ==")
# Q8: H={1,-1} = {0,4}; T = {1,-i,-j,-k} = {0,5,6,7}; (-i)(-j)=ij=k
print("(-i)(-j) =", idx(els[5] @ els[6]), "(index 3 = k, not in T={0,5,6,7})")
# non-split: does any transversal of Z(Q8), Z(D4), 2Z/4 form a subgroup?
for name, (gens, n), Hdesc in [("Q8", groups["Q8"], 2), ("D4", groups["D4"], 2), ("C4", groups["C4"], 2)]:
    G = gen(gens, n)
    Z = frozenset(z for z in G if all(mul(z, g) == mul(g, z) for g in G))
    Hs = [H for H in subgroups(G, n) if len(H) == 2 and is_normal(H, G)]
    for H in Hs:
        cos = right_cosets(H, G); e = tuple(range(n))
        choices = [[e] if e in C else sorted(C) for C in cos]
        nsub = sum(all(mul(a, b) in set(Tt) for a in Tt for b in Tt) for Tt in itertools.product(*choices))
        print(f"{name}: normal H of order 2 {'(=centre)' if H <= Z else ''}: {len(choices and list(itertools.product(*choices)))} transversals, {nsub} are subgroups (0 <=> extension non-split)")
