#!/usr/bin/env python3
"""Independent check of MacBeth uid449 'The mirror law (mLA)' (content 70b090c).
Re-implemented from the printed statements, not from MacBeth's scripts.

Part A: all monoids of order <= 3 (labelled: any element may be the unit), all normalizers
        K (K1, U, LA).  Tabulate (mLA), (RA#), #-assoc, closure, left-congruence.
Part A': idempotent unit-fixing K WITHOUT (LA): does (mLA) alone force #-associativity?
Part B: the two printed separators of Thm 3.1, verbatim, and the corrected second table.
Part C: groups S3, S4, S3xC2: all (H,T), 1 in T.  (LA); (mLA) <=> H normal;
        assoc <=> f phi f^-1 in H; criterion <=> phi(T,T) subset core_G(H).
Exact finite enumeration only.
"""
import itertools, sys

# ---------------- Part A ----------------
def monoids(n):
    els = range(n)
    for flat in itertools.product(els, repeat=n * n):
        t = [flat[i * n:(i + 1) * n] for i in range(n)]
        units = [e for e in els if all(t[e][x] == x and t[x][e] == x for x in els)]
        if not units:
            continue
        if all(t[t[a][b]][c] == t[a][t[b][c]] for a in els for b in els for c in els):
            yield t, units[0]

def analyse(t, e, K, n):
    els = range(n)
    m = lambda a, b: t[a][b]
    img = [a for a in els if K[a] == a]
    LA = all(K[m(K[x], c)] == K[m(x, c)] for x in els for c in els)
    mLA = all(K[m(c, K[x])] == K[m(c, x)] for x in els for c in els)
    RA = all(K[m(f, K[m(g, h)])] == K[m(m(f, g), h)] for f in img for g in img for h in img)
    sh = lambda f, g: K[m(f, g)]
    assoc = all(sh(sh(f, g), h) == sh(f, sh(g, h)) for f in img for g in img for h in img)
    closed = all(K[m(f, g)] == m(f, g) for f in img for g in img)
    lcong = all(K[m(c, a)] == K[m(c, b)] for a in els for b in els if K[a] == K[b] for c in els)
    rcong = all(K[m(a, c)] == K[m(b, c)] for a in els for b in els if K[a] == K[b] for c in els)
    return dict(LA=LA, mLA=mLA, RA=RA, assoc=assoc, closed=closed, lcong=lcong, rcong=rcong)

print("=== Part A: monoids of order <= 3, normalizers (K1,U,LA) ===")
for n in (1, 2, 3):
    nmon = 0
    stats = dict(norm=0, mLA=0, mLA_notRA=0, mLA_notassoc=0, RA_ne_assoc=0,
                 mLA_ne_lcong=0, LA_ne_rcong=0, sep_RA_notmLA=0,
                 sep_closed_notmLA=0, sep_mLA_notclosed=0, assoc_neither=0)
    seps = []
    noLA = dict(K=0, mLA=0, mLA_notassoc=0)
    for t, e in monoids(n):
        nmon += 1
        for K in itertools.product(range(n), repeat=n):
            if K[e] != e or any(K[K[a]] != K[a] for a in range(n)):
                continue
            r = analyse(t, e, K, n)
            if not r['LA']:
                noLA['K'] += 1
                if r['mLA']:
                    noLA['mLA'] += 1
                    if not r['assoc']:
                        noLA['mLA_notassoc'] += 1
                continue
            stats['norm'] += 1
            stats['mLA'] += r['mLA']
            stats['mLA_notRA'] += r['mLA'] and not r['RA']
            stats['mLA_notassoc'] += r['mLA'] and not r['assoc']
            stats['RA_ne_assoc'] += r['RA'] != r['assoc']
            stats['mLA_ne_lcong'] += r['mLA'] != r['lcong']
            stats['LA_ne_rcong'] += r['LA'] != r['rcong']
            if r['RA'] and not r['mLA']:
                stats['sep_RA_notmLA'] += 1
                seps.append((t, e, K))
            stats['sep_closed_notmLA'] += r['closed'] and not r['mLA']
            stats['sep_mLA_notclosed'] += r['mLA'] and not r['closed']
            stats['assoc_neither'] += r['assoc'] and not r['mLA'] and not r['closed']
    print(f"order {n}: {nmon} labelled monoids; normalizers stats {stats}")
    print(f"   idempotent unit-fixing K failing (LA): {noLA}")
    for s in seps:
        print("   separator (RA#) & not (mLA):", s)

# ---------------- Part B ----------------
print("\n=== Part B: printed separators of Thm 3.1 ===")
def check_table(name, t, e, K):
    n = len(t)
    ass = all(t[t[a][b]][c] == t[a][t[b][c]] for a in range(n) for b in range(n) for c in range(n))
    print(f"{name}: table associative = {ass}")
    if ass:
        print("   ", analyse(t, e, K, n))
    else:
        bad = [(a, b, c) for a in range(n) for b in range(n) for c in range(n)
               if t[t[a][b]][c] != t[a][t[b][c]]]
        print("    non-associative triples:", bad)
# separator 1: M={0,1,2}, unit 0, rows = left factor
t1 = [[0, 1, 2], [1, 1, 2], [2, 1, 2]]
check_table("Sep1 (printed)", t1, 0, (0, 0, 2))
# separator 2 as PRINTED in uid449: e=0,a=1,b=2; a.a=b, a.b=b.a=b.b=a
t2p = [[0, 1, 2], [1, 2, 1], [2, 1, 1]]
check_table("Sep2 (as printed: b.b=a)", t2p, 0, (0, 1, 1))
# separator 2 as in the closure-criterion note: b.b=b
t2c = [[0, 1, 2], [1, 2, 1], [2, 1, 2]]
check_table("Sep2 (corrected: b.b=b)", t2c, 0, (0, 1, 1))

# ---------------- Part C ----------------
print("\n=== Part C: groups ===")
def perm_group(gens, deg):
    idp = tuple(range(deg))
    G = {idp}; frontier = [idp]
    while frontier:
        new = []
        for g in frontier:
            for s in gens:
                h = tuple(s[g[i]] for i in range(deg))  # composition; any fixed convention
                if h not in G:
                    G.add(h); new.append(h)
        frontier = new
    return sorted(G)

def make_group(elems, mulf):
    idx = {g: i for i, g in enumerate(elems)}
    N = len(elems)
    mt = [[idx[mulf(a, b)] for b in elems] for a in elems]
    e = next(i for i in range(N) if all(mt[i][j] == j for j in range(N)))
    inv = [next(j for j in range(N) if mt[i][j] == e) for i in range(N)]
    return N, mt, e, inv

def subgroups(N, mt, e, inv):
    # all subgroups via closure of subsets generated by <=2 generators (enough for S4, S3xC2)
    subs = set()
    def gen(gs):
        S = {e}; fr = [e]
        while fr:
            nw = []
            for x in fr:
                for g in gs:
                    y = mt[x][g]
                    if y not in S:
                        S.add(y); nw.append(y)
            fr = nw
        return frozenset(S)
    for a in range(N):
        for b in range(N):
            subs.add(gen([a, b]))
    return sorted(subs, key=len)

def run_group(name, N, mt, e, inv):
    subs = subgroups(N, mt, e, inv)
    tot = LAfail = mLA_mm = crit_mm = core_mm = 0
    assoc_nonsub_nonnormal = 0
    for H in subs:
        if len(H) in (1, N):
            continue
        normal = all(mt[mt[g][h]][inv[g]] in H for g in range(N) for h in H)
        core = frozenset(x for x in H if all(mt[mt[inv[g]][x]][g] in H for g in range(N)))
        # right cosets Hg
        cosets = []; seen = set()
        for g in range(N):
            if g in seen: continue
            c = frozenset(mt[h][g] for h in H); seen |= c; cosets.append(sorted(c))
        cosets.sort(key=lambda c: e not in c)
        assert e in cosets[0]
        cos_of = {}
        for ci, c in enumerate(cosets):
            for x in c: cos_of[x] = ci
        for choice in itertools.product(*cosets[1:]):
            T = (e,) + choice
            tot += 1
            K = [T[cos_of[x]] for x in range(N)]
            LA = all(K[mt[K[x]][c]] == K[mt[x][c]] for x in range(N) for c in range(N))
            if not LA: LAfail += 1
            mLA = all(K[mt[c][K[x]]] == K[mt[c][x]] for x in range(N) for c in range(N))
            if mLA != normal: mLA_mm += 1
            sh = lambda f, g: K[mt[f][g]]
            assoc = all(sh(sh(f, g), h) == sh(f, sh(g, h)) for f in T for g in T for h in T)
            phi = {(g, h): mt[mt[g][h]][inv[sh(g, h)]] for g in T for h in T}
            assert all(p in H for p in phi.values())
            crit = all(mt[mt[f][phi[g, h]]][inv[f]] in H for f in T for g in T for h in T)
            corec = all(p in core for p in phi.values())
            if assoc != crit: crit_mm += 1
            if crit != corec: core_mm += 1
            Tset = set(T)
            issub = all(mt[a][b] in Tset for a in T for b in T)
            if assoc and not issub and not normal: assoc_nonsub_nonnormal += 1
    print(f"{name}: {len(subs)} subgroups; (H,T) pairs (H proper nontrivial) = {tot}; "
          f"(LA) fails {LAfail}; [(mLA)<=>normal] mismatches {mLA_mm}; "
          f"[assoc<=>criterion] mismatches {crit_mm}; [criterion<=>phi in core] mismatches {core_mm}; "
          f"assoc & non-subgroup & non-normal H: {assoc_nonsub_nonnormal}")

S3 = perm_group([(1, 0, 2), (1, 2, 0)], 3)
run_group("S3", *make_group(S3, lambda a, b: tuple(b[a[i]] for i in range(3))))
S3C2 = [(p, z) for p in S3 for z in (0, 1)]
run_group("S3xC2", *make_group(S3C2, lambda a, b: (tuple(b[0][a[0][i]] for i in range(3)), (a[1] + b[1]) % 2)))
S4 = perm_group([(1, 0, 2, 3), (1, 2, 3, 0)], 4)
run_group("S4", *make_group(S4, lambda a, b: tuple(b[a[i]] for i in range(4))))

# printed S3 witness: H=<(1 2)>, T=A3
N, mt, e, inv = make_group(S3, lambda a, b: tuple(b[a[i]] for i in range(3)))
idx = {g: i for i, g in enumerate(S3)}
H = {idx[(0, 1, 2)], idx[(1, 0, 2)]}
T = [idx[(0, 1, 2)], idx[(1, 2, 0)], idx[(2, 0, 1)]]
cos_of = {}
for t in T:
    for h in H: cos_of[mt[h][t]] = t
K = [cos_of[x] for x in range(N)]
assert sorted(cos_of) == list(range(N))
LA = all(K[mt[K[x]][c]] == K[mt[x][c]] for x in range(N) for c in range(N))
mLA = all(K[mt[c][K[x]]] == K[mt[c][x]] for x in range(N) for c in range(N))
closed = all(mt[a][b] in T for a in T for b in T)
print(f"S3 witness H=<(12)>, T=A3: LA={LA} mLA={mLA} closed={closed}")

# ---------------- Part D ----------------
# Referee's own candidate unifier (not in uid449): theta = largest two-sided congruence inside ~K
#   a theta b  <=>  K(x a y) = K(x b y) for all x, y.
# (C_theta): K(f.g) theta f.g for all f, g in Im K  ("Im K closed modulo theta").
# (mLA) <=> theta = ~K ; closure => (C_theta) trivially ; (C_theta) => (RA#) always (proof in review).
# Test: #-associative <=> (C_theta), over all monoids of order <= 4 (unit fixed at 0).
print("\n=== Part D: assoc <=> (C_theta) (closure modulo the syntactic congruence of ker K) ===")
def monoids_unit0(n):
    els = range(n)
    inner = range(1, n)
    for flat in itertools.product(els, repeat=(n - 1) * (n - 1)):
        t = [[0] * n for _ in els]
        for x in els:
            t[0][x] = x; t[x][0] = x
        k = 0
        for a in inner:
            for b in inner:
                t[a][b] = flat[k]; k += 1
        if all(t[t[a][b]][c] == t[a][t[b][c]] for a in inner for b in inner for c in inner):
            yield t
for n in (2, 3, 4):
    nm = nn = mism = assoc_cnt = Cth_cnt = neither_assoc = 0
    for t in monoids_unit0(n):
        nm += 1
        els = range(n)
        for K in itertools.product(els, repeat=n):
            if K[0] != 0 or any(K[K[a]] != K[a] for a in els):
                continue
            if not all(K[t[K[x]][c]] == K[t[x][c]] for x in els for c in els):
                continue
            nn += 1
            img = [a for a in els if K[a] == a]
            theta = lambda a, b: all(K[t[t[x][a]][y]] == K[t[t[x][b]][y]] for x in els for y in els)
            Cth = all(theta(K[t[f][g]], t[f][g]) for f in img for g in img)
            sh = lambda f, g: K[t[f][g]]
            assoc = all(sh(sh(f, g), h) == sh(f, sh(g, h)) for f in img for g in img for h in img)
            mLA = all(K[t[c][K[x]]] == K[t[c][x]] for x in els for c in els)
            closed = all(K[t[f][g]] == t[f][g] for f in img for g in img)
            assoc_cnt += assoc; Cth_cnt += Cth
            if assoc != Cth:
                mism += 1
                if mism <= 3:
                    print("   mismatch:", t, K, "assoc", assoc, "Cth", Cth)
            if assoc and not mLA and not closed:
                neither_assoc += 1
    print(f"order {n}: {nm} monoids (unit 0), {nn} normalizers, assoc {assoc_cnt}, (C_theta) {Cth_cnt}, "
          f"mismatches {mism}; assoc with neither (mLA) nor closure: {neither_assoc}")
