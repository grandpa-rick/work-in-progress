# Brute force: braces on A = Z2 x Z4, ideals D ~ Z/4 that are trivial braces,
# (mult-split, add-split, sb-split) over H = G/D ~ Z2, by exhausting all
# normalised sections s: H -> G (i.e. all theta), for MacBeth UID 292 W4 claim.
import itertools
A = [(x, y) for x in range(2) for y in range(4)]
idx = {a: i for i, a in enumerate(A)}
def add(a, b): return ((a[0]+b[0]) % 2, (a[1]+b[1]) % 4)
def neg(a): return ((-a[0]) % 2, (-a[1]) % 4)
# automorphisms as tuples of images of A
auts = []
for i1, i2 in itertools.product(A, A):
    f = {}
    for (x, y) in A:
        v = (0, 0)
        for _ in range(x): v = add(v, i1)
        for _ in range(y): v = add(v, i2)
        f[(x, y)] = v
    # well-defined hom: 2*i1 = 0
    if add(i1, i1) != (0, 0): continue
    if len(set(f.values())) == 8:
        auts.append(tuple(f[a] for a in A))
auts = sorted(set(auts)); assert len(auts) == 8
def ap(f, a): return f[idx[a]]
def comp(f, g): return tuple(f[idx[g[i]]] for i in range(8))
ID = tuple(A)
# DFS for lambda: A -> Aut with lambda_{a+lambda_a(b)} = lambda_a lambda_b
sols = []
def propagate(lam):
    changed = True
    while changed:
        changed = False
        for a in A:
            if lam[a] is None: continue
            for b in A:
                if lam[b] is None: continue
                c = add(a, ap(lam[a], b)); v = comp(lam[a], lam[b])
                if lam[c] is None: lam[c] = v; changed = True
                elif lam[c] != v: return False
    return True
def dfs(lam):
    if not propagate(lam): return
    free = [a for a in A if lam[a] is None]
    if not free:
        sols.append(dict(lam)); return
    for f in auts:
        l2 = dict(lam); l2[free[0]] = f; dfs(l2)
lam0 = {a: None for a in A}; lam0[(0, 0)] = ID
dfs(lam0)
sols = [dict(t) for t in {tuple(sorted(s.items())) for s in sols}]
print("braces (labelled) on Z2xZ4:", len(sols))
Z4subs = []
for g in A:
    S = {(0, 0)}; v = g
    while v != (0, 0): S.add(v); v = add(v, g)
    if len(S) == 4 and frozenset(S) not in Z4subs and all(
        (lambda o: o)(1) for _ in [0]): 
        # cyclic?
        Z4subs.append(frozenset(S))
Z4subs = list(set(Z4subs))
print("cyclic order-4 subgroups:", [sorted(s) for s in Z4subs])
def circ(l, a, b): return add(a, ap(l[a], b))
count = {}
witness = None
for l in sols:
    for D in Z4subs:
        # ideal: lambda_a(D) = D, D circ-normal
        if not all(ap(l[a], d) in D for a in A for d in D): continue
        cinv = {a: next(b for b in A if circ(l, a, b) == (0, 0)) for a in A}
        if not all(circ(l, circ(l, a, d), cinv[a]) in D for a in A for d in D): continue
        trivial = all(ap(l[d], e) == e for d in D for e in D)
        if not trivial: continue
        # sections s: H={0,1} -> G, s(0)=0, s(1) in coset outside D
        outs = [k for k in A if k not in D]
        beta0 = [k for k in outs if add(k, k) == (0, 0)]          # additive hom section
        tau0 = [k for k in outs if circ(l, k, k) == (0, 0)]       # multiplicative hom section
        both = [k for k in outs if k in beta0 and k in tau0]
        key = (bool(tau0), bool(beta0), bool(both))
        insoc = all(ap(l[d], a) == a for d in D for a in A)
        count[(key, insoc)] = count.get((key, insoc), 0) + 1
        if key == (True, True, False) and witness is None:
            witness = (l, D, beta0, tau0)
print("(mult,add,sb) , D<=Soc : count of labelled (brace, D) pairs")
for k, v in sorted(count.items()): print(k, v)
if witness:
    l, D, b0, t0 = witness
    print("witness D =", sorted(D), "add-hom sections s(1) in", b0, "mult-hom s(1) in", t0)
    for a in A: print(a, "lambda:", [ap(l[a], e) for e in [(1, 0), (0, 1)]])

# ---- cross-check 1: naive enumeration of lambda maps (completeness of DFS)
naive = 0
others = [a for a in A if a != (0, 0)]
for choice in itertools.product(range(8), repeat=7):
    lam = {(0, 0): ID}
    for a, c in zip(others, choice): lam[a] = auts[c]
    ok = True
    for a in A:
        la = lam[a]
        for b in A:
            if lam[add(a, ap(la, b))] != comp(la, lam[b]): ok = False; break
        if not ok: break
    if ok: naive += 1
print("naive count:", naive)
# ---- cross-check 2: per-D breakdown for D = {0} x Z4
D0 = frozenset((0, y) for y in range(4))
for l in sols:
    if not all(ap(l[a], d) in D0 for a in A for d in D0): continue
    cinv = {a: next(b for b in A if circ(l, a, b) == (0, 0)) for a in A}
    if not all(circ(l, circ(l, a, d), cinv[a]) in D0 for a in A for d in D0): continue
    if not all(ap(l[d], e) == e for d in D0 for e in D0): continue
    outs = [k for k in A if k not in D0]
    b0 = [k for k in outs if add(k, k) == (0, 0)]
    t0 = [k for k in outs if circ(l, k, k) == (0, 0)]
    print("D={0}xZ4: add-hom s(1):", b0, " mult-hom s(1):", t0,
          " lambda_(1,0),(0,1) on gens:", [[ap(l[g], e) for e in [(1,0),(0,1)]] for g in [(1,0),(0,1)]])
