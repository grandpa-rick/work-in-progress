"""Referee checks for MacBeth crown Stage 2 (4136ecf). Rick, 2026-10-10."""
import itertools
from collections import Counter

def trees(n):
    """All rooted trees on nodes 0..n-1 as parent arrays (parent[i]<i), root 0 (labelled, increasing)."""
    if n == 1:
        yield (None,)
        return
    for par in itertools.product(*[range(i) for i in range(1, n)]):
        yield (None,) + par

def ancestors(par, v):
    out = [v]
    while par[v] is not None:
        v = par[v]; out.append(v)
    return out

def leaves(par):
    n = len(par); kids = {p for p in par if p is not None}
    return [v for v in range(n) if v not in kids]

def lam(par):
    c = Counter()
    for l in leaves(par):
        for v in ancestors(par, l): c[v] += 1
    return c

def embeddings(p, q):
    """rooted subtree embeddings T_p -> T_q: injective, root->root, parent-preserving (=> prefix-closed image)."""
    n, m = len(p), len(q)
    for img in itertools.permutations(range(m), n):
        if img[0] != 0: continue
        if all(q[img[v]] == img[p[v]] for v in range(1, n)):
            yield img

# ---- (A) over-count factor 2^{sum(lambda-1)} : container fibre over Psi(T) vs arboreal fibre
print("(A) over-count")
for n in range(1, 7):
    bad = 0; seen = 0
    for par in trees(n):
        L = lam(par)
        U = [(l, v) for l in leaves(par) for v in ancestors(par, l)]
        assert len(U) == sum(L.values())
        # prefix-consistent colourings of U
        cons = 0
        for bits in range(2 ** len(U)):
            A = {U[i] for i in range(len(U)) if bits >> i & 1}
            if all(((l, v) in A) == ((l2, v) in A) for (l, v) in U for (l2, v2) in U if v2 == v):
                cons += 1
        seen += 1
        if not (cons == 2 ** n and 2 ** len(U) == 2 ** n * 2 ** sum(L[v] - 1 for v in range(n))):
            bad += 1
    print(f"  n={n}: {seen} labelled trees, failures {bad}")
V = (None, 0, 0); Y = (None, 0, 1, 1)
for name, t in [("V", V), ("Y", Y)]:
    L = lam(t); U = sum(L.values())
    print(f"  {name}: |nodes|={len(t)} |U|={U}  2^n={2**len(t)} 2^|U|={2**U}")

# ---- (B) (*) by brute force in MacBeth's A: objects (T,R), morphisms = rooted-subtree embeddings j with R subset j^{-1}R'
print("(B) cartesian <=> R = j^{-1}R' <=> pathwise embedding  (trees up to N nodes)")
N = 4
objs = []
shapes = []
for n in range(1, N + 1):
    for par in trees(n):
        # dedupe up to iso crudely: keep all labelled (harmless for the check)
        shapes.append(par)
for par in shapes:
    for bits in range(2 ** len(par)):
        objs.append((par, frozenset(v for v in range(len(par)) if bits >> v & 1)))
emb = {}
for p in shapes:
    for q in shapes:
        emb[(p, q)] = list(embeddings(p, q))
def hom(X, Y_):
    (p, R), (q, R2) = X, Y_
    return [j for j in emb[(p, q)] if all(j[v] in R2 for v in R)]
def is_cartesian(X, Y_, j):
    p, R = X; q, R2 = Y_
    for Z in objs:
        s, Q = Z
        for k in hom(Z, Y_):
            for h in emb[(s, p)]:
                if tuple(j[h[v]] for v in range(len(s))) != k: continue
                if not all(h[v] in R for v in Q):  # lift over h must exist (uniqueness automatic: P faithful)
                    return False
    return True
def chains(par):
    return [ancestors(par, v) for v in range(len(par))]
def pathwise(X, Y_, j):
    p, R = X; q, R2 = Y_
    return all(all((v in R) == (j[v] in R2) for v in c) for c in chains(p))
import random
random.seed(1)
cnt = Counter()
pairs = [(X, Y_) for X in objs for Y_ in objs if len(X[0]) <= len(Y_[0])]
random.shuffle(pairs)
for X, Y_ in pairs:  # exhaustive
    for j in hom(X, Y_):
        R2 = Y_[1]; pre = frozenset(v for v in range(len(X[0])) if j[v] in R2)
        c = is_cartesian(X, Y_, j); e = (X[1] == pre); pw = pathwise(X, Y_, j)
        cnt[(c, e, pw)] += 1
print("  (cartesian, R==j^-1R', pathwise) counts:", dict(cnt))

# ---- (C) Definition 6 legs vs relation reflection on the V-tree
print("(C) Def 6 nearest-ancestor retraction rho_f vs j^{-1}")
T  = (None, 0)          # r, a
Tp = (None, 0, 0)       # r, a, b
j = (0, 1)
def rho(v):  # nearest ancestor of v (in T') lying in image of j
    img = set(j)
    while v not in img: v = Tp[v]
    return j.index(v)
for bits in range(8):
    R2 = {v for v in range(3) if bits >> v & 1}
    jinv = {v for v in range(2) if j[v] in R2}
    direct = {rho(v) for v in R2}
    if jinv != direct:
        print(f"  R'={sorted(R2)}: j^-1 R'={sorted(jinv)}  Sigma_rho R'={sorted(direct)}  (differ)")
print("  rho^* goes Sub(T)->Sub(T'): wrong direction for reflection; e.g. rho^{-1}({r}) =",
      sorted(v for v in range(3) if rho(v) == 0))

# ---- (D) Psi on morphisms: branches of T' restricted to T need not be branches of T
print("(D) Psi on the V-tree inclusion T={r,a} -> T'={r,a,b}: branch rb restricts to",
      [v for v in ancestors(Tp, 2) if v in set(j)], "(not a branch of T: leaves(T) =", leaves(T), ")")
