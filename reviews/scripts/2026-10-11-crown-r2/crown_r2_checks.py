"""Second-round referee checks, MacBeth crown Stage 2 (content 4e7264d, PDF from stamp commit cf2c49d).
Rick, 2026-10-11.  Run: python3 crown_r2_checks.py

Conventions (as in the note, Def. 2, Def. 3, Sec. 2.1):
  tree      = parent array, root 0, parent[i] < i (all labelled increasing trees).
  T-morphism = injective node map j, root->root, parent-preserving (so the image is a rooted down-set).
  Cont morphism (S,P) -> (S',P') = (phi: S->S', {rho_s : P'_{phi s} -> P_s}).
  Composite of (phi,rho):A->B then (psi,sig):B->C is (psi.phi, s |-> rho_s . sig_{phi s}).
"""
import itertools
from collections import Counter

def trees(n):
    if n == 1:
        yield (None,); return
    for par in itertools.product(*[range(i) for i in range(1, n)]):
        yield (None,) + par

def anc(par, v):
    out = [v]
    while par[v] is not None:
        v = par[v]; out.append(v)
    return out

def leaves(par):
    kids = {p for p in par if p is not None}
    return [v for v in range(len(par)) if v not in kids]

def branches(par):          # branch = frozenset of its chain, indexed by leaf
    return {l: frozenset(anc(par, l)) for l in leaves(par)}

def embeddings(p, q):
    n, m = len(p), len(q)
    for img in itertools.permutations(range(m), n):
        if img[0] != 0: continue
        if all(q[img[v]] == img[p[v]] for v in range(1, n)):
            yield img

def compose(jf, jg):        # g o f as node maps, f: T->T', g: T'->T''
    return tuple(jg[x] for x in jf)

def subsets(n):
    for bits in range(2 ** n):
        yield frozenset(i for i in range(n) if bits >> i & 1)

def pre(j, Rp):             # j^{-1} R'
    return frozenset(v for v in range(len(j)) if j[v] in Rp)

def img(j, R):
    return frozenset(j[v] for v in R)

OBJ = [t for n in range(1, 5) for t in trees(n)]                     # 10 objects
HOM = {(a, b): list(embeddings(OBJ[a], OBJ[b])) for a in range(len(OBJ)) for b in range(len(OBJ))}
nmor = sum(len(v) for v in HOM.values())
print(f"objects (labelled trees <=4 nodes): {len(OBJ)}, T-morphisms: {nmor}")

# ---------------------------------------------------------------- (1) Psi' is a functor T^op -> Cont
# Psi'(T) = (1,{nodes T});  Psi'(f) = (id_1, rho = j) : Psi'(T') -> Psi'(T),  rho : nodes T -> nodes T'.
def cont_comp(m1, m2):      # m1: A->B, m2: B->C (one-shape: phi trivial); rho_total = rho1 . rho2
    rho1, rho2 = m1, m2
    return tuple(rho2[x] for x in rho1)
bad_id = bad_comp = checked = 0
for a in range(len(OBJ)):
    idj = tuple(range(len(OBJ[a])))
    if HOM[(a, a)].count(idj) != 1: bad_id += 1
    # Psi'(id) has backward map id
for a, b, c in itertools.product(range(len(OBJ)), repeat=3):
    for jf in HOM[(a, b)]:
        for jg in HOM[(b, c)]:
            jgf = compose(jf, jg)
            # Psi'(g o f) : Psi'(T'') -> Psi'(T)  must equal  Psi'(f) o_Cont Psi'(g)
            # Psi'(g): Psi'(T'')->Psi'(T') rho=jg ; Psi'(f): Psi'(T')->Psi'(T) rho=jf
            lhs = jgf
            # careful: in Cont, composite of X:Psi'(T'')->Psi'(T') (rho=jg) then Y:Psi'(T')->Psi'(T) (rho=jf)
            # has rho = rho_X . rho_Y = jg . jf  (as functions nodes T -> nodes T'')
            rhs = tuple(jg[jf[v]] for v in range(len(jf)))
            checked += 1
            if lhs != rhs: bad_comp += 1
print(f"(1) Psi' : T^op -> Cont  identity failures {bad_id}, composable pairs {checked}, composition failures {bad_comp}")
# covariant reading check: would (id, j) be a Cont morphism Psi'(T) -> Psi'(T')? needs rho: nodes T' -> nodes T.
print("    (typing) (id_1,j) needs backward map P_target -> P_source; j: nodes T -> nodes T' => source Psi'(T'), target Psi'(T): contravariant. OK")

# ---------------------------------------------------------------- (2) legs of truncated Fam(cod^op) along (id,rho)
# total-category arrow (A,U) -> (B,V) over (id,rho), rho: P_B -> P_A, exists iff rho(V) subset U.
def arrow(rho, U, V):
    return img(rho, V) <= U
bad_op = bad_cart = bad_match = bad_fun = nlegs = 0
for a, b in itertools.product(range(len(OBJ)), repeat=2):
    nT, nTp = len(OBJ[a]), len(OBJ[b])
    for j in HOM[(a, b)]:
        # Psi'(f): A = Psi'(T') -> B = Psi'(T), rho = j : nodes T -> nodes T'
        for Up in subsets(nTp):                  # predicate over A = Psi'(T')
            nlegs += 1
            # opcartesian pushforward = initial target in fibre over B, fibre order: V -> V' iff V' subset V
            targets = [V for V in subsets(nT) if arrow(j, Up, V)]
            init = [V for V in targets if all(W <= V for W in targets)]
            if init != [pre(j, Up)]: bad_op += 1
            # P-cartesian lift in A at R'=Up along f: brute force universal property inside OBJ
            cands = [R for R in subsets(nT) if img(j, R) <= Up]
            def is_cart(R):
                for s in range(len(OBJ)):
                    for h in HOM[(s, a)]:
                        k = compose(h, j)
                        for Q in subsets(len(OBJ[s])):
                            if img(k, Q) <= Up and not img(h, Q) <= R:
                                return False
                return True
            carts = [R for R in cands if is_cart(R)]
            if carts != [pre(j, Up)]: bad_match += 1
        for V in subsets(nT):                    # cartesian reindexing of V (over B) back to A
            sources = [U for U in subsets(nTp) if arrow(j, U, V)]
            term = [U for U in sources if all(U <= W for W in sources)]   # terminal in (Sub,⊇)
            if term != [img(j, V)]: bad_cart += 1
print(f"(2) {nlegs} (f,R') pairs: opcartesian leg != j^-1 : {bad_op};  P-cartesian lift != j^-1 R' : {bad_match};"
      f"  cartesian leg != direct image : {bad_cart}")
# indexed poset functoriality: (jg jf)^-1 = jf^-1 jg^-1 and monotone
for a, b, c in itertools.product(range(len(OBJ)), repeat=3):
    for jf in HOM[(a, b)]:
        for jg in HOM[(b, c)]:
            jgf = compose(jf, jg)
            for R in subsets(len(OBJ[c])):
                if pre(jgf, R) != pre(jf, pre(jg, R)): bad_fun += 1
print(f"    indexed poset T^op->Pos, f |-> j^-1: strict functoriality failures {bad_fun}")

# ---------------------------------------------------------------- (3) no-go for Psi on morphisms
def search(objs, variance, compat):
    """Count shape-level functors on the full subcategory on objs (branch sets), all phi chosen per morphism
    subject to compat; variance 'co' (phi: br T -> br T') or 'contra' (phi: br T' -> br T)."""
    homs = [(a, b, j) for a in range(len(objs)) for b in range(len(objs)) for j in embeddings(objs[a], objs[b])]
    BR = [branches(t) for t in objs]
    choices = []
    for (a, b, j) in homs:
        src, tgt = (a, b) if variance == 'co' else (b, a)
        dom, cod = list(BR[src]), list(BR[tgt])
        opts = []
        for vals in itertools.product(cod, repeat=len(dom)):
            phi = dict(zip(dom, vals))
            if all(compat(BR, a, b, j, x, phi[x]) for x in dom):
                opts.append(phi)
        choices.append(opts)
    if any(len(o) == 0 for o in choices):
        return 0, [homs[i] for i, o in enumerate(choices) if not o]
    idx = {h: i for i, h in enumerate(homs)}
    count = 0
    for pick in itertools.product(*choices):
        ok = True
        for (a, b, j), pf in zip(homs, pick):
            if a == b and j == tuple(range(len(objs[a]))) and any(pf[x] != x for x in pf): ok = False; break
        if not ok: continue
        for (a, b, jf), pf in zip(homs, pick):
            for (b2, c, jg), pg in zip(homs, pick):
                if b2 != b: continue
                ph = pick[idx[(a, c, compose(jf, jg))]]
                comp = {x: pg[pf[x]] for x in pf} if variance == 'co' else {x: pf[pg[x]] for x in pg}
                if comp != ph: ok = False; break
            if not ok: break
        if ok: count += 1
    return count, None

# compatibilities.  BR[t][leaf] = chain.
def co_strong(BR, a, b, j, x, y):      # x in br(T), y in br(T'): j(chain_x) subset chain_y
    return img(j, BR[a][x]) <= BR[b][y]
def contra_strong(BR, a, b, j, x, y):  # x in br(T'), y in br(T): j(chain_y) subset chain_x   (note, l.419)
    return img(j, BR[a][y]) <= BR[b][x]
def contra_restrict(BR, a, b, j, x, y):  # natural restriction: j^-1(chain_x) subset chain_y
    return pre(j, BR[b][x]) <= BR[a][y]
def free(BR, a, b, j, x, y):
    return True

P2 = (None, 0); V = (None, 0, 0); T4 = (None, 0, 1, 1); STAR3 = (None, 0, 0, 0)
print("(3) no-go searches (shape level; a functor on T restricts to any full subcategory)")
print("   covariant, compatible, {ra, r-a-{c,d}}:", search([P2, T4], 'co', co_strong)[0], "functors")
print("   contravariant, note's compatibility, {ra, V}:", search([P2, V], 'contra', contra_strong))
print("   contravariant, restriction compatibility j^-1(chain b') in chain phi(b'), {ra, V}:",
      search([P2, V], 'contra', contra_restrict)[0], "functors")
print("   contravariant, restriction compatibility, {V, star3}:", search([V, STAR3], 'contra', contra_restrict)[0], "functors")
print("   covariant, NO compatibility, {ra, r-a-{c,d}}:", search([P2, T4], 'co', free)[0], "functors")
print("   contravariant, NO compatibility, {ra, V}:", search([P2, V], 'contra', free)[0], "functors")

# ---------------------------------------------------------------- (4) an unrestricted functor extending Psi exists
# Skeleton: one tree per iso class, <=5 nodes.  Isos = automorphisms -> (id, id).  Non-iso f: T->T' ->
# covariant: phi = const c_{T'} (first leaf), rho_b = const root : chain_{c_{T'}} -> chain_b.
def canon(par):
    n = len(par); kids = {v: [w for w in range(n) if par[w] == v] for v in range(n)}
    def c(v): return "(" + "".join(sorted(c(w) for w in kids[v])) + ")"
    return c(0)
skel = {}
for n in range(1, 6):
    for t in trees(n):
        skel.setdefault(canon(t), t)
S = list(skel.values())
BR = [branches(t) for t in S]
def F(a, b, j, variance):
    if a == b:                                     # skeleton: endo = automorphism
        return ({x: x for x in BR[a]}, {x: {v: v for v in BR[a][x]} for x in BR[a]})
    if variance == 'co':
        c = min(BR[b]); phi = {x: c for x in BR[a]}
        rho = {x: {v: 0 for v in BR[b][c]} for x in BR[a]}     # chain_{phi x} (in T') -> chain_x : root
    else:
        c = min(BR[a]); phi = {x: c for x in BR[b]}
        rho = {x: {v: 0 for v in BR[a][c]} for x in BR[b]}     # chain_{phi x} (in T) -> chain_x (in T') : root
    return phi, rho
def ccomp(m1, m2):          # m1: A->B, m2: B->C in Cont
    (p1, r1), (p2, r2) = m1, m2
    return ({s: p2[p1[s]] for s in p1}, {s: {v: r1[s][r2[p1[s]][v]] for v in r2[p1[s]]} for s in p1})
nS = len(S); H = {(a, b): list(embeddings(S[a], S[b])) for a in range(nS) for b in range(nS)}
for variance in ('co', 'contra'):
    bad = pairs = 0
    for a, b, c in itertools.product(range(nS), repeat=3):
        for jf in H[(a, b)]:
            for jg in H[(b, c)]:
                pairs += 1
                if variance == 'co':
                    lhs = F(a, c, compose(jf, jg), 'co'); rhs = ccomp(F(a, b, jf, 'co'), F(b, c, jg, 'co'))
                else:
                    lhs = F(a, c, compose(jf, jg), 'contra'); rhs = ccomp(F(b, c, jg, 'contra'), F(a, b, jf, 'contra'))
                if lhs != rhs: bad += 1
    print(f"(4) unrestricted {variance}variant functor extending Psi on objects, skeleton <=5 nodes ({nS} objects): "
          f"{pairs} composable pairs, failures {bad}")

# ---------------------------------------------------------------- (5) k=2 on the V-tree
nV = 3; tot = 0; bad = 0
for jT in [t for t in OBJ if len(t) <= 3]:
    for j in embeddings(jT, V):
        j2 = {(u, v): (j[u], j[v]) for u in range(len(jT)) for v in range(len(jT))}
        for bits in range(2 ** (nV * nV)):
            Rp = {(u, v) for u in range(nV) for v in range(nV) if bits >> (3 * u + v) & 1}
            refl = {(u, v) for (u, v) in j2 if (j[u], j[v]) in Rp}
            if {x for x in j2 if j2[x] in Rp} != refl: bad += 1
            tot += 1
print(f"(5) k=2: V-tree has 2^(3^2) = {2**9} binary relations (the note says 2^(n^2) = 16 'on the V-tree');"
      f" (j^2)^-1 = relation reflection on {tot} (embedding, relation) pairs, failures {bad}")

# ---------------------------------------------------------------- (6) the branch relation R_f (note, l.441-443)
# (b,b') in R_f iff j(chain_b) subset chain_b'.  Note claims: left-total, neither left- nor right-functional.
lt = lu = ru_fail = tot = 0
for a, b in itertools.product(range(len(OBJ)), repeat=2):
    Bs, Bt = branches(OBJ[a]), branches(OBJ[b])
    for j in HOM[(a, b)]:
        tot += 1
        R = {(x, y) for x in Bs for y in Bt if img(j, Bs[x]) <= Bt[y]}
        if all(any((x, y) in R for y in Bt) for x in Bs): lt += 1
        if all(sum((x, y) in R for x in Bs) <= 1 for y in Bt): lu += 1      # each b' has <= 1 partner b
        if any(sum((x, y) in R for y in Bt) >= 2 for x in Bs): ru_fail += 1
print(f"(6) R_f over {tot} T-morphisms: left-total in {lt}; each b' related to <=1 b (R_f^op a partial map) in {lu};"
      f" some b with >=2 partners in {ru_fail}")

# ---------------------------------------------------------------- (7) strengthened no-go: only naturality on automorphisms
# Hypothesis: Psi acts on automorphisms by the canonical action on branches; NO condition on other morphisms.
def co_nat_auts(BR, a, b, j, x, y):
    return co_strong(BR, a, b, j, x, y) if a == b else True
def contra_nat_auts(BR, a, b, j, x, y):
    return contra_restrict(BR, a, b, j, x, y) if a == b else True
print("(7) natural on automorphisms only (objects pairwise non-isomorphic, so a==b <=> automorphism):")
print("   covariant, {ra, r-a-{c,d}}:", search([P2, T4], 'co', co_nat_auts)[0], "functors")
print("   contravariant, {V, star3}:", search([V, STAR3], 'contra', contra_nat_auts)[0], "functors")
