"""Enumerate all skew braces (G,+,o) on a given finite abelian additive group G,
via lambda: G -> Aut(G,+) with lambda_{a + lambda_a(b)} = lambda_a lambda_b.
a o b = a + lambda_a(b).  Backtracking with propagation.
Elements of G = Z/n1 x ... x Z/nk encoded as tuples -> ints.
"""
import itertools


class AbGroup:
    def __init__(self, mods):
        self.mods = mods
        self.elts = list(itertools.product(*[range(m) for m in mods]))
        self.idx = {e: i for i, e in enumerate(self.elts)}
        self.n = len(self.elts)
        n = self.n
        self.add = [[self.idx[tuple((x + y) % m for x, y, m in zip(a, b, mods))]
                     for b in self.elts] for a in self.elts]
        self.neg = [self.idx[tuple((-x) % m for x, m in zip(a, mods))] for a in self.elts]
        self.zero = self.idx[tuple(0 for _ in mods)]
        self.auts = self._auts()

    def _auts(self):
        # automorphisms determined by images of standard generators e_i
        k = len(self.mods)
        gens = []
        for i in range(k):
            e = [0] * k
            e[i] = 1
            gens.append(self.idx[tuple(e)])
        order = [self._order(g) for g in range(self.n)]
        auts = []
        cands = [[g for g in range(self.n) if self.mods[i] % order[g] == 0] for i in range(k)]
        for imgs in itertools.product(*cands):
            perm = [None] * self.n
            ok = True
            for i_, e in enumerate(self.elts):
                v = self.zero
                for c, g in zip(e, imgs):
                    for _ in range(c):
                        v = self.add[v][g]
                perm[i_] = v
            if len(set(perm)) == self.n:
                auts.append(tuple(perm))
        return auts

    def _order(self, g):
        o, v = 1, g
        while v != self.zero:
            v = self.add[v][g]
            o += 1
        return o


def enumerate_braces(G, limit=None):
    n = G.n
    A = G.auts
    aidx = {a: i for i, a in enumerate(A)}
    ident = aidx[tuple(range(n))]
    comp = [[aidx[tuple(A[x][A[y][t]] for t in range(n))] for y in range(len(A))]
            for x in range(len(A))]
    out = []

    def propagate(lam):
        changed = True
        while changed:
            changed = False
            assigned = [a for a in range(n) if lam[a] is not None]
            for a in assigned:
                La = A[lam[a]]
                for b in assigned:
                    c = G.add[a][La[b]]
                    v = comp[lam[a]][lam[b]]
                    if lam[c] is None:
                        lam[c] = v
                        changed = True
                        assigned.append(c)
                    elif lam[c] != v:
                        return False
        return True

    def rec(lam):
        if limit and len(out) >= limit:
            return
        try:
            a = lam.index(None)
        except ValueError:
            out.append(tuple(lam))
            return
        for v in range(len(A)):
            l2 = list(lam)
            l2[a] = v
            if propagate(l2):
                rec(l2)

    lam = [None] * n
    lam[G.zero] = ident
    rec(lam)
    return out, A
