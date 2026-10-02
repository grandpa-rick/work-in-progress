"""Day 217: independent enumeration of skew braces with abelian additive group A
(regular subgroups of Hol(A)), for re-deriving MacBeth's o_nu lock / torsor claims.
Written from scratch (no MacBeth code)."""
import itertools, sys
from functools import lru_cache


class AbGroup:
    def __init__(self, mods):
        self.mods = tuple(mods)
        self.elts = list(itertools.product(*[range(m) for m in mods]))
        self.idx = {e: i for i, e in enumerate(self.elts)}
        self.n = len(self.elts)
        self.add = [[self.idx[tuple((a + b) % m for a, b, m in zip(x, y, self.mods))]
                     for y in self.elts] for x in self.elts]
        self.zero = self.idx[tuple(0 for _ in mods)]
        self.neg = [next(j for j in range(self.n) if self.add[i][j] == self.zero) for i in range(self.n)]

    def order(self, i):
        k, x = 1, i
        while x != self.zero:
            x = self.add[x][i]; k += 1
        return k

    def mul(self, k, i):
        x = self.zero
        for _ in range(k % self.n if k >= 0 else 0):
            x = self.add[x][i]
        return x

    def automorphisms(self):
        gens = []
        for j, m in enumerate(self.mods):
            e = [0] * len(self.mods); e[j] = 1; gens.append(self.idx[tuple(e)])
        auts = []
        cands = [[x for x in range(self.n) if m % self.order(x) == 0] for m in self.mods]
        for imgs in itertools.product(*cands):
            f = []
            for e in self.elts:
                y = self.zero
                for c, im in zip(e, imgs):
                    y = self.add[y][self.mul(c, im)]
                f.append(y)
            if len(set(f)) == self.n:
                auts.append(tuple(f))
        return auts


def regular_subgroups(A):
    """Return list of lam: tuple indexed by a of automorphism tuples, one per regular subgroup."""
    auts = A.automorphisms()
    comp = lambda f, g: tuple(f[g[x]] for x in range(A.n))
    ident = tuple(range(A.n))

    def hmul(p, q):  # (a,f)(b,g) = (a+f(b), fg)
        return (A.add[p[0]][p[1][q[0]]], comp(p[1], q[1]))

    def closure(gens):
        S = {(A.zero, ident)}
        frontier = list(S)
        firsts = {A.zero: ident}
        while frontier:
            new = []
            for x in frontier:
                for g in gens:
                    y = hmul(x, g)
                    if y not in S:
                        if y[0] in firsts and firsts[y[0]] != y[1]:
                            return None
                        firsts[y[0]] = y[1]
                        S.add(y); new.append(y)
                        if len(S) > A.n:
                            return None
            frontier = new
        return firsts

    found = set()
    results = []

    def rec(gens, firsts):
        if len(firsts) == A.n:
            key = tuple(firsts[a] for a in range(A.n))
            if key not in found:
                found.add(key); results.append(key)
            return
        a = min(x for x in range(A.n) if x not in firsts)
        for f in auts:
            g2 = gens + [(a, f)]
            fs = closure(g2)
            if fs is not None:
                rec(g2, fs)

    rec([], {A.zero: ident})
    return results


def brace_ops(A, lam):
    circ = [[A.add[a][lam[a][b]] for b in range(A.n)] for a in range(A.n)]
    zero = A.zero
    inv = [next(b for b in range(A.n) if circ[a][b] == zero) for a in range(A.n)]
    return circ, inv


def check_brace(A, lam, circ):
    for a in range(A.n):
        for b in range(A.n):
            ab = circ[a][b]
            for c in range(A.n):
                # a o (b+c) = a o b - a + a o c
                if circ[a][A.add[b][c]] != A.add[A.add[ab][A.neg[a]]][circ[a][c]]:
                    return False
    return True
