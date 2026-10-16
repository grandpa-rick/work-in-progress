# Lattice (N^I)-graded braided ABS with twisted target A=k_phi[N^I], u_a u_b = q^{phi(a,b)} u_{a+b}
from fractions import Fraction as Fr
from collections import defaultdict
from itertools import product as iprod
import random
import sympy as sp

def add(a, b): return tuple(x+y for x, y in zip(a, b))
def form(a, b, C): return sum(a[i]*C[i][j]*b[j] for i in range(len(a)) for j in range(len(b)))

class QQ:
    def __init__(s, C, B, q): s.C, s.B, s.q = C, B, q
    def br(s, a, b): return s.q ** form(a, b, s.C)
    def ph(s, a, b): return s.q ** form(a, b, s.B)
    def qqsh(s, a, b):
        memo = {}
        def tot(w):
            r = tuple(0 for _ in s.C)
            for x in w: r = add(r, x)
            return r
        def rec(i, j):
            if (i, j) in memo: return memo[(i, j)]
            out = defaultdict(Fr)
            if i == len(a) and j == len(b): out[()] = Fr(1)
            else:
                if i < len(a):
                    for w, c in rec(i+1, j).items(): out[(a[i],)+w] += c
                if j < len(b):
                    wt = s.br(b[j], tot(a[i:]))
                    for w, c in rec(i, j+1).items(): out[(b[j],)+w] += c*wt
                if i < len(a) and j < len(b):
                    wt = s.br(b[j], tot(a[i+1:])) * s.ph(a[i], b[j])
                    for w, c in rec(i+1, j+1).items(): out[(add(a[i], b[j]),)+w] += c*wt
            memo[(i, j)] = out; return out
        return rec(0, 0)
    def mul(s, X, Y):
        out = defaultdict(Fr)
        for a, ca in X.items():
            for b, cb in Y.items():
                for w, c in s.qqsh(a, b).items(): out[w] += ca*cb*c
        return {k: v for k, v in out.items() if v != 0}
    def Psi(s, word, t):
        """word: tuple of simple-root indices (free braided algebra, primitive F_i of degree e_i);
        zeta(F_i)=t_i u_{e_i}; returns coefficients on M_alpha (alpha = tuple of lattice vectors)"""
        n = len(word); I = len(s.C)
        e = [tuple(1 if k == i else 0 for k in range(I)) for i in range(I)]
        if n == 0: return {(): Fr(1)}
        out = defaultdict(Fr)
        for l in range(1, n+1):
            for b in iprod(range(l), repeat=n):
                if len(set(b)) < l: continue
                wt = Fr(1)
                for i in range(n):
                    for j in range(i+1, n):
                        if b[i] > b[j]: wt *= s.br(e[word[i]], e[word[j]])
                # within each block, zeta(subword) = prod t * q^{sum_{k<l in block} phi(e_k,e_l)} u_{deg}
                al = []
                for k in range(l):
                    pos = [i for i in range(n) if b[i] == k]
                    d = tuple(0 for _ in range(I)); c = Fr(1)
                    for x, i in enumerate(pos):
                        for i2 in pos[x+1:]: c *= s.ph(e[word[i]], e[word[i2]])
                        c *= t[word[i]]; d = add(d, e[word[i]])
                    al.append(d); wt *= c
                out[tuple(al)] += wt
        return {k: v for k, v in out.items() if v != 0}
    def PsiLin(s, X, t):
        out = defaultdict(Fr)
        for w, c in X.items():
            for k, v in s.Psi(w, t).items(): out[k] += c*v
        return {k: v for k, v in out.items() if v != 0}
