# Day 234 shared helpers for the printed-statement sweep. Everything here is re-implemented from the
# PRINTED definitions in fpsac2027-draft.pdf (see printed.txt); ground truth comes from engine.py.
import sys, pickle, os, itertools
from functools import lru_cache
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from engine import *
import sympy as sp

ST, TT = sp.symbols('s t')

def load_table():
    d = pickle.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'ctable_n7.pkl'), 'rb'))
    return {lam: {mu: R(sp.sympify(c)) for mu, c in row.items()} for lam, row in d.items()}

def c(T, lam, mu):
    return T[lam].get(mu, R(0))

def nfun(lam): return sum(i * l for i, l in enumerate(lam))
def conj(lam): return tuple(sum(1 for l in lam if l > i) for i in range(lam[0])) if lam else ()
def srt(p): return tuple(sorted([x for x in p if x > 0], reverse=True))

def dom(mu, lam):  # mu ⊵ lam
    if sum(mu) != sum(lam): return False
    a = b = 0
    for i in range(max(len(mu), len(lam))):
        a += mu[i] if i < len(mu) else 0
        b += lam[i] if i < len(lam) else 0
        if a < b: return False
    return True

NEG = float('-inf')

@lru_cache(None)
def kappa(lam, mu):
    """Printed Definition (Sec. 3): largest r with lam = ⊔ lam^i, mu = ⊔ mu^i (multisets), blocks nonempty,
    mu^i ⊵ lam^i (so equal sizes); -inf if mu not ⊵ lam."""
    lam = srt(lam); mu = srt(mu)
    if not lam and not mu: return 0
    if not lam or not mu: return NEG
    best = NEG
    first, rest = lam[0], lam[1:]
    seen = set()
    for r in range(len(rest) + 1):
        for S in itertools.combinations(range(len(rest)), r):
            blk = srt((first,) + tuple(rest[i] for i in S))
            rem = tuple(rest[i] for i in range(len(rest)) if i not in S)
            if (blk, rem) in seen: continue
            seen.add((blk, rem))
            sz = sum(blk)
            for r2 in range(1, len(mu) + 1):
                for U in itertools.combinations(range(len(mu)), r2):
                    mb = srt(tuple(mu[i] for i in U))
                    if sum(mb) != sz or not dom(mb, blk): continue
                    mrem = tuple(mu[i] for i in range(len(mu)) if i not in U)
                    k2 = kappa(srt(rem), srt(mrem))
                    if k2 != NEG: best = max(best, 1 + k2)
    return best

def set_partitions(lst):
    if not lst:
        yield []; return
    first, rest = lst[0], lst[1:]
    for p in set_partitions(rest):
        for i in range(len(p)):
            yield p[:i] + [[first] + p[i]] + p[i + 1:]
        yield [[first]] + p

def qint(m, q=None):
    q = t if q is None else q
    return sum((q**i for i in range(m)), R(0))

def is_coarsening(lam, mu): return kappa(lam, mu) == len(mu)

class Tally:
    def __init__(self, name): self.name = name; self.ok = 0; self.bad = 0; self.msgs = []
    def add(self, cond, msg=''):
        if cond: self.ok += 1
        else:
            self.bad += 1
            if len(self.msgs) < 5: self.msgs.append(msg)
    def line(self):
        return f"{self.name}: {self.ok} pass / {self.bad} fail" + ("  e.g. " + " | ".join(self.msgs) if self.msgs else "")

def report(main, neg):
    v = 'PASS' if main.bad == 0 and main.ok > 0 else 'FAIL'
    nv = 'neg-control FIRES' if neg.bad > 0 else 'neg-control DID NOT FIRE'
    print(f"RESULT [{main.name}] {v}: cases={main.ok + main.bad} fails={main.bad} | {nv} ({neg.bad}/{neg.ok + neg.bad} perturbed cases rejected)")
    for m in main.msgs: print('   ', m)
