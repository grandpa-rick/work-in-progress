"""Kostka-Foulkes K_{nu mu}(t) via Lascoux-Schutzenberger charge; returns coeff lists (ints)."""
import itertools, sys
from functools import lru_cache
sys.path.insert(0, '/home/agent/projects/scripts/day214')
from ds_pointeval import parts, dominates, nstat
def conj(k): return tuple(sum(1 for p in k if p > i) for i in range(k[0])) if k else ()
def charge_std(word):  # word = list of (pos, letter) letters 1..m distinct
    pos = {l: p for p, l in word}; idx = 0; tot = 0
    for r in range(2, len(word)+1):
        if pos[r] > pos[r-1]: idx += 1
        tot += idx
    return tot
def charge(w):
    w = list(enumerate(w)); tot = 0
    while w:
        m = max(l for _, l in w); sub = []; L = len(w)
        # start at right end, scan leftwards cyclically for 1, then 2, ...
        cur = L  # position index in list w (exclusive)
        for r in range(1, m+1):
            j = cur - 1
            for _ in range(L):
                if j < 0: j = L - 1
                if w[j][1] == r and all(w[j] is not s for s in sub): break
                j -= 1
            sub.append(w[j]); cur = j
        tot += charge_std([(p, l) for p, l in sub])
        ids = {id(s) for s in sub}; w = [x for x in w if id(x) not in ids]
    return tot
def ssyt(shape, content):
    """generate SSYT as list of rows"""
    n = sum(shape); letters = []
    for i, c in enumerate(content): letters += [i+1]*c
    res = []
    def rec(T, lets):
        if not lets: res.append([r[:] for r in T]); return
        l = lets[0]
        # place letter l (all copies of l form horizontal strip); place one at a time in increasing row order constraint: simpler brute force
        for i in range(len(shape)):
            if len(T[i]) < shape[i] and (i == 0 or len(T[i-1]) > len(T[i])) and (i == 0 or T[i-1][len(T[i])] < l) and (not T[i] or T[i][-1] <= l):
                # avoid duplicates: copies of same letter placed in nonincreasing row order
                if hasattr(rec, 'last') and False: pass
                T[i].append(l); rec2(T, lets[1:], l, i); T[i].pop()
    def rec2(T, lets, prev, prow):
        if not lets: res.append([r[:] for r in T]); return
        l = lets[0]
        for i in range(len(shape)):
            if l == prev and i > prow: continue
            if len(T[i]) < shape[i] and (i == 0 or len(T[i-1]) > len(T[i])) and (i == 0 or T[i-1][len(T[i])] < l) and (not T[i] or T[i][-1] <= l):
                T[i].append(l); rec2(T, lets[1:], l, i); T[i].pop()
    rec2([[] for _ in shape], letters, None, 99)
    return res
@lru_cache(None)
def KF(nu, mu):
    out = [0]*(nstat(mu)+1)
    for T in ssyt(nu, mu):
        w = [x for row in reversed(T) for x in row]; out[charge(w)] += 1
    while len(out) > 1 and out[-1] == 0: out.pop()
    return tuple(out)
if __name__ == '__main__':
    print(KF((2,1),(1,1,1)), KF((3,),(1,1,1)), KF((2,2),(2,1,1)), KF((3,1),(2,1,1)), KF((4,2),(2,2,1,1)))
    # sanity: K(1) = Kostka, sum_nu f^nu K_{nu,1^n}(t) = [n]!
