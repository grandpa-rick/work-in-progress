# Wake 252 Case-B gate, fast version for m=6 (and m<=5 regression):
#   H_kappa  from LS-parity model  H(j) = F_kappa(2j)  (Wake 251, verified m<=5), interpolated through j=0..|kappa|
#   [O_D]    from Buch-Ravikumar K-Pieri (brpieri.py; reproduces kclass.py 46/46 on m=4,5)
# Consistency: b(j) := H(j)-H(j-1) == sum_rho c_rho H_rho(j) on a half-integer grid (joint check of LS model + BR).
import sys, itertools, collections
from fractions import Fraction as Fr
from functools import lru_cache
sys.path.insert(0, '/home/agent/projects/work-in-progress/scripts/wake251/orderpoly')
from orderpoly import shifted, boxless, longest_chain, ideals
import brpieri
N = int(sys.argv[1]) if len(sys.argv) > 1 else 6
def countF(lam, Nmax):   # F(0..Nmax): order-preserving g: P->{0..N}, even on diagonal; DP over ideals, all N at once
    P = shifted(list(lam)); J = ideals(P, boxless); full = frozenset(P)
    sub = {I: [K for K in J if I <= K] for I in J}
    vec = {I: 1 for I in J}   # I = {g<=0}; 0 is even so diagonal boxes allowed
    out = [vec.get(full, 0)]
    for v in range(1, Nmax + 1):
        nv = collections.defaultdict(int)
        for I, c in vec.items():
            for K in sub[I]:
                if v % 2 == 1 and any(b[0] == b[1] for b in K - I): continue
                nv[K] += c
        vec = nv; out.append(vec.get(full, 0))
    return out
@lru_cache(None)
def Hvals(k):
    if not k: return (1,)
    F = countF(k, 2 * sum(k)); return tuple(F[2 * j] for j in range(sum(k) + 1))
def H(k, x):
    v = Hvals(k); n = len(v); s = Fr(0)
    for i in range(n):
        t = Fr(v[i])
        for l in range(n):
            if l != i: t *= Fr(x - l, i - l)
        s += t
    return s
def r(k): return longest_chain(shifted(list(k)), boxless) if k else 0
def rp(k): return r(k[1:])
def caseB(k): return any(k[i] <= r(k[i + 1:]) + 1 for i in range(len(k) - 1))
def nm(t): return ''.join(map(str, t)) or '0'
if __name__ == '__main__':
    Ms = list(range(1, 4 * N + 4, 2)); rows = []
    for L in range(1, N + 1):
      for c in itertools.combinations(range(1, N + 1), L):
        k = tuple(sorted(c, reverse=True)); K = brpieri.pieri_dim(k, N)
        pts = [Fr(t, 2) for t in range(-4 * N - 4, 2 * sum(k) + 4)]
        b = lambda x: H(k, x) - H(k, x - 1)
        kok = all(b(x) == sum(cc * H(rr, x) for rr, cc in K.items()) for x in pts)
        a = (rp(k) + 1) // 2; top = 2 * a + 1; need = list(range(3, top, 2))
        zH = [M for M in Ms if H(k, Fr(-M, 2)) == 0]; zb = [M for M in Ms if b(Fr(-M, 2)) == 0]
        prop = all(M in zb for M in need); law = zH == list(range(3, top + 1, 2))
        bad = {nm(rr): [M for M in need if H(rr, Fr(-M, 2)) != 0] for rr in K}; bad = {q: v for q, v in bad.items() if v}
        ind = all(rp(rr) >= rp(k) - 2 for rr in K)
        rows.append(dict(k=k, B=caseB(k), need=need, prop=prop, law=law, kok=kok, tbt=not bad, ind=ind, K=K))
        print(k, 'B' if caseB(k) else 'A', "r'", rp(k), 'zH', zH, 'zb', zb, 'need', need, 'prop', prop, 'law', law,
              'Kcls', kok, 'Ind', ind, 'K', {nm(q): v for q, v in K.items()}, 'nonvanishing terms', bad, flush=True)
    print('BR K-class x LS-model reproduces b:', sum(R['kok'] for R in rows), '/', len(rows))
    for lab, f in (('A', lambda R: not R['B']), ('B', lambda R: R['B'])):
        S = [R for R in rows if f(R)]; nt = [R for R in S if R['need']]
        print(f"Case {lab}: {len(S)} kappa, {len(nt)} with nonempty need; propagation {sum(R['prop'] for R in S)}/{len(S)}; "
              f"law {sum(R['law'] for R in S)}/{len(S)}; every K-term vanishes {sum(R['tbt'] for R in nt)}/{len(nt)}; "
              f"(Ind) {sum(R['ind'] for R in S)}/{len(S)}")
