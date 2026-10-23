# Wake 252 Case-B gate: hyperplane-section boundary b_k(j)=H_k(j)-H_k(j-1) vs inclusion-exclusion over
# the codim-1 Schubert divisors X_nu (nu = kappa minus one box, strict).  Bruhat order on LG(n,2n) Schubert
# varieties = containment of strict partitions; it is a distributive lattice (meet = componentwise min), and
# X_a ∩ X_b = X_{a∧b} (reduced, Ramanathan).  Two independent expansions of chi(O_∂(j)):
#   (IE)  sum over nonempty sets S of covers-below of (-1)^{|S|+1} H_{∧S}(j)          (Mayer-Vietoris)
#   (MOB) sum_{rho<kappa} -mu(rho,kappa) H_rho(j)                                       (Knutson-type union formula)
import sys, itertools
from fractions import Fraction as Fr
from functools import lru_cache
sys.path.insert(0, '/home/agent/projects/work-in-progress/scripts/prove249d')
sys.path.insert(0, '/home/agent/projects/work-in-progress/scripts/wake251/orderpoly')
from dem import chi
from orderpoly import shifted, boxless, longest_chain
N = int(sys.argv[1]) if len(sys.argv) > 1 else 5
@lru_cache(None)
def Hvals(k):  # k: decreasing tuple; values H(-|k|..0) (negative weights: small Demazure characters, fast)
    if not k: return (1,)
    return tuple(chi(set(k), N, j)[0] for j in range(-sum(k), 1))
def H(k, x):  # Lagrange evaluation at rational x
    v = Hvals(k); n = len(v); s = Fr(0); o = n - 1   # node i <-> j = i - o
    for i in range(n):
        t = Fr(v[i])
        for l in range(n):
            if l != i: t *= Fr(x - (l - o), i - l)
        s += t
    return s
def r(k): return longest_chain(shifted(list(k)), boxless) if k else 0
def rp(k): return r(k[1:])
def covers_below(k):
    out = []
    for i in range(len(k)):
        nu = list(k); nu[i] -= 1
        nu = tuple(x for x in nu if x > 0)
        if all(nu[t] > nu[t + 1] for t in range(len(nu) - 1)) and len(set(nu)) == len(nu): out.append(nu)
    return out
def leq(a, b): return len(a) <= len(b) and all(a[i] <= b[i] for i in range(len(a)))
def meet(a, b): return tuple(x for x in (min(p, q) for p, q in zip(a, b)) if x > 0)
def below(k):
    return [tuple(sorted(c, reverse=True)) for L in range(len(k) + 1) for c in itertools.combinations(range(1, N + 1), L)
            if leq(tuple(sorted(c, reverse=True)), k)]
def mobius_coeffs(k):  # c_rho = -mu(rho,k), rho<k
    els = below(k); mu = {k: 1}
    for rho in sorted(els, key=sum, reverse=True):
        if rho == k: continue
        mu[rho] = -sum(mu[s] for s in mu if s != rho and leq(rho, s))
    return {rho: -mu[rho] for rho in els if rho != k and mu[rho]}
def ie_terms(k):
    cv = covers_below(k); terms = []
    for L in range(1, len(cv) + 1):
        for S in itertools.combinations(cv, L):
            m = S[0]
            for s in S[1:]: m = meet(m, s)
            terms.append(((-1) ** (L + 1), m))
    return terms
def caseB(k):  # proof file §6quater remark: some kappa_i <= r(kappa_{>i})+1
    return any(k[i] <= r(k[i + 1:]) + 1 for i in range(len(k) - 1))
def caseB1(k): return len(k) >= 2 and k[0] <= rp(k) + 1   # kappa_1-only variant
def main():
    Ms = range(1, 4 * N + 4, 2)
    stats = dict(xcheck=0, tot=0)
    rows = []
    for L in range(1, N + 1):
      for c in itertools.combinations(range(1, N + 1), L):
        k = tuple(sorted(c, reverse=True)); stats['tot'] += 1
        mob = mobius_coeffs(k); ie = ie_terms(k)
        # cross-check at integers and half-integers
        pts = [Fr(t, 2) for t in range(-4 * N - 4, 2 * sum(k) + 4)]
        ok = all(H(k, x) - H(k, x - 1) == sum(c * H(rho, x) for rho, c in mob.items()) == sum(s * H(m, x) for s, m in ie) for x in pts)
        stats['xcheck'] += ok
        a = (rp(k) + 1) // 2; top = 2 * a + 1
        zH = [M for M in Ms if H(k, Fr(-M, 2)) == 0]
        zb = [M for M in Ms if H(k, Fr(-M, 2)) - H(k, Fr(-M, 2) - 1) == 0]
        need = list(range(3, top, 2))   # b(-M/2)=0 for M=3..top-2  (i.e. M+2 in [5,top])
        prop = all(M in zb for M in need)
        lawH = zH == list(range(3, top + 1, 2))
        thmC = (H(k, Fr(-3, 2)) == 0) == (len(k) >= 2)
        # term-by-term: every Mobius-support term vanishes at each needed M
        mob_tbt = all(H(rho, Fr(-M, 2)) == 0 for rho in mob for M in need)
        ie_tbt = all(H(m, Fr(-M, 2)) == 0 for _, m in ie for M in need)
        badM = sorted({rho for rho in mob for M in need if H(rho, Fr(-M, 2)) != 0})
        rows.append((k, caseB(k), caseB1(k), rp(k), zH, zb, need, prop, lawH, thmC, mob_tbt, ie_tbt, mob, badM))
        print(k, 'B' if caseB(k) else 'A', "r'=", rp(k), 'zH', zH, 'zb', zb, 'need b at', need, 'prop', prop, 'law', lawH,
              'thmC', thmC, 'xcheck', ok, 'MOBtbt', mob_tbt, 'IEtbt', ie_tbt, 'mob', {''.join(map(str, p)) or '0': v for p, v in mob.items()},
              'nonvanishing mob terms', badM, flush=True)
    print('cross-check', stats)
    for lab, f in (('A', lambda R: not R[1]), ('B', lambda R: R[1])):
        S = [R for R in rows if f(R)]; nontriv = [R for R in S if R[6]]
        print(f'Case {lab}: {len(S)} kappa, {len(nontriv)} with nonempty need; prop holds {sum(R[7] for R in S)}/{len(S)}; '
              f'law {sum(R[8] for R in S)}/{len(S)}; MOB term-by-term {sum(R[10] for R in nontriv)}/{len(nontriv)}; '
              f'IE term-by-term {sum(R[11] for R in nontriv)}/{len(nontriv)}')
    print('caseB1 (kappa_1-only) count', sum(R[2] for R in rows), 'caseB count', sum(R[1] for R in rows))
if __name__ == '__main__': main()
