# Wake 257: independent check of Conjecture conj:T(b) in lg-og-negative-degree/main.tex.
# T computed from the paper's definitions only: A_o, A_e (antichain matrices, Sec. "odd-zero law"),
# A_oA_e = I + C (Lemma layers (a) = matrix of eq:star), T = sum_k binom(1/2,k) C^k (C nilpotent).
# Rules: OLD = tex as printed before fix; NEW = corrected printed text; CTRL = deliberately wrong.
import itertools, sys
from fractions import Fraction as Fr

def strict(m):
    return [tuple(sorted(c, reverse=True)) for L in range(m+1) for c in itertools.combinations(range(1, m+1), L)]
def D(k):  # boxes (a,b), 0<=a<l, a<=b<=a+k_{a+1}-1
    return frozenset((a, b) for a, ka in enumerate(k) for b in range(a, a+ka))
le = lambda x, y: x[0] <= y[0] and x[1] <= y[1]
def antichain(th): return all(not (x != y and le(x, y)) for x in th for y in th)
def comps(th):
    th = set(th); out = []
    while th:
        st = [th.pop()]; cc = {st[0]}
        while st:
            r, c = st.pop()
            for q in ((r+1, c), (r-1, c), (r, c+1), (r, c-1)):
                if q in th: th.remove(q); cc.add(q); st.append(q)
        out.append(frozenset(cc))
    return out
has2x2 = lambda cc: any({(r+1, c), (r, c+1), (r+1, c+1)} <= cc for r, c in cc)
def ne_end(cc):
    cmax = max(c for r, c in cc); return min((r, c) for r, c in cc if c == cmax)

def rule(th, variant):
    v = Fr(1)
    for cc in comps(th):
        dg = [b for b in cc if b[0] == b[1]]
        if not dg:
            if len(cc) == 1: continue          # off-diagonal single box: weight 1
            return Fr(0)                        # off-diagonal, >=2 boxes: 0
        if variant == 'OLD':
            # printed: exactly one diagonal box -> ribbon D, weight +-2^{-|D|}, + iff left nbr of NE end in D.
            # printed text silent on >=2 diagonal boxes and on non-ribbons: we read "is a ribbon" as an
            # assertion, so evaluate the sign formula regardless (literal), and 0 for >=2 diagonal (unstated -> 0).
            if len(dg) > 1: return Fr(0)
            e = ne_end(cc); sg = 1 if (e[0], e[1]-1) in cc else -1
            v *= sg * Fr(1, 2**len(cc))
        elif variant == 'NEW':
            # new printed text: >=2 diagonal boxes -> 0; one diagonal box, not a ribbon -> 0;
            # ribbon -> eps*(-1/2)^{|D|}, eps=+1 iff left neighbour of NE end box lies in D.
            if len(dg) >= 2 or has2x2(cc): return Fr(0)
            e = ne_end(cc); eps = 1 if (e[0], e[1]-1) in cc else -1
            v *= eps * Fr(-1, 2)**len(cc)
        elif variant == 'CTRL':
            # wrong on purpose: NEW but without the 2x2 zero case and with unsigned weight
            if len(dg) >= 2: return Fr(0)
            v *= Fr(1, 2**len(cc))
    return v

def run(m):
    KS = sorted(strict(m), key=sum); n = len(KS); idx = {k: i for i, k in enumerate(KS)}
    Ds = [D(k) for k in KS]
    Ao = [[Fr(int(Ds[j] <= Ds[i] and antichain(Ds[i]-Ds[j]) and all(a != b for a, b in Ds[i]-Ds[j]))) for j in range(n)] for i in range(n)]
    Ae = [[Fr(int(Ds[j] <= Ds[i] and antichain(Ds[i]-Ds[j]))) for j in range(n)] for i in range(n)]
    mul = lambda A, B: [[sum(A[i][l]*B[l][j] for l in range(n) if A[i][l]) for j in range(n)] for i in range(n)]
    M = mul(Ao, Ae)
    C = [[M[i][j]-(i == j) for j in range(n)] for i in range(n)]
    T = [[Fr(int(i == j)) for j in range(n)] for i in range(n)]; P = [r[:] for r in T]; b = Fr(1)
    for k in range(1, n+1):
        P = mul(P, C); b = b*(Fr(1, 2)-(k-1))/k
        if all(x == 0 for r in P for x in r): break
        T = [[T[i][j]+b*P[i][j] for j in range(n)] for i in range(n)]
    # sanity: T^2 = A_oA_e, T unipotent lower-triangular w.r.t. inclusion
    assert mul(T, T) == M
    assert all(T[i][j] == 0 for i in range(n) for j in range(n) if not Ds[j] <= Ds[i])
    pairs = [(i, j) for i in range(n) for j in range(n) if Ds[j] <= Ds[i]]
    res = {}
    for var in ('OLD', 'NEW', 'CTRL'):
        ok = sum(rule(Ds[i]-Ds[j], var) == T[i][j] for i, j in pairs)
        res[var] = ok
    # catalogue: components with >=2 diagonal boxes / non-ribbon one-diag, single-component cases
    multi = [(KS[i], KS[j], T[i][j]) for i, j in pairs if len(comps(Ds[i]-Ds[j])) == 1
             and sum(a == b for a, b in Ds[i]-Ds[j]) >= 2]
    nonrib = [(KS[i], KS[j], T[i][j]) for i, j in pairs if len(comps(Ds[i]-Ds[j])) == 1
              and sum(a == b for a, b in Ds[i]-Ds[j]) == 1 and has2x2(Ds[i]-Ds[j])]
    # part (a): T s_0 = s_{1/2}: s_0 = (-1)^{|k|}H(0) = (-1)^{|k|}; s_{1/2} anchor (-1)^{|k|}2^{-k_1}
    s0 = [Fr((-1)**sum(k)) for k in KS]
    a_ok = all(sum(T[i][j]*s0[j] for j in range(n)) == Fr((-1)**sum(KS[i]), 2**(KS[i][0] if KS[i] else 0)) for i in range(n))
    print(f"m={m}: pairs={len(pairs)}  OLD {res['OLD']}/{len(pairs)}  NEW {res['NEW']}/{len(pairs)}  CTRL {res['CTRL']}/{len(pairs)}"
          f"  | 1-comp >=2diag: {len(multi)} pairs, T values {sorted(set(map(str,(t for *_,t in multi))))}"
          f"  | 1-comp 1diag non-ribbon: {len(nonrib)}, T values {sorted(set(map(str,(t for *_,t in nonrib))))}  | T s0 = s_1/2 anchor: {a_ok}")
    ex = [(KS[i], KS[j], str(T[i][j]), str(rule(Ds[i]-Ds[j], 'OLD'))) for i, j in pairs if rule(Ds[i]-Ds[j], 'OLD') != T[i][j]][:4]
    print("   OLD failures e.g. (kappa, nu, T, OLD):", ex)

if __name__=="__main__":
  for m in map(int, sys.argv[1:] or [3,4,5,6]): run(m)
