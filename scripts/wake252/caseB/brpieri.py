# Buch-Ravikumar (arXiv 1005.2605, Crelle 668 (2012)) K-Pieri for LG(n,2n): Cor 5.4 + eq.(11)
#   O^lam . O^1 = sum_nu c^nu_{lam,1} O^nu,  c^nu_{lam,p} = C(nu/lam, p)
#   C(th,p) = sum_{th' <= phi <= th} (-1)^{|th|-|phi|} h(N'(phi), d(phi)-p),  h(a,b)=sum_{j<=b} (-1)^j 2^{a-j} C(a,j)
# Codimension indexing X^lam; our dimension-indexed X_kappa = X^{kappa^vee}, kappa^vee = [n] \ kappa.
# [O_D] = [O_{X_kappa}] . O^1 (D = hyperplane section; Schubert divisor X^1 is a hyperplane section).
import sys, itertools, ast, re
from math import comb
def boxes(l): return {(i, i + j) for i, li in enumerate(l) for j in range(li)}   # (row, col), diagonal: row==col
def h(a, b): return sum((-1) ** j * 2 ** (a - j) * comb(a, j) for j in range(0, min(a, b) + 1)) if b >= 0 else 0
def se_corners(th): return {(r, c) for (r, c) in th if (r + 1, c) not in th and (r, c + 1) not in th}
def se_rim(th): return {(r, c) for (r, c) in th if not any((r2 > r and c2 > c) for (r2, c2) in th)}
def comps(th):
    th = set(th); out = []
    while th:
        st = [th.pop()]; cc = set(st)
        while st:
            r, c = st.pop()
            for nb in ((r + 1, c), (r - 1, c), (r, c + 1), (r, c - 1)):
                if nb in th: th.remove(nb); cc.add(nb); st.append(nb)
        out.append(cc)
    return out
def Nprime(th): return sum(1 for cc in comps(th) if not any(r == c for r, c in cc))
def C(th, p):
    cor = se_corners(th); base = th - cor; tot = 0
    for L in range(len(cor) + 1):
        for S in itertools.combinations(sorted(cor), L):
            phi = base | set(S)
            tot += (-1) ** (len(th) - len(phi)) * h(Nprime(phi), len(se_rim(phi)) - p)
    return tot
def dual(k, n): return tuple(sorted(set(range(1, n + 1)) - set(k), reverse=True))
def pieri_dim(k, n):
    """[O_{X_kappa}].O^1 in the dimension-indexed basis {X_rho}."""
    lam = dual(k, n); bl = boxes(lam); out = {}
    for L in range(n + 1):
        for c in itertools.combinations(range(1, n + 1), L):
            nu = tuple(sorted(c, reverse=True)); bn = boxes(nu)
            if bl <= bn and bl != bn:
                v = C(bn - bl, 1)
                if v: out[dual(nu, n)] = v
    return out
if __name__ == '__main__':
    f = sys.argv[1]; n = int(sys.argv[2]); ok = bad = 0
    for L in open(f):
        if not L.startswith('('): continue
        k = ast.literal_eval(re.match(r'(\(.*?\))', L).group(1))
        K = ast.literal_eval(re.search(r"'?K'? (\{.*?\})", L).group(1))
        K = {(tuple(int(ch) for ch in s) if s != '0' else ()): v for s, v in K.items()}
        P = pieri_dim(k, n)
        if P == K: ok += 1
        else: bad += 1; print('MISMATCH', k, K, P)
    print(f'n={n}: Buch-Ravikumar Cor 5.4 reproduces computed [O_D]: {ok}/{ok+bad}')
