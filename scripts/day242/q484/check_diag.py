import sympy as sp, sys
from engine import *
def rootprod(X, N, I, w, qq=q):
    # type C: eq typecfirstblock/typecsecondblock (Motegi L1017-1046)
    Ic = [j for j in range(1, N+1) if j not in I]
    W_ = lambda k: w[k-1]
    p = sp.Integer(1)
    for a in range(len(I)):
        for b in range(a, len(I)):
            p *= (1 - 1/(W_(I[a])*W_(I[b])))
    for j in Ic:
        for il in I:
            if il > j: p *= (1 - W_(j)/W_(il))
    for a in range(len(Ic)):
        for b in range(a, len(Ic)):
            p *= (1 - qq*W_(Ic[a])*W_(Ic[b]))
    for j in Ic:
        for il in I:
            if il < j: p *= (1 - qq*W_(j)/W_(il))
    return p
X = 'C'
for N in (1, 2, 3):
    ok = bad = 0
    for I in subsets(N):
        d, w = W(X, N, I, I)
        ok += sp.cancel(d - rootprod(X, N, I, w)) == 0
        # planted: q -> -q in ascent block must break
        bad += sp.cancel(d - rootprod(X, N, I, w, -q)) == 0
    print(N, 'diag match', ok, '/', 2**N, ' planted(q->-q) matches', bad)
