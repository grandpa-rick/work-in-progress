# LS-path model for LG: H_kappa(j) =? #{g: shifted(kappa) -> {0..2j} order-preserving, g even on diagonal boxes}
import sys, itertools
sys.path.insert(0, '/home/agent/projects/work-in-progress/scripts/wake251/orderpoly')
from match import Hs_for, j
from orderpoly import shifted, young, boxless
import sympy as sp
from functools import lru_cache

def count(P, N, parity):  # order-preserving maps P->{0..N}; parity(box) in {None,0,1}; DP over order ideals
    from orderpoly import ideals
    J = ideals(P, boxless); idx = {I: i for i, I in enumerate(J)}
    vec = {I: 1 for I in J}  # I_0 = {g<=0}: any ideal satisfying parity at v=0
    vec = {I: 1 for I in J if all(parity(b) in (None, 0) for b in I)}
    for v in range(1, N + 1):
        nv = {}
        for I, c in vec.items():
            for K in J:
                if I <= K and all(parity(b) in (None, v % 2) for b in K - I):
                    nv[K] = nv.get(K, 0) + c
        vec = nv
    return vec.get(frozenset(P), 0)

def model(kap, jj, P_of=shifted, diag_even=True):
    lam = sorted(kap, reverse=True); P = P_of(lam)
    par = (lambda b: 0 if b[0] == b[1] else None) if diag_even else (lambda b: None)
    return count(P, 2 * jj, par)

if __name__ == '__main__':
    for m in (3, 4, 5):
        Hs = Hs_for(m); ok = 0; bad = []
        okneg = 0
        for k, H in Hs.items():
            d = sum(k)
            vals = [model(k, x) for x in range(d + 2)] 
            good = all(H.subs(j, x) == v for x, v in enumerate(vals))
            ok += good
            if not good: bad.append((k, vals, [H.subs(j, x) for x in range(len(vals))]))
            # negative control: same model on the UNSHIFTED Young diagram (diag boxes (r,r) still even)
            vn = [model(k, x, P_of=young) for x in range(3)]
            okneg += all(H.subs(j, x) == v for x, v in enumerate(vn))
        print(f'm={m}: LS-model (shifted, diag even, N=2j) matches {ok}/{len(Hs)} on j=0..min(d+1,2); Young-diagram control matches {okneg}/{len(Hs)}; bad {bad[:3]}')
