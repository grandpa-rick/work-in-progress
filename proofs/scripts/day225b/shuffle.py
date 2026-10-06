import sympy as sp, itertools
t, w = sp.symbols('t w')
def Sh(A, B):
    tot = 0
    for pos in itertools.combinations(range(A+B), A):
        S = set(pos); word = ['U' if m in S else 'V' for m in range(A+B)]
        ia = ib = 0; idx = []
        for L in word:
            if L == 'U': idx.append(('U', ia)); ia += 1
            else: idx.append(('V', ib)); ib += 1
        P = 1
        for m in range(A+B):
            for l in range(m):
                (L1, k1), (L2, k2) = idx[l], idx[m]
                if L1 == L2: continue
                if L1 == 'U': al, be = k1, k2; P *= (w*t**(be-al)-1)/(w*t**(be-al)-t)
                else: be, al = k1, k2; P *= (1-w*t**(be-al))/(1-w*t**(be-al+1))
        tot += P
    return sp.factor(sp.cancel(tot))
for A in range(1, 4):
    for B in range(1, 4):
        print(A, B, Sh(A, B))
