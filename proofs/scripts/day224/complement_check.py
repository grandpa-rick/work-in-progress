"""Check Box Complement: c_{N^l - lam, N^l - mu} = s^{N C(l,2) - (l-1)|lam|} c_{lam,mu}, on full (s,t) data n<=6."""
import pickle, sympy
C = {}
for n in range(2, 7):
    f = '../wake222/gprime/cst_n%d.pkl' % n if n < 6 else '../wake222/gprime/cst_n6_skip1n.pkl'
    C.update(pickle.load(open(f, 'rb')))
s, t = sympy.symbols('s t')
ok = bad = 0
for (lam, mu), c in sorted(C.items()):
    c = sympy.sympify(c)
    for l in range(len(lam), len(lam)+3):
        if len(mu) > l: continue
        lp = tuple(lam)+(0,)*(l-len(lam)); mp = tuple(mu)+(0,)*(l-len(mu))
        for N in range(max(lp+mp), max(lp+mp)+3):
            lh = tuple(sorted([N-v for v in lp if N-v], reverse=True)); mh = tuple(sorted([N-v for v in mp if N-v], reverse=True))
            if sum(lh) == 0 or (lh, mh) not in C: continue
            c2 = sympy.sympify(C[(lh, mh)])
            sig = N*l*(l-1)//2 - (l-1)*sum(lam)
            good = sympy.simplify(c2 - s**sig*c) == 0
            ok += good; bad += not good
            if not good: print('FAIL', lam, mu, l, N, c, c2)
print('complement: ok', ok, 'bad', bad)
