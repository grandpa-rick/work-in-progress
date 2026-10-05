from lead2 import *
import pickle, itertools
D = pickle.load(open('../wake223/leads_n6.pkl', 'rb'))
t0 = Fr(3)
ok = bad = 0
for (lam, mu), (v, Ld) in D.items():
    if len(lam) != 3 or v != 2: continue
    vals = set()
    for perm in set(itertools.permutations(lam)):
        vals.add(second(*perm, t0).get(mu, 0))
    want = Ld.subs(sympy.symbols('t'), 3)
    good = len(vals) == 1 and list(vals)[0] == want
    ok += good; bad += not good
    print(lam, mu, vals, want, good)
print(ok, bad)
