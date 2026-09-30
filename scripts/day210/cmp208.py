import pickle, sympy as sp, time
from tc_extract import coeffs, s, t
q = sp.Symbol('q')
for k, f in [(2,'data_k2_N4.pkl'), (3,'data_k3_N3.pkl')]:
    D = pickle.load(open('../day208/'+f,'rb'))
    for (a,b), d in sorted(D.items()):
        t0=time.time(); ex = coeffs(k, b, a)
        # detect variable convention in data
        fs = set().union(*[sp.sympify(v).free_symbols for v in d.values()])
        dd = {m: sp.sympify(v) for m,v in d.items()}
        if q in fs: dd = {m: v.subs(q, 1/s) for m,v in dd.items()}
        bad = [m for m in set(ex)|set(dd) if sp.cancel(ex.get(m,0)-dd.get(m,0)) != 0]
        print(k,(a,b),'OK' if not bad else 'MISMATCH '+str(bad), f'{time.time()-t0:.1f}s', flush=True)
