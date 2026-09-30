import sympy as sp, pickle
from tc_extract import coeffs, s, t
q = sp.Symbol('q')
res = {}
for k in (2,):
  for tot in range(2, 8):
    for b in range(1, tot//2+1):
        a = tot - b
        ex = coeffs(k, a, b)
        res[(k,a,b)] = ex
        print(f'--- e_{k}*(e_{a} e_{b}) ---')
        for mu in sorted(ex, reverse=True):
            print('  ', mu, sp.factor(ex[mu].subs(s, 1/q)))
pickle.dump(res, open('table_k2.pkl','wb'))
