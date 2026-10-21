import re, sympy as sp
from dem import chi
x=sp.Symbol('m'); c=sp.Symbol('c')
for fns,n in [(['../prove249c/hpC3x.txt'],3),(['../prove249c/hpC4x_0.txt','../prove249c/hpC4x_1.txt'],4)]:
    ok=tot=0
    for fn in fns:
        for line in open(fn):
            m=re.match(r'\[\] \[\] (\[.*?\]) d (\d+) .* P (.*)$',line.strip())
            if not m: continue
            kap=eval(m.group(1)); P=sp.sympify(m.group(3)); d=sum(kap)
            h=sp.interpolate([(k,chi(set(kap),n,k)[0]) for k in range(d+2)],x)
            Q=sp.expand((-1)**d*h.subs(x,-c/2-1)); tot+=1; ok+= sp.expand(P-Q)==0
            # even-c floor: value >=1 for c=2k, k>=n
            assert all(Q.subs(c,2*k)>=1 for k in range(n,n+6))
    print('n',n,'engine P == (-1)^|k| chi(X_k,O(-c/2-1)):',ok,'/',tot)
