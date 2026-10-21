import sys, sympy as sp, re
c=sp.Symbol('c')
for fns,n in [(['../prove249c/hpC3x.txt'],3),(['../prove249c/hpC4x_0.txt','../prove249c/hpC4x_1.txt'],4)]:
    rows=[]
    for fn in fns:
        for line in open(fn):
            m=re.match(r'(\[.*?\]) (\[.*?\]) (\[.*?\]) d (\d+) deg (\d+).* P (.*)$',line.strip())
            if not m: continue
            lam,mu,kap=[tuple(eval(m.group(i))) for i in (1,2,3)]; d=int(m.group(4)); P=sp.sympify(m.group(6))
            U=2*n-2
            rows.append((lam,mu,kap,d,[sp.sign(P.subs(c,U+k)) for k in (1,2,3)]))
    sig={}; incons=0
    for lam,mu,kap,d,s in rows:
        if lam==() :
            if len(set(s))!=1: print('not one-signed?',lam,mu,kap,s)
            sig.setdefault(mu,set()).add(s[0])
    print('n',n,'sigma(mu) from lam=():',{k:v for k,v in sig.items()})
    ok=bad=0; badl=[]
    for lam,mu,kap,d,s in rows:
        if lam in sig and mu in sig and len(sig[lam])==1 and len(sig[mu])==1:
            pred=list(sig[lam])[0]*list(sig[mu])[0]
            if d>0: continue
            if s[0]==pred: ok+=1
            else: bad+=1; badl.append((lam,mu,kap,d))
    print(' product rule ok',ok,'bad',bad,badl[:8])
