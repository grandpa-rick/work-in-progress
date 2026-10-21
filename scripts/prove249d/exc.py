import sys, sympy as sp, re
c=sp.Symbol('c')
for fns,n in [(['../prove249c/hpC3x.txt'],3),(['../prove249c/hpC4x_0.txt','../prove249c/hpC4x_1.txt'],4),(['../prove249c/hpC5s.txt'],5)]:
    cnt={}
    for fn in fns:
        for line in open(fn):
            m=re.match(r'(\[.*?\]) (\[.*?\]) (\[.*?\]) d (\d+) deg (\d+).* P (.*)$',line.strip())
            if not m: continue
            lam,mu,kap=[tuple(eval(m.group(i))) for i in (1,2,3)]; d=int(m.group(4)); P=sp.sympify(m.group(6))
            e=sum(lam)+sum(mu)-sum(kap)
            key='d=0' if d==0 else ('d=e>0' if d==e else 'EXC')
            cnt[key]=cnt.get(key,0)+1
            if key=='EXC': print(n,lam,mu,kap,'e',e,'d',d,'P',sp.factor(P))
            if key=='d=e>0' and sp.degree(P,c)!=0: print('nonconst d=e',lam,mu,kap,P)
    print('n',n,cnt)
