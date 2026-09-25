import json, glob, os, sympy as sp
q, t, u = sp.symbols('q t u')
D5='/home/agent/projects/proofs/scripts/day205'; D6='/home/agent/projects/proofs/scripts/day206'
def load(k=4):
    data={}
    for fn in glob.glob(f'{D5}/k{k}_r*_m*.json')+[f for f in glob.glob(f'{D6}/k{k}_r*_m*.json') if '_s' not in os.path.basename(f)]:
        b=os.path.basename(fn)[:-5].split('_'); r,m=int(b[1][1:]),int(b[2][1:])
        data[(r,m)]=sp.sympify(json.load(open(fn))[str((r+k,))])
    return data
def qint(n,x=t): return sum(x**i for i in range(n))
if __name__=='__main__':
    data=load()
    rs=sorted({r for r,m in data})
    for r in rs:
        ms=sorted(m for rr,m in data if rr==r)
        base=data[(r,ms[0])]
        st=[sp.cancel(data[(r,m)]-base)==0 for m in ms[1:]]
        print(f'r={r} ms={ms} stable={st}')
    for r in rs:
        tau=data[(r,min(m for rr,m in data if rr==r))]
        f=sp.factor_list(sp.cancel(tau))
        print(f'\nr={r}: const={f[0]}')
        for fac,e in f[1]:
            d=sp.Poly(fac,q,t).degree_list() if fac.free_symbols else None
            print(f'   ({"deg(q,t)="+str(d)}) ^{e}:', fac if len(str(fac))<80 else str(fac)[:80]+'...')
