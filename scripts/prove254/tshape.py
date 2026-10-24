from half import *
import sys
exec(open('tmat.py').read().split("for k in KS:\n    if k: print")[0].replace("N=int(sys.argv[1])","N=int(sys.argv[1])"))
from collections import Counter
def comps(th):
    th=set(th); out=[]
    while th:
        st=[th.pop()]; cc={st[0]}
        while st:
            (r,c)=st.pop()
            for q in ((r+1,c),(r-1,c),(r,c+1),(r,c-1)):
                if q in th: th.remove(q); cc.add(q); st.append(q)
        out.append(frozenset(cc))
    return out
def canon(cc):  # shape up to translation along diagonal? keep diag info
    r0=min(r for r,c in cc); c0=min(c for r,c in cc)
    return tuple(sorted((r-r0,c-c0,int(r==c)) for r,c in cc))
tab=Counter(); allv={}
for k in KS:
  for r in KS:
    v=T[idx[k]][idx[r]]
    th=boxes(k)-boxes(r)
    if not boxes(r)<=boxes(k): continue
    key=tuple(sorted(canon(c) for c in comps(th)))
    allv.setdefault(key,set()).add(v)
bad=[(k,v) for k,v in allv.items() if len(v)>1]
print('component-multiset determines T?', len(bad)==0, len(allv))
for k,v in bad[:15]: print(k,v)
single={k[0]:list(v)[0] for k,v in allv.items() if len(k)==1}
ok=bad2=miss=0
for key,v in allv.items():
    v=list(v)[0]
    if len(key)<=1: continue
    if all(c in single for c in key):
        p=Fr(1)
        for c in key: p*=single[c]
        if p==v: ok+=1
        else: bad2+=1; print('nonmult',key,v,p)
    else: miss+=1
print('multiplicative',ok,'bad',bad2,'missing singles',miss)
nz={k:v for k,v in single.items() if v!=0}
for k,v in sorted(nz.items(),key=lambda t:len(t[0])): print(len(k),v,k)
