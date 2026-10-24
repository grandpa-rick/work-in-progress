# conjectured closed form T' of the half-step T=(I+C)^{1/2}
import sys
exec(open('tshape.py').read().split("tab=Counter()")[0])
def Tp(th):
    v=Fr(1)
    for cc in comps(th):
        dg=[b for b in cc if b[0]==b[1]]
        if not dg:
            if len(cc)>1: return Fr(0)
            continue
        if len(dg)>1: return Fr(0)
        # ribbon? no 2x2
        if any((r+1,c) in cc and (r,c+1) in cc and (r+1,c+1) in cc for r,c in cc): return Fr(0)
        # diagonal box must be the SW end: path from it NE
        d=dg[0]
        # NE end = box with max column then min row
        cmax=max(c for r,c in cc); e=min((r,c) for r,c in cc if c==cmax)
        if len(cc)==1: sg=-1
        else: sg= 1 if (e[0],e[1]-1) in cc else -1
        v*=sg*Fr(-1,2)**len(cc)
    return v
ok=bad=0
for k in KS:
  for r in KS:
    if not boxes(r)<=boxes(k): continue
    th=boxes(k)-boxes(r); t=T[idx[k]][idx[r]]
    if Tp(th)==t: ok+=1
    else:
        bad+=1
        if bad<10: print('X',k,r,sorted(th),t,Tp(th))
print('ok',ok,'bad',bad)
