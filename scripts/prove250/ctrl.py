# negative controls for chk1: wrong e^alpha sign; N without w_L^2; sigma-lemma with z instead of z^2
import sys
sys.argv=[sys.argv[0],sys.argv[1]]
exec(open('chk1.py').read().split('# (i)')[0])
bad1=bad2=bad3=0
for nup in subsets(n):
    for s in range(1,n+1):
        if s<n:
            if not ((s+1) in nup and s not in nup): continue
            nu=sorted([s if x==s+1 else x for x in nup])
        else:
            if 1 not in nup: continue
            nu=[x for x in nup if x!=1]
        ea=ema(s,w0)  # WRONG: e^{-alpha}
        for M in subsets(n):
            if gam(nup,sL(s,M),sw(s,w0))!=ea*gam(nup,M,w0)-(ea-1)*gam(nu,M,w0): bad1+=1
            if gam(nup,M,sw(s,w0))!=(1/ema(s,w0))*gam(nup,M,w0)-(1/ema(s,w0)-1)*gam(nu,M,w0): bad3+=1  # no label action
def N2(L,w): return cc(L,w)/(F(L,L,w)*gam(comp(L),comp(L),w))
for s in range(1,n+1):
    for L in subsets(n):
        if N2(sL(s,L),sw(s,w0))!=N2(L,w0): bad2+=1
print('ctrl wrong e^a:',bad1,' ctrl no label action:',bad3,' ctrl N w/o w_L^2:',bad2)
