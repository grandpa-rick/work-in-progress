import sys; sys.path.insert(0,'../prove247')
import mp; mp.K=12
from mp import S,Sys
import scert
scert.NV=2
n=2;a=1;kap=(1,2);c=2
key=('B',(),(),(),(),kap,(),0,(0,0,0))
for lim in (0,1):
    print('lim',lim,[ (o,sl,sg,nk) for o,sl,sg,nk in scert.summand_opts('B',2,1,None,key,lim)])
# actual: in-term and out-term separately, u=eps or 1/eps
for lim in (0,1):
    s=Sys(2,3); s.w[1]=S(1 if lim==0 else -1,[1]+[0]*(mp.K-1))
    x=s.w[0]; one=S.const(1)
    for M in [(1,),(1,2)]:
        f={M:S.const(1)}
        if M==(1,2):
            u=s.w[1]; g=u*u*(one-x/u)/(one-x*u); f={M:(u**c)/g}
        v=s.coord(f,kap).norm()
        print('lim',lim,'M',M,'order',v.v)
