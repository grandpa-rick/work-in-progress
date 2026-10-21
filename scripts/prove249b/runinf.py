import sys,time,collections
import scert8 as S
S.LIMS=(1,0); S.IOTA=(sys.argv[2]=='1'); N=int(sys.argv[1])
if len(sys.argv)>3: S.LIMS=(1,)
bad=[];tot=0;t=time.time()
for x in S.allA(N):
    tot+=1
    if not S.start('A',N,*x): bad.append(x)
print('A N',N,'LIMS',S.LIMS,'IOTA',S.IOTA,'uncert',len(bad),'/',tot,round(time.time()-t),'s',S.PSTAT,flush=True)
print(collections.Counter(x[2] for x in bad).most_common(6))
