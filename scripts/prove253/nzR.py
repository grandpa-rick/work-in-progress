import sys
from fractions import Fraction as Fr
from check import sp, dem, interp, r
N=int(sys.argv[1]); ok=tot=0; rest=[]
for k in sp(N):
    if not k: continue
    d=sum(k); vals=[dem(k,N,j,'C')[0] for j in range(d+1)]
    kp=k[1:]; Wp=(r(kp)+1)//2
    for a in range(1,N+1):
        v=interp(vals,Fr(-2*a-1,2)); M=2*a+1
        if 2*a<=k[0]:
            tot+=1; pred=0 if a<=Wp else (-1)**(d-k[0])
            ok+= (v==0 and pred==0) or (v!=0 and pred!=0 and (v>0)==(pred>0))
        elif a>Wp:
            rest.append((k,M,'+' if v*(-1)**d>0 else ('-' if v!=0 else '0'), 'k1-par',k[0]%2))
print(f'm={N}: Prop-R regime 2a<=k1: zero iff a<=ceil(r(k\')/2), sign (-1)^(|k|-k1): {ok}/{tot}')
import collections
print('beyond regime (M>=k1+2, M<2m+3), sign rel (-1)^|k|:',collections.Counter((x[2],x[1]-x[0][0]) for x in rest))
