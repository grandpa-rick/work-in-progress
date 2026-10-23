import sys
from fractions import Fraction as Fr
from check import sp, dem, interp
N=int(sys.argv[1]); cache={}
def S(k,x):
    if not k: return Fr(1)
    if k not in cache: cache[k]=[dem(k,N,j,'C')[0] for j in range(sum(k)+1)]
    return (-1)**sum(k)*interp(cache[k],x)
for k in sp(N):
    if not k: continue
    k1=k[0]; kp=k[1:]
    if k1%2==0:
        b=k1//2; print('even',k,'s_{b+1}(k)=',S(k,-(b+1)),' s_b(k\')=',S(kp,-b), ' H_k(-(k1+2)/2-int)?')
    else:
        b=(k1+1)//2; print('odd ',k,'s_{b+1}(k)=',S(k,-(b+1)),' s_b(k\')=',S(kp,-b),' s_{b+1/2}(k)=',S(k,Fr(-2*b-1,2)), 's_{b-1/2}(kp)',S(kp,Fr(-2*b+1,2)))
