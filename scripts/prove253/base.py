# B(kappa), kappa_1 = 2b+1 odd > kappa'_1: need s_{b+1/2}(kappa') >= 2^{-(2b+1)} s_{b+1}(kappa'), s_x(k)=(-1)^|k| H_k(-x)
import sys
from fractions import Fraction as Fr
from check import sp, dem, interp, r
N=int(sys.argv[1]); worst=[]
for k in sp(N):
    d=sum(k); vals=[dem(k,N,j,'C')[0] for j in range(d+1)] if k else [1]
    s=lambda x:(-1)**d*interp(vals,x) if k else Fr(1)
    k1=k[0] if k else 0
    for b in range(0,N+3):
        if 2*b+1<=k1: continue
        L=s(Fr(-2*b-1,2)); Rr=s(-(b+1))
        worst.append((float(L*2**(2*b+1)/Rr) if Rr else float('inf'), k, b, L, Rr))
worst.sort(key=lambda x:x[0])
for w in worst[:12]: print(f'{w[0]:.4f}', w[1], 'b',w[2], 's_half',w[3], 's_int',w[4])
print('violations (ratio<1):', sum(1 for w in worst if w[0]<1), 'of', len(worst))
