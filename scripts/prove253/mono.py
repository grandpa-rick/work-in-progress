import sys
sys.argv=[sys.argv[0],sys.argv[1]]
from check import sp, dem, interp, r, N
from fractions import Fraction as Fr
for typ in 'CB':
    ok=tot=0; bad=[]
    for k in sp(N):
        if not k: continue
        d=sum(k); vals=[dem(k,N,j,typ)[0] for j in range(d+1)]
        s=[(-1)**d*interp(vals,-j) for j in range(1,d+8)]
        W=(r(k)+1)//2 if typ=='C' else r(k)
        good=all(x>=0 for x in s) and all(s[i+1]>=s[i] for i in range(len(s)-1)) and all(s[j-1]==0 for j in range(1,W+1)) and s[W]>0
        tot+=1; ok+=good
        if not good: bad.append((k,s[:W+3]))
    print(f'm={N} {typ}: s_j>=0, nondecreasing, =0 iff j<=W: {ok}/{tot}', bad[:3])
