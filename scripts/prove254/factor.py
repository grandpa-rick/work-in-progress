import sys
exec(open('tprime.py').read().split("ok=bad=0")[0])
def anti(th): return all(not(x!=y and le(x,y)) for x in th for y in th)
Ao=[[Fr(0)]*n for _ in range(n)]; Ae=[[Fr(0)]*n for _ in range(n)]
for k in KS:
  for r in KS:
    if boxes(r)<=boxes(k):
      th=boxes(k)-boxes(r)
      if anti(th):
        Ae[idx[k]][idx[r]]=Fr(1)
        if all(a!=b for a,b in th): Ao[idx[k]][idx[r]]=Fr(1)
def inv(A):  # unitriangular wrt size order
    X=[[Fr(0)]*n for _ in range(n)]
    for i in range(n):
        X[i][i]=Fr(1)
        for j in range(i-1,-1,-1):
            X[i][j]=-sum(A[i][l]*X[l][j] for l in range(j,i))
    return X
AC=mul(Ao,Ae); assert all(AC[i][j]==(I[i][j]+C[i][j]) for i in range(n) for j in range(n)), 'AoAe != I+C'
print('I+C = Ao Ae OK; commute?', mul(Ao,Ae)==mul(Ae,Ao))
U=mul(inv(Ao),T); U2=mul(T,inv(Ao)); W=mul(inv(Ae),T); W2=mul(T,inv(Ae))
for name,M in [('Ao^-1 T',U),('T Ao^-1',U2),('Ae^-1 T',W),('T Ae^-1',W2)]:
    print('==',name)
    for k in KS[:12]:
        if k: print(' ',k,{(r if r else '0'):str(M[idx[k]][idx[r]]) for r in KS if M[idx[k]][idx[r]]!=0 and r!=k})
