from mp import *
import random
def det(Mx):
    n=len(Mx); A=[r[:] for r in Mx]; d=1
    for i in range(n):
        p=next((r for r in range(i,n) if A[r][i]%P),None)
        if p is None: return 0
        if p!=i: A[i],A[p]=A[p],A[i]; d=-d
        d=d*A[i][i]%P; iv=inv(A[i][i])
        for r in range(i+1,n):
            f=A[r][i]*iv%P
            for c in range(i,n): A[r][c]=(A[r][c]-f*A[i][c])%P
    return d%P
for n in (5,):
  s=Sys(n,7); wv=[x.c[0] for x in s.w]; R=random.Random(3); rats=set()
  for t in range(4):
    g=[R.randrange(1,P) for _ in range(n)]
    f={}
    for M in s.subs:
        v=1
        for m in M: v=v*g[m-1]%P
        f[M]=S.const(v)
    l=s.coord(f,tuple(range(1,n+1))).coef(0)
    D=det([[(g[i]*inv(pow(wv[i],k,P))-pow(wv[i],k,P))%P for k in range(n)] for i in range(n)])
    rats.add(l*inv(D)%P)
  print(n,len(rats))
