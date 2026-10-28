import sys; sys.path.insert(0,'prove262'); from q46 import *
m=int(sys.argv[1]); a=[2*(m-i)+1 for i in range(m)]  # 2*rho_check = (2m-1,...,1)
for k in strict(m):
    R=[r for r in strict(m) if sub(r,k)]; pred=len({2*sum(r)-len(r) for r in R})
    rk,n,_=rank_for(k,a,m,qv=7); print(k,n,rk,pred,'OK' if rk==pred else 'MISMATCH',flush=True)
