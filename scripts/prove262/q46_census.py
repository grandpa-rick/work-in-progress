# Test: rank{spec_a ch_rho : rho⊆kappa} == #{(eps_rho,a): rho⊆kappa} (bins), and pi_q => |kappa|+1; x=1 => |kappa|+1
import sys, random, itertools, sympy as sp
sys.path.insert(0,'prove262'); from q46 import *
m=int(sys.argv[1]); random.seed(int(sys.argv[2]) if len(sys.argv)>2 else 1)
def eps(r): return [(-1 if (m+1-i) in r else 1) for i in range(1,m+1)]
def regular(a): return all(x!=0 for x in a) and all(a[i]!=a[j] and a[i]!=-a[j] for i in range(m) for j in range(i+1,m))
grads=[list(range(m,0,-1))]
while len(grads)<5:
    a=[random.randint(-9,9) for _ in range(m)]
    if regular(a): grads.append(a)
# a few regular gradings ON collision hyperplanes
grads+= ([[3,2,1],[5,3,2],[4,3,1],[7,4,3]] if m==3 else [[5,3,2,1],[6,4,3,1],[7,4,2,1],[9,5,3,1]])
bad=0
for a in grads:
    assert regular(a),a
    for k in strict(m):
        R=[r for r in strict(m) if sub(r,k)]
        bins=len({sum(e*x for e,x in zip(eps(r),a)) for r in R})
        rk,n,_=rank_for(k,a,m,qv=7)   # numeric q=7 gives a lower bound; equality with bins (an upper bound) certifies
        if rk!=bins: bad+=1; print('MISMATCH a=',a,k,rk,bins)
        if a==list(range(m,0,-1)) and rk!=sum(k)+1: bad+=1; print('PI_q MISMATCH',k,rk)
    print('a=',a,'done')
for k in strict(m):
    R=[r for r in strict(m) if sub(r,k)]
    M=sp.Matrix([[sum(chw(r,j,m).values()) for r in R] for j in range(2*sum(k)+3)])
    if M.rank()!=sum(k)+1: bad+=1; print('x=1 MISMATCH',k)
print('m=',m,'mismatches:',bad)
