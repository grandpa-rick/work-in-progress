# Lemma I (inversion duality): Psi_{k,th}(M)(1/w) = (-1)^{|k|} e_M^{2r} W_k^{-1} Psi_{k,thv}(M)(w), thv_j = k_j-2j+1-th_j (j 1-based)
import sys, random, itertools
sys.path.insert(0,'../prove249'); 
from fractions import Fraction as Fr
exec(open('../prove249/q532.py').read().split('print("== (1)')[0])
random.seed(7)
def chk(kap,th,M):
    x=pt(max(M)+1); xi={k:1/v for k,v in x.items()}
    r=len(kap); thv=tuple(kap[j]-2*(j+1)+1-th[j] for j in range(r))
    eM=1
    for m in M: eM*=x[m]
    Wk=1
    for kj in kap:
        for k in range(1,kj): Wk*=x[k]
    L=Psi(kap,th,M,xi); R=(-1)**sum(kap)*eM**(2*r)/Wk*Psi(kap,thv,M,x)
    return L==R, thv
for kap,th,M in [((1,),(0,),(1,2,3)),((2,),(1,),(1,2,4)),((1,3),(0,2),(1,2,3,4)),((2,4),(1,-1),(2,3,4,5)),((1,2,5),(0,0,3),(1,2,4,5,6)),((2,3),(4,4),(1,2,3,5)),((1,3,4),(2,-1,5),(1,3,4,5))]:
    print(kap,th,M,chk(kap,th,M))
# control: wrong thv (no +1)
x=pt(6); xi={k:1/v for k,v in x.items()}
kap,th,M=(2,4),(1,-1),(2,3,4,5); thv=tuple(kap[j]-2*(j+1)-th[j] for j in range(2))
eM=1
for m in M: eM*=x[m]
Wk=x[1]*x[1]*x[2]*x[3]
print('control (thv without +1) differs:',Psi(kap,th,M,xi)!=(-1)**6*eM**4/Wk*Psi(kap,thv,M,x))
