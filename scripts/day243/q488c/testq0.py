import sympy as sp,pickle
from engine import subsets
from q0formula import phi0
N=3; A0,A1=pickle.load(open('A01_N3.pkl','rb')); w=sp.symbols('w1:4'); bad=0
for M in subsets(N):
  for I in subsets(N):
    if sp.cancel(phi0(I,M,w)-A0[M,I])!=0: bad+=1; print('bad',M,I)
print('q0 formula mismatches',bad)
