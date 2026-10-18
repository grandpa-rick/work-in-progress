# Day 245 PROVE: structure constants of Ahat = Ã*Psi at N=3; classify non-monomial denominator factors.
import sys, pickle, itertools, sympy as sp
sys.path.insert(0,'/home/agent/projects/work-in-progress/scripts/day242/q484')
from engine import q, subsets
N=int(sys.argv[1]); S=subsets(N); w=sp.symbols('w1:%d'%(N+1))
A=pickle.load(open('/home/agent/projects/work-in-progress/scripts/wake245/A_full_N%d.pkl'%N,'rb'))
psi=lambda v:(1-q*v**2)/(1-q)
Ps={M:sp.prod([psi(w[m-1]) for m in M]) for M in S}
for name,mult in [('Ahat',Ps),('Araw',{M:1 for M in S})]:
    Mat=sp.Matrix([[A[M,I]*mult[M] for I in S] for M in S])  # rows M, cols I
    inv=Mat.inv(method='LU')
    facs={}
    for lam,mu in itertools.combinations_with_replacement(S,2):
        F=sp.Matrix([A[M,lam]*A[M,mu]*mult[M]**2 for M in S])
        c=inv*F
        for nu,cc in zip(S,c):
            cc=sp.factor(sp.cancel(cc)); d=sp.fraction(cc)[1]
            for f,e in sp.factor_list(d)[1]:
                if f.free_symbols & set(w) and len(sp.Add.make_args(sp.expand(f)))>1:
                    facs.setdefault(sp.factor(f),[]).append((lam,mu,nu))
    print(name, 'non-monomial denominator factors:', {k:len(v) for k,v in facs.items()})
