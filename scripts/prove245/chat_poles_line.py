# Day 245: exact univariate test of root-divisor regularity of structure constants, Ahat vs raw Ã.
# q, spectator w's = random rationals; one variable t free; divisor point t0. Pole iff denominator(t0)=0.
import sys, pickle, itertools, random, sympy as sp
sys.path.insert(0,'/home/agent/projects/work-in-progress/scripts/day242/q484')
from engine import q, subsets
N=int(sys.argv[1]); seed=int(sys.argv[2]) if len(sys.argv)>2 else 1
S=subsets(N); w=sp.symbols('w1:%d'%(N+1)); t=sp.Symbol('t')
A=pickle.load(open('/home/agent/projects/work-in-progress/scripts/wake245/A_full_N%d.pkl'%N,'rb'))
random.seed(seed)
divs=[]
for a in range(1,N+1):
    divs += [('w%d=+1'%a, a, 1), ('w%d=-1'%a, a, -1)]
for a,b in itertools.combinations(range(1,N+1),2):
    divs += [('H%d%d'%(a,b), b, ('inv',a)), ('w%d=w%d'%(a,b), b, ('eq',a))]
for name,var,pt in divs:
    sub={q:sp.Rational(random.randint(2,40),random.randint(41,97))}
    for k in range(1,N+1):
        if k!=var: sub[w[k-1]]=sp.Rational(random.randint(2,60),random.randint(61,97))*random.choice([1,-1])
    sub[w[var-1]]=t
    t0 = pt if not isinstance(pt,tuple) else (1/sub[w[pt[1]-1]] if pt[0]=='inv' else sub[w[pt[1]-1]])
    res={}
    for lab,usepsi in [('Ahat',True),('Araw',False)]:
        mult={M:(sp.prod([(1-q*w[m-1]**2)/(1-q) for m in M]) if usepsi else 1) for M in S}
        Mat=sp.Matrix([[sp.cancel((A[M,I]*mult[M]).subs(sub)) for I in S] for M in S])
        bad=0
        for lam,mu in itertools.combinations_with_replacement(S,2):
            F=sp.Matrix([Mat[i,S.index(lam)]*Mat[i,S.index(mu)] for i in range(len(S))])
            c=Mat.LUsolve(F)  # triangular-ish, exact
            for cc in c:
                d=sp.fraction(sp.cancel(sp.together(cc)))[1]
                if sp.Poly(d,t).eval(t0)==0: bad+=1
        res[lab]=bad
    print(name, res, flush=True)
