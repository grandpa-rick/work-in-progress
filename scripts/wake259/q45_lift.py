# Wake 259 Task 2 (Q45 kill test): does the (★)/peel recursion lift to T-equivariant characters?
# ch_k(j) = Demazure character D_w e^{j omega_m} = chi_T(X_k, O(j)) (Demazure formula, valid for all integers j).
# Solve ch_k(j+1) = sum_{rho ⊆ k} c_{k rho}(x) ch_rho(j) over Q(x) on rows j=0..|R|-1, verify on all rows
# j in [-4, |K|+2] (incl. negative j = strict side), test Laurent-ness, specialise x=1 and compare with
# S M^{-1} S and M = I+V (paper (★) / lem:peel), S = diag((-1)^{|k|}).
import sys, itertools
import sympy as sp
from sympy.polys.matrices import DomainMatrix
sys.path.insert(0,'../prove249d'); from dem import Di, word

def strict(m): return [tuple(sorted(c,reverse=True)) for L in range(m+1) for c in itertools.combinations(range(1,m+1),L)]
def sub(r,k): return len(r)<=len(k) and all(r[i]<=k[i] for i in range(len(r)))
def D(k): return [(a,b) for a in range(len(k)) for b in range(a,a+k[a])]
def V(theta):  # strictly order-preserving maps theta -> {e<o}, diagonal boxes -> e
    th=list(theta); cnt=0
    for vals in itertools.product((0,1),repeat=len(th)):
        f=dict(zip(th,vals))
        if any(f[x]!=0 for x in th if x[0]==x[1]): continue
        if all(not(x!=y and x[0]<=y[0] and x[1]<=y[1]) or f[x]<f[y] for x in th for y in th): cnt+=1
    return cnt

m=int(sys.argv[1]); X=sp.symbols(f'x1:{m+1}')
F=sp.QQ.frac_field(*X)
def ch(k,j):
    f={tuple([j]*m):1}
    for i in reversed(word(set(k),m)): f=Di(f,i,m)
    return F.convert(sp.Add(*[c*sp.Mul(*[X[t]**e[t] for t in range(m)]) for e,c in f.items()]))
K=strict(m); idx={k:i for i,k in enumerate(K)}
jr=list(range(-4,len(K)+3))
CH={(k,j):ch(k,j) for k in K for j in jr+[jr[-1]+1]}
Cx=[[F.zero]*len(K) for _ in K]; allok=True
lines=[]
for k in K:
    R=[r for r in K if sub(r,k)]
    A=DomainMatrix([[CH[(r,j)] for r in R] for j in range(len(R))],(len(R),len(R)),F)
    b=DomainMatrix([[CH[(k,j+1)]] for j in range(len(R))],(len(R),1),F)
    sol=[s[0] for s in A.lu_solve(b).to_Matrix().tolist()] if False else None
    solM=A.lu_solve(b); sol=[solM[i,0].element for i in range(len(R))]
    ok=all(sum((s*CH[(r,j)] for s,r in zip(sol,R)),F.zero)==CH[(k,j+1)] for j in jr)
    sols=[F.to_sympy(s) for s in sol]
    laur=all(sp.Poly(sp.denom(sp.together(s)),*X).is_monomial for s in sols)
    allok&=ok and laur
    lines.append(f'kappa={k}: j-independent solution verified on j in [{jr[0]},{jr[-1]}]: {ok}; all coefficients Laurent: {laur}')
    for s,r in zip(sols,R):
        Cx[idx[k]][idx[r]]=s
        if s!=0: lines.append(f'    c[{k},{r}] = {sp.factor(s)}')
Sg=sp.diag(*[(-1)**sum(k) for k in K])
M1=sp.Matrix(len(K),len(K),lambda a,c: (1 if a==c else (V(set(D(K[a]))-set(D(K[c]))) if sub(K[c],K[a]) else 0)))
C1=sp.Matrix(len(K),len(K),lambda a,c: sp.sympify(Cx[a][c]).subs({x:1 for x in X}) if Cx[a][c]!=0 else 0)
print(f'm={m}; strict partitions (row order): {K}')
print('\n'.join(lines))
print('C(1) =',C1.tolist())
print('I+V (lem:peel matrix) =',M1.tolist())
print('C(1) == S (I+V)^{-1} S :', C1==Sg*M1.inv()*Sg)
print('S C(1)^{-1} S == I+V  :', Sg*C1.inv()*Sg==M1)
print('ALL j-independent Laurent lifts:',allok)
# strict-side lift: n_k(j) = (-1)^{|k|} ch_k(-j) satisfies n(j+1) = P(x) n(j), P(x) = S C(x)^{-1} S
Cs=sp.Matrix(len(K),len(K),lambda a,c: Cx[a][c] if Cx[a][c]!=0 else 0)
Px=(Sg*Cs.inv()*Sg).applyfunc(lambda e: sp.factor(sp.cancel(e)))
print('Strict-side lift P(x) = S C(x)^{-1} S (rows kappa, cols rho, nonzero entries):')
pos=True; laurP=True
for a in range(len(K)):
    for c in range(len(K)):
        e=Px[a,c]
        if e!=0:
            num=sp.Poly(sp.numer(sp.together(e)),*X); den=sp.Poly(sp.denom(sp.together(e)),*X)
            laurP&=den.is_monomial; pos = pos and all(cf>0 for cf in num.coeffs()) and den.LC()>0
            print(f'    P[{K[a]},{K[c]}] = {e}   (at x=1: {e.subs({x:1 for x in X})}, V = {M1[a,c]})')
print('P(x) Laurent:',laurP,'; all coefficients positive:',pos,'; P(1)==I+V:',Px.subs({x:1 for x in X})==M1)
