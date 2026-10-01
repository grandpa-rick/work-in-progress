import sys; m=int(sys.argv[1]); d=int(sys.argv[2]); sys.argv=['x',str(m),str(d)]
exec(open('nsym_b.py').read().split('res={}')[0])
def fac(r): return None
def eig2(deg):
    Ys=[mat(lambda F,i=i:Y(i,F,False),deg,deg) for i in range(1,m+1)]
    lam=sp.Symbol('L'); vecs=[]
    cand=[sorted(sp.Poly(Yi.charpoly(lam).as_expr(),lam).ground_roots().keys()) for Yi in Ys]
    import itertools as it
    for ys in it.product(*cand):
        M=sp.Matrix.vstack(*[Ys[i]-ys[i]*sp.eye(Ys[0].rows) for i in range(m)])
        for v in M.nullspace(): vecs.append((v,list(ys)))
    assert len(vecs)==Ys[0].rows
    return vecs
for deg in range(d):
    E0=eig2(deg);E1=eig2(deg+1)
    P0=sp.Matrix.hstack(*[v for v,_ in E0]);P1=sp.Matrix.hstack(*[v for v,_ in E1])
    for name,op,dd in (('pi',lambda F:pi(F,False),0),('pibul',lambda F:pi(F,True),1),('X1',lambda F:sp.expand(X[0]*F),1),('Xm',lambda F:sp.expand(X[-1]*F),1)):
        Mt=(P1 if dd else P0).inv()*mat(op,deg,deg+dd)*P0
        nnz=[sum(1 for a in range(Mt.rows) if Mt[a,b]!=0) for b in range(Mt.cols)]
        print(deg,name,'nonzeros per column',nnz)
    # spectral check for X1: support mu has y(mu)= y(lam) with y1 scaled by s
    Mt=P1.inv()*mat(lambda F:sp.expand(X[0]*F),deg,deg+1)*P0
    good=True
    for b,(v,yl) in enumerate(E0):
        for a,(w,ym) in enumerate(E1):
            if Mt[a,b]!=0:
                exp=[s*yl[0]]+yl[1:]
                if sorted(ym)!=sorted(exp): good=False
    print(deg,'X1 spectral support claim:',good)
