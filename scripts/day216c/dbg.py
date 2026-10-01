import sys; sys.argv=['x','3','1']
exec(open('nsym.py').read().split('res={}')[0])
lam=sp.Symbol('L')
for i in range(1,m+1):
    Yi=mat(lambda F,i=i:Y(i,F,False),1,1); print(i, sp.factor(Yi.charpoly(lam).as_expr()))
for i in range(1,m+1):
    for j in range(1,m+1):
        A=mat(lambda F,i=i:Y(i,F,False),1,1);B=mat(lambda F,j=j:Y(j,F,False),1,1); print(i,j,(A*B-B*A).is_zero_matrix)
