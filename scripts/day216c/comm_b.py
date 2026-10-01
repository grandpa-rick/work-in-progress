import sys; m=int(sys.argv[1]); sys.argv=['x',str(m),'1']
exec(open('nsym_b.py').read().split('res={}')[0])
import random
random.seed(1)
Fs=[sp.prod([x**random.randint(0,2) for x in X]) + 3*X[0]**2*X[-1] - X[1] for _ in range(3)]
for i in range(1,m+1):
  for j in range(1,m+1):
    ok=[]
    for F in Fs:
        c=sp.expand(Y(j,sp.expand(X[i-1]*F),False)-X[i-1]*Y(j,F,False))
        ok.append(c)
    # compare with multiples of X_i Y_i
    print(i,j,'[Y_j,X_i]==0:',all(c==0 for c in ok))
  for F in Fs:
    c=sp.expand(sum(Y(j,sp.expand(X[i-1]*F),False)-X[i-1]*Y(j,F,False) for j in range(1,m+1)))
    b=sp.expand(X[i-1]*Y(i,F,False))
    print(i,'ratio [e1Y,X_i]/(X_iY_i):',sp.cancel(c/b))
# X relation
for i in range(1,m):
    F=Fs[0]
    print('Xrel',i,sp.expand(T(i,sp.expand(X[i-1]*T(i,F)))-t*X[i]*F)==0, sp.expand(T(i,sp.expand(X[i]*T(i,F)))-t*X[i-1]*F)==0)
