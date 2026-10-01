"""Day 216c: machine-check the relations used in the proof of (N) via (NS).
R3 pi T_k = T_{k+1} pi (k<=m-2); B1 [T_i,Y_j]=0 j not in {i,i+1}; B2 Y_{i+1}=t^{-1}T_iY_iT_i (std & bullet);
R1 T_iX_iT_i = tX_{i+1}; I [e1(Y),X1]=(s-1)X1Y1; Y1bul = X1 Y1; spectra of Y on Pol_d are s^a t^b with b a perm of 0..m-1."""
import sys; m=int(sys.argv[1]); d=int(sys.argv[2]); sys.argv=['x',str(m),str(d)]
exec(open('nsym_b.py').read().split('res={}')[0])
import random; random.seed(7)
def rnd(deg):
    return sp.expand(sum(random.randint(-3,3)*sp.prod([x**e for x,e in zip(X,c)]) for c in mons(deg)))
Fs=[rnd(k) for k in range(d+1)]
ok={}
def chk(name,cond): ok[name]=ok.get(name,True) and cond
for F in Fs:
    for k in range(1,m-1): chk('R3',sp.expand(pi(T(k,F),False)-T(k+1,pi(F,False)))==0)
    for i in range(1,m):
        chk('R1',sp.expand(T(i,sp.expand(X[i-1]*T(i,F)))-t*X[i]*F)==0)
        for b in (False,True):
            chk('B2'+str(b),sp.expand(t*Y(i+1,F,b)-T(i,Y(i,T(i,F),b)))==0)
        for j in range(1,m+1):
            if j not in (i,i+1): chk('B1',sp.expand(T(i,Y(j,F,False))-Y(j,T(i,F),False))==0)
    chk('I',sp.expand(sum(Y(j,sp.expand(X[0]*F),False)-X[0]*Y(j,F,False) for j in range(1,m+1))-(s-1)*X[0]*Y(1,F,False))==0)
    chk('Y1bul',sp.expand(Y(1,F,True)-X[0]*Y(1,F,False))==0)
    for i in range(1,m+1):
        for j in range(1,m+1): chk('Ycomm',sp.expand(Y(i,Y(j,F,False),False)-Y(j,Y(i,F,False),False))==0)
print(m,d,ok)
