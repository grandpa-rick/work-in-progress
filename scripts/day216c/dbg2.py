import sys; sys.argv=['x','3','1']
exec(open('nsym.py').read().split('res={}')[0])
def comm(bul,deg):
    Ms=[mat(lambda F,i=i:Y(i,F,bul),deg,deg+(1 if bul else 0)) for i in range(1,m+1)]
    return Ms
# bullet: compose Y_i Y_j from deg to deg+2
for deg in (0,1):
  for i in range(1,4):
    for j in range(i+1,4):
        F=[sp.prod([x**e for x,e in zip(X,c)]) for c in mons(deg)]
        print(deg,i,j,all(sp.expand(Y(i,Y(j,f,True),True)-Y(j,Y(i,f,True),True))==0 for f in F))
