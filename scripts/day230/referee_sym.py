exec(open('referee_v2.py').read().split('ok=bad=0')[0])
for (a,b,c,x,y) in [(3,1,2,4,2),(2,2,3,5,2),(3,2,2,6,1),(4,1,2,4,3)]:
    print((a,b,c,x,y), sp.simplify(Phi(a,b,c,x,y)-Phi(a,b,c,y,x))==0)
