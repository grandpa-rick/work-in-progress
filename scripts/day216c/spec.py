import sys; m=int(sys.argv[1]); d=int(sys.argv[2]); sys.argv=['x',str(m),str(d)]
exec(open('nsym.py').read().split('res={}')[0])
allok=True
for deg in range(d+1):
    E=eig(deg); seen=set()
    for v,ys in E:
        ab=[fac(y) for y in ys]; bs=sorted(b for a,b in ab)
        allok &= bs==list(range(m)); seen.add(tuple(a for a,b in ab))
    allok &= len(seen)==len(E)  # s-exponent vector = lambda determines eigenvector (simple spectrum)
    allok &= seen==set(mons(deg))
print(m,d,'spectral form + simple spectrum + a=lambda:',allok)
