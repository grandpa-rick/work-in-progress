from braces import *
import itertools
A=AbGroup((8,)); auts=A.automorphisms()
cnt=0
for lam in itertools.product(auts,repeat=8):
    if lam[0]!=tuple(range(8)): continue
    ok=all(lam[A.add[a][lam[a][b]]]==tuple(lam[a][lam[b][x]] for x in range(8)) for a in range(8) for b in range(8))
    cnt+=ok
print('brute Z/8 lambda-maps satisfying lam_{a o b}=lam_a lam_b:',cnt, 'search:',len(regular_subgroups(A)))
