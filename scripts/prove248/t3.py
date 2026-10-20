from psi import *
for kap,M in [((2,),(2,3)),((2,),(3,4)),((3,),(3,4)),((2,3),(2,3,4)),((2,4),(2,3,4)),((1,3),(1,2,3)),((2,3),(3,4,5)),((1,),(2,3))]:
    p=Psi(kap,(0,)*len(kap),M); print(kap,M,p, 'free:',sorted(p.free_symbols,key=str))
